from abc import ABCMeta
from dataclasses import dataclass
from datetime import date
import logging
import os
from typing import List


from genie import create_case_lists, process_functions
import numpy as np
import pandas as pd
import synapseclient
from synapseclient import File

from geniesp import metafiles
from geniesp.extract import Extract, get_synid_data
from geniesp.config import BpcConfig
from geniesp.utils import (
    _convert_to_int,
    fill_cancer_dx_start_date,
    create_regimens,
    get_drug_mapping,
    remap_os_values,
    remap_pfs_values,
    change_days_to_years,
    hack_remap_laterality
)

# All cbioportal file formats written in BPC
CBIO_FILEFORMATS_ALL = [
    "data_timeline_treatment.txt",
    "data_timeline_cancer_diagnosis.txt",
    "data_timeline_pathology.txt",
    "data_timeline_sample_acquisition.txt",
    "data_timeline_performance_status.txt",
    "data_timeline_medonc.txt",
    "data_timeline_imaging.txt",
    "data_timeline_sequencing.txt",
    "data_timeline_labtest.txt",
    "data_clinical_supp_survival.txt",
    "data_clinical_supp_survival_treatment.txt",
    "data_clinical_sample.txt",
    "data_clinical_patient.txt",
    "data_mutations_extended.txt",
    "data_gene_matrix.txt",
    "data_cna_hg19.seg",
    "data_fusions.txt",
    "data_sv.txt",
    "data_CNA.txt",
]


@dataclass
class Transforms(metaclass=ABCMeta):
    """Timeline data class."""
    # timeline_infodf: pd.DataFrame
    extract: Extract
    bpc_config: BpcConfig
    filepath: str = None

    def map_to_cbioportal_format(
        self
    ) -> dict:
        """Extracts the sample, patient and timeline data frame

        Args:
            syn (Synapse): Synapse connection
            mappingdf (pd.DataFrame): Mapping dataframe
            sampletype (str): sample type label
            cohort (str, optional): cohort label. Defaults to "NSCLC".

        Returns:
            dict: dictionary with two keys ('df' and 'used') corresponding to data frame
            of data for sample type and a list of Synapse IDs used
        """
        print(self.extract.timeline_infodf)
        mappingdf = self.extract.timeline_infodf[self.extract.timeline_infodf["data_type"] != "portal_value"]
        # Group by dataset because different datasets could have the
        # same variable
        datasets = mappingdf.groupby("dataset")
        finaldf = pd.DataFrame()
        dataset_map = self.extract.get_derived_variable_files()
        sample_type = self.extract.sample_type
        for dataset, df in datasets:
            # obtain columns to subset df
            cols = df["code"][df["sampleType"] == sample_type]
            cols = cols.tolist()
            if "record_id" not in cols:
                cols.append("record_id")
            # Must add path_rep_number for sample and sample acquisition file
            if sample_type in ["SAMPLE", "TIMELINE-SAMPLE"]:
                cols.append("path_rep_number")
            # Must add path_proc_number to sample file
            if sample_type == "SAMPLE":
                cols.append("path_proc_number")
            # Only get specific cohort and subset cols
            table = dataset_map['derived_variable_entities'][dataset]
            tabledf = pd.read_csv(table.path, low_memory=False)
            tabledf = tabledf[tabledf["cohort"] == self.bpc_config.cohort]
            tabledf = tabledf[cols]

            # Append to final dataframe if empty
            if finaldf.empty:
                finaldf = pd.concat([finaldf, tabledf])
            else:
                # Records missing pathology reports still have to be present
                # So a left merge has to happen.  This logic also assumes that
                # The pathology report dataset mapping info isn't first.
                # This also assumes that the TIMELINE-PATHOLOGY file only
                # uses columns from the pathology-report dataset
                if df["dataset"].iloc[0] == "Pathology-report level dataset":
                    finaldf = finaldf.merge(
                        tabledf,
                        on=["record_id", "path_proc_number", "path_rep_number"],
                        how="left",
                    )
                    del finaldf["path_rep_number"]
                else:
                    finaldf = finaldf.merge(tabledf, on="record_id")

        return finaldf

    def transforms(self, timelinedf, filter_start):
        # Obtain portal value (EVENT_TYPE)
        subset_infodf = self.extract.timeline_infodf
        print(subset_infodf)
        portal_value_idx = subset_infodf["data_type"] == "portal_value"
        # HACK: Passing in data mapping without portal_values should
        # still generate a file.
        portal_value = subset_infodf["code"][portal_value_idx].values[0]
        subset_infodf = subset_infodf[subset_infodf["data_type"] != "portal_value"]

        timelinedf["EVENT_TYPE"] = portal_value
        mapping = subset_infodf["cbio"].to_dict()
        # Must add in PATIENT_ID
        mapping["record_id"] = "PATIENT_ID"
        timelinedf = timelinedf.rename(columns=mapping)
        timelinedf["STOP_DATE"] = ""
        # timeline file must be in this order
        cols_to_order = ["PATIENT_ID", "START_DATE", "STOP_DATE", "EVENT_TYPE"]
        cols_to_order.extend(timelinedf.columns.drop(cols_to_order).tolist())
        timelinedf = self.retract_samples_and_patients(timelinedf)
        # Remove all null START_DATE rows if requested
        if filter_start:
            timelinedf = timelinedf[~timelinedf["START_DATE"].isnull()]
        return timelinedf[cols_to_order].drop_duplicates()

    def custom_transform(self, timelinedf):
        return timelinedf

    def create_timeline_file(
        self,
        filter_start: bool = True,
    ) -> dict:
        """Create timeline files straight from derived variables.

        Args:
            timeline_infodf (pd.DataFrame): cBio mapping information relevant to the timeline
            timeline_type (str): timeline label
            filter_start (bool, optional): whether to filter out rows with null START_DATEs. Defaults to True.

        Returns:
            dict: 'df' mapped dataframe, 'used' list of entities
        """
        timelinedf = self.map_to_cbioportal_format()
        timelinedf = self.transforms(timelinedf=timelinedf, filter_start=filter_start)
        timelinedf = self.custom_transform(timelinedf=timelinedf)
        return timelinedf

    def retract_samples_and_patients(self, df: pd.DataFrame) -> pd.DataFrame:
        """Make sure samples and patients exist in the main genie
        clinical samples with the exception of the removal of
        retraction database samples

        Args:
            df (pd.DataFrame): GENIE data with a patient or sample ID column

        Returns:
            pd.DataFrame: GENIE data with retracted patients or samples removed.
        """

        if df.get("SAMPLE_ID") is not None:
            to_keep_samples_idx = df["SAMPLE_ID"].isin(
                self.extract.to_keep_samples_and_patients["SAMPLE_ID"]
            )
            df = df[to_keep_samples_idx]
        elif df.get("PATIENT_ID") is not None:
            to_keep_patient_idx = df["PATIENT_ID"].isin(
                self.extract.to_keep_samples_and_patients["PATIENT_ID"]
            )
            df = df[to_keep_patient_idx]
        return df

    def write(
        self, df: pd.DataFrame
    ):
        """Write the file

        Args:
            df (pd.DataFrame): data frame to store
            filepath (str): path to written file
            used_entities (list, optional): Synapse IDs used to generate the file. Defaults to [].
        """

        df_text = process_functions.removePandasDfFloat(df)
        with open(self.filepath, "w") as file_f:
            file_f.write(df_text)

@dataclass
class TimelinePerformanceTransform(Transforms):
    """TimelinePerformance data class."""
    filepath = "data_timeline_performance_status.txt"
    def custom_transform(
        self, timelinedf
    ) -> dict:
        """Get TIMELINE-PERFORMANCE file data.

        Args:
            df_map (pd.DataFrame): variable to cBioPortal mapping info
            df_file (pd.DataFrame): data file to Synapse ID mapping

        Returns:
            dict: TIMELINE-PERFORMANCE data
        """
        # HACK: Due to remapping logic, we will re-create RESULT column with correct
        has_md_karnof = ~timelinedf['MD_KARNOF'].fillna('Not').str.startswith(("Not" ,"not"))
        has_md_ecog = ~timelinedf['MD_ECOG'].fillna('Not').str.startswith(("Not" ,"not"))
        # Only add in values for SCORE_TYPE and RESULT when MD_KARNOF
        # and ECOG are present
        timelinedf['SCORE_TYPE'] = ""
        timelinedf['SCORE_TYPE'][has_md_karnof] = "KARNOFSKY"
        timelinedf['SCORE_TYPE'][has_md_ecog] = "ECOG"
        timelinedf['RESULT'] = ""
        timelinedf['RESULT'][has_md_karnof] = timelinedf['MD_KARNOF'][has_md_karnof]
        timelinedf['RESULT'][has_md_ecog] = timelinedf['MD_ECOG'][has_md_ecog]
        # CbioPortal doesn't want any rows without MD_KARNOF or MD_ECOG
        timelinedf = timelinedf[timelinedf['RESULT'] != ""]
        timelinedf['RESULT'] = [
            _convert_to_int(val.split(":")[0])
            for val in timelinedf['RESULT']
        ]
        return timelinedf


@dataclass
class TimelineTreatmentTransform(Transforms):
    """TimelinePerformance data class."""

    def map_to_cbioportal_format(
        self
    ) -> dict:
        subset_infodf = self.extract.timeline_infodf[self.extract.timeline_infodf["sampleType"] == self.extract.sample_type]
        # Exclude heme onc columns
        subset_infodf = subset_infodf[
            ~subset_infodf["data_type"].isin(["portal_value", "heme"])
        ]
        dataset = subset_infodf["dataset"].unique()[0]
        dataset_map = self.extract.get_derived_variable_files()
        timelinedf = pd.read_csv(dataset_map['derived_variable_entities'][dataset].path, low_memory=False)
        # Only take lung cohort
        timelinedf = timelinedf[timelinedf["cohort"] == self.bpc_config.cohort]
        # Only take samples where redcap_ca_index is Yes
        timelinedf = timelinedf[timelinedf["redcap_ca_index"] == "Yes"]
        # Flatten multiple columns values into multiple rows
        multiple_cols_idx = subset_infodf["code"].str.contains("[*]")
        final_timelinedf = pd.DataFrame()
        for _, row in subset_infodf[multiple_cols_idx].iterrows():
            code = row["code"]
            # Make sure to use regex to get colname + integer
            new_code = code.replace("*", "[\d]")
            cols = timelinedf.columns.str.contains(new_code)
            wanted_cols = timelinedf.columns[cols].tolist()
            wanted_cols.extend(["record_id", "regimen_drugs", "regimen_number"])
            # melt function creates multiple rows from multiple columns
            melted_df = pd.melt(
                timelinedf[wanted_cols],
                id_vars=["record_id", "regimen_drugs", "regimen_number"],
                value_name=row["cbio"],
            )
            del melted_df["variable"]
            if final_timelinedf.empty:
                final_timelinedf = pd.concat([final_timelinedf, melted_df])
            else:
                final_timelinedf[row["cbio"]] = melted_df[row["cbio"]]
        final_timelinedf["TREATMENT_TYPE"] = "Systemic Therapy"
        # Remove all START_DATE is NULL
        final_timelinedf = final_timelinedf[~final_timelinedf["START_DATE"].isnull()]
        non_multi_cols = subset_infodf[~multiple_cols_idx]["code"].tolist()
        non_multi_cols.append("record_id")

        # Merge final timeline
        final_timelinedf = final_timelinedf.merge(
            timelinedf[non_multi_cols],
            on=["record_id", "regimen_drugs", "regimen_number"],
        )

        # Make sure all events types are treatment
        final_timelinedf["EVENT_TYPE"] = "TREATMENT"

        # Make sure AGENT is not null and doesn't have parenthesis
        agents = []
        final_timelinedf = final_timelinedf[~final_timelinedf["AGENT"].isnull()]
        for index, agent in enumerate(final_timelinedf["AGENT"]):
            if "(" in agent:
                agents.append(agent.split("(")[0].strip())
            else:
                agents.append(agent.split(",")[0].strip())

        final_timelinedf["AGENT"] = agents
        # Map timeline treatment columns
        mapping = subset_infodf["cbio"].to_dict()
        # Must add in PATIENT_ID
        mapping["record_id"] = "PATIENT_ID"
        final_timelinedf = final_timelinedf.rename(columns=mapping)
        # timeline file must be in this order
        cols_to_order = [
            "PATIENT_ID",
            "START_DATE",
            "STOP_DATE",
            "EVENT_TYPE",
            "TREATMENT_TYPE",
            "AGENT",
        ]
        cols_to_order.extend(final_timelinedf.columns.drop(cols_to_order).tolist())
        final_timelinedf = self.retract_samples_and_patients(final_timelinedf)

        return final_timelinedf[cols_to_order].drop_duplicates()

    def transforms(self, timelinedf, filter_start):
        return timelinedf


@dataclass
class TimelineTreatmentRadTransform(Transforms):
    """TimelinePerformance data class."""

    def custom_transform(
        self, timelinedf
    ) -> dict:
        """Get TIMELINE-PERFORMANCE file data.

        Args:
            df_map (pd.DataFrame): variable to cBioPortal mapping info
            df_file (pd.DataFrame): data file to Synapse ID mapping

        Returns:
            dict: TIMELINE-PERFORMANCE data
        """
        rad_df = timelinedf
        rad_df["STOP_DATE"] = rad_df["START_DATE"] + rad_df["TEMP"]
        rad_df = rad_df[rad_df["INDEX_CANCER"] == "Yes"]
        rad_df["EVENT_TYPE"] = "TREATMENT"
        rad_df["TREATMENT_TYPE"] = "Radiation Therapy"
        del rad_df["INDEX_CANCER"]
        del rad_df["TEMP"]
        return rad_df


class TimelineDxTransform(Transforms):

    def custom_transform(
        self, timelinedf
    ) -> pd.DataFrame:
        timelinedf = fill_cancer_dx_start_date(timelinedf)
        timelinedf = timelinedf[
            ~timelinedf["START_DATE"].isnull()
        ]
        return timelinedf


class TimelineTransform(Transforms):
    pass


class TimelineSampleTransform(Transforms):

    def custom_transform(
        self, timelinedf
    ) -> pd.DataFrame:
        # TODO: Can add getting of samples with NULL start dates in
        # self.create_fixed_timeline_files
        null_dates_idx = timelinedf["START_DATE"].isnull()
        if null_dates_idx.any():
            logging.warning(
                "timeline sample with null START_DATE: {}".format(
                    ", ".join(timelinedf["SAMPLE_ID"][null_dates_idx])
                )
            )
            timelinedf = timelinedf[~null_dates_idx]
        return timelinedf


class TimelineSequenceTransform(Transforms):

    def custom_transform(
        self, timelinedf
    ) -> pd.DataFrame:
        # HACK: Manually calculate the START_DATE based on criteria defined
        # in GEN-94
        # sequence_data["df"]["START_DATE"]
        seq_df =timelinedf
        index_seq_df = seq_df[seq_df["INDEX_CANCER"] == "Yes"]
        index_seq_df["START_DATE"] = (
            index_seq_df["DPT_REPORT_DAYS"] - index_seq_df["CA_DX_DAYS"]
        )
        del seq_df["START_DATE"]
        seq_df = seq_df.merge(index_seq_df[["SAMPLE_ID", "START_DATE"]], on="SAMPLE_ID")
        # Delete unneeded columns
        del seq_df["DPT_REPORT_DAYS"]
        del seq_df["INDEX_CANCER"]
        del seq_df["CA_DX_DAYS"]
        # reorder Columns to match requirement
        cols_to_order = ["PATIENT_ID", "START_DATE", "STOP_DATE", "EVENT_TYPE"]
        cols_to_order.extend(seq_df.columns.drop(cols_to_order).tolist())
        seq_df = seq_df[cols_to_order]
        seq_df.drop_duplicates(inplace=True)
        return seq_df


# TODO: Create a ClinicalTransform class
class SurvivalTransform(Transforms):

    def configure_clinicaldf(
        self, clinicaldf: pd.DataFrame, redcap_to_cbiomappingdf: pd.DataFrame
    ) -> pd.DataFrame:
        """Create clinical file from sponsored project mapped dataframe

        Args:
            clinicaldf (pd.DataFrame): clinical information
            redcap_to_cbiomappingdf (pd.DataFrame): cBio mapping info

        Raises:
            ValueError: All column names must be in mapping dataframe
            ValueError: Must have no null patient ids
            ValueError: Must have no null sample ids

        Returns:
            pd.DataFrame: configured clinical information
        """
        if not clinicaldf.columns.isin(redcap_to_cbiomappingdf["code"]).all():
            raise ValueError("All column names must be in mapping dataframe")
        mapping = redcap_to_cbiomappingdf["cbio"].to_dict()
        print(mapping)
        clinicaldf = clinicaldf.rename(columns=mapping)
        clinicaldf = clinicaldf.drop_duplicates()
        if sum(clinicaldf["PATIENT_ID"].isnull()) > 0:
            raise ValueError("Must have no null patient ids")
        # Remove white spaces for PATIENT/SAMPLE ID
        clinicaldf["PATIENT_ID"] = [
            patient.strip() for patient in clinicaldf["PATIENT_ID"]
        ]

        if clinicaldf.get("SAMPLE_ID") is not None:
            # This line should not be here
            clinicaldf = clinicaldf[~clinicaldf["SAMPLE_ID"].isnull()]
            if sum(clinicaldf["SAMPLE_ID"].isnull()) > 0:
                raise ValueError("Must have no null sample ids")
            clinicaldf["SAMPLE_ID"] = [
                sample.strip() for sample in clinicaldf["SAMPLE_ID"]
            ]
        # JIRA: GEN-10- Must standardize SEQ_ASSAY_ID values to be uppercase
        if clinicaldf.get("SEQ_ASSAY_ID") is not None:
            clinicaldf["SEQ_ASSAY_ID"] = clinicaldf["SEQ_ASSAY_ID"].str.upper()

        clinicaldf["SP"] = self.bpc_config.cohort

        for col in clinicaldf:
            num_missing = sum(clinicaldf[col].isnull())
            if num_missing > 0:
                logging.warning(f"Number of missing {col}: {num_missing}")

        return clinicaldf

    def transforms(self, timelinedf, filter_start) -> dict:

        df_final_survival = self.configure_clinicaldf(timelinedf, self.extract.timeline_infodf)

        # Only take rows where cancer index is null
        df_final_survival = df_final_survival[
            df_final_survival["CANCER_INDEX"].isnull()
        ]
        # Remove cancer index column
        del df_final_survival["CANCER_INDEX"]
        # remove a row if patient ID is duplicated and PFS_I_ADV_STATUS is null or empty
        # tested on current survival data file and produces unique patient list
        if "PFS_I_ADV_STATUS" in df_final_survival.columns:
            pfs_not_null_idx = ~df_final_survival["PFS_I_ADV_STATUS"].isnull()
            pfs_not_blank_idx = df_final_survival["PFS_I_ADV_STATUS"] != ""
            nondup_patients_idx = ~df_final_survival["PATIENT_ID"].duplicated(
                keep=False
            )
            df_final_survival = df_final_survival[
                (pfs_not_null_idx & pfs_not_blank_idx) | (nondup_patients_idx)
            ]
        # Only patients and samples that exist in the
        # sponsored project uploads are going to be pulled into the SP project
        subset_survivaldf = self.retract_samples_and_patients(df_final_survival)

        del subset_survivaldf["SP"]
        subset_survivaldf = remap_os_values(df=subset_survivaldf)
        subset_survivaldf = remap_pfs_values(df=subset_survivaldf)
        cols_to_order = ["PATIENT_ID"]
        cols_to_order.extend(subset_survivaldf.columns.drop(cols_to_order).tolist())

        return subset_survivaldf[cols_to_order]

    def write(
        self,
        df: pd.DataFrame,
    ) -> str:
        """Writes out the clinical file

        Args:
            clinicaldf (pd.DataFrame): clinical information
            redcap_to_cbiomappingdf (pd.DataFrame): cBio mapping info
            filetype (str): file type label

        Raises:
            ValueError: sample type must be patient, sample, supp_survival or supp_survival_treatment

        Returns:
            str: file path to clinical info
        """
        # TODO fix this
        # if filetype not in [
        #     "patient",
        #     "sample",
        #     "supp_survival",
        #     "supp_survival_treatment",
        # ]:
        #     raise ValueError(
        #         "sample type must be patient, sample, supp_survival or "
        #         "supp_survival_treatment"
        #     )
        # Must have this for the dict mappings after
        redcap_to_cbiomappingdf = self.extract.survival_info_df
        redcap_to_cbiomappingdf.index = redcap_to_cbiomappingdf["cbio"]
        # HACK this is to remove the index cancer for the patient file
        # redcap_to_cbiomappingdf = redcap_to_cbiomappingdf[redcap_to_cbiomappingdf['cbio'] != "INDEX_CANCER"]
        label_map = redcap_to_cbiomappingdf["labels"].to_dict()
        description_map = redcap_to_cbiomappingdf["description"].to_dict()
        coltype_map = redcap_to_cbiomappingdf["colType"].to_dict()
        # Priority column will determine which columns are shown on cBioPortal
        # Must columns should be shown on cBioPortal will have a 1
        # but survival columns should all be 0
        priority_map = redcap_to_cbiomappingdf["priority"].to_dict()
        labels = [str(label_map[col]) for col in df]
        descriptions = [str(description_map[col]) for col in df]
        coltype = [str(coltype_map[col]) for col in df]
        priority = [str(int(priority_map[col])) for col in df]

        # clin_path = os.path.join(
        #     self._SPONSORED_PROJECT, f"data_clinical_{filetype}.txt"
        # )

        with open(self.filepath, "w+") as clin_file:
            clin_file.write("#{}\n".format("\t".join(labels)))
            clin_file.write("#{}\n".format("\t".join(descriptions)))
            clin_file.write("#{}\n".format("\t".join(coltype)))
            # TODO attributes in the supp file are PATIENT, so must
            # specify that.
            if "SURVIVAL" in self.filepath or "REGIMEN" in self.filepath:
                clin_file.write("#{}\n".format("\t".join(["PATIENT"] * len(labels))))
            clin_file.write("#{}\n".format("\t".join(priority)))
            clin_file.write(
                process_functions.removeStringFloat(
                    df.to_csv(index=False, sep="\t")
                )
            )
        return self.filepath

class SurvivalTreatmentTransform(SurvivalTransform):

    def map_to_cbioportal_format(self):
        # TODO: Make sure these are in the extract class
        drug_mapping = get_drug_mapping(
            syn=self.extract.syn,
            cohort=self.bpc_config.cohort,
            synid_table_prissmm=self.bpc_config.prissmm_synid,
        )
        regimens_data = create_regimens(
            self.extract.syn,
            self.extract.timeline_infodf,
            mapping=drug_mapping,
            top_x_regimens=20,
            cohort=self.bpc_config.cohort,
        )

        df_survival_treatment = regimens_data["df"]
        df_survival_treatment = remap_os_values(df=df_survival_treatment)
        df_survival_treatment = remap_pfs_values(df=df_survival_treatment)
        cols_to_order = ["PATIENT_ID"]
        cols_to_order.extend(df_survival_treatment.columns.drop(cols_to_order).tolist())
        # Retract patients from survival treatment file
        df_survival_treatment = self.retract_samples_and_patients(df_survival_treatment)

        return df_survival_treatment[cols_to_order]

    def transforms(self, timelinedf, filter_start):
        return timelinedf


class SampleTransform(SurvivalTransform):

    def transforms(self, timelinedf, filter_start) -> dict:
        del timelinedf["path_proc_number"]
        df_sample = self.configure_clinicaldf(timelinedf, self.extract.timeline_infodf)
        df_sample_subset = self.retract_samples_and_patients(df_sample)

        del df_sample_subset["SP"]
        days_to_years_col = [
            "AGE_AT_SEQ_REPORT_YEARS",
            "CPT_ORDER_INT",
            "CPT_REPORT_INT",
        ]
        for col in days_to_years_col:
            # not all columns could exist, so check if column exists
            if col in df_sample_subset:
                years = df_sample_subset[col].apply(change_days_to_years)
                df_sample_subset[col] = years
        # Use np.floor because np handles NaN values
        df_sample_subset["AGE_AT_SEQUENCING"] = df_sample_subset[
            "AGE_AT_SEQUENCING"
        ].apply(np.floor)
        # Remove CPT_SEQ_DATE because the values are incorrect
        del df_sample_subset["CPT_SEQ_DATE"]
        # Obtain this information from the main GENIE cohort
        df_sample_subset = df_sample_subset.merge(
            self.extract.to_keep_samples_and_patients[["SAMPLE_ID", "SEQ_YEAR"]],
            on="SAMPLE_ID",
            how="left",
        )
        df_sample_subset.rename(columns={"SEQ_YEAR": "CPT_SEQ_DATE"}, inplace=True)
        df_sample_subset.sort_values("PDL1_POSITIVE_ANY", ascending=False, inplace=True)
        df_sample_subset.drop_duplicates("SAMPLE_ID", inplace=True)

        return df_sample_subset


class PatientTransform(SurvivalTransform):

    def transforms(self, timelinedf, filter_start) -> dict:
        df_patient = timelinedf[timelinedf["redcap_ca_index"] == "Yes"]
        df_patient.drop(columns="redcap_ca_index", inplace=True)
        df_patient_final = self.configure_clinicaldf(df_patient, self.extract.timeline_infodf)

        df_patient_subset = self.retract_samples_and_patients(df_patient_final)

        # Fix patient duplicated values due to cancer index DOB
        # Take the larger DX_LASTALIVE_INT_MOS value for all records
        # subset_patientdf.sort_values("DX_LASTALIVE_INT_MOS", inplace=True,
        #                              ascending=False)
        df_patient_subset.drop_duplicates("PATIENT_ID", inplace=True)
        duplicated = df_patient_subset.PATIENT_ID.duplicated()
        if duplicated.any():
            logging.warning(
                "DUPLICATED PATIENT_IDs: {}".format(
                    ",".join(df_patient_subset["PATIENT_ID"][duplicated])
                )
            )

        del df_patient_subset["SP"]
        cols_to_order = ["PATIENT_ID"]
        cols_to_order.extend(df_patient_subset.columns.drop(cols_to_order).tolist())

        df_patient_subset = hack_remap_laterality(df_patient_subset=df_patient_subset)

        return df_patient_subset[cols_to_order]


def configure_mafdf(mafdf: pd.DataFrame, keep_samples: list) -> pd.DataFrame:
    """Configures a maf dataframe

    Args:
        mafdf (pd.DataFrame): Chunk of maf dataframe
        keep_samples (list):  Samples to keep in the maf file

    Returns:
        pd.DataFrame: Configured maf dataframe
    """
    keep_mafdf = mafdf[mafdf["Tumor_Sample_Barcode"].isin(keep_samples.tolist())]
    if not keep_mafdf.empty:
        fillnas = [
            "t_depth",
            "t_ref_count",
            "t_alt_count",
            "n_depth",
            "n_ref_count",
            "n_alt_count",
        ]
        for col in fillnas:
            keep_mafdf[col].loc[keep_mafdf[col] == "."] = ""
        keep_mafdf["Validation_Status"] = ""
    return keep_mafdf


@dataclass
class MainGenie:
    """BPC redcap to cbioportal export"""
    bpc_config: BpcConfig
    extract: Extract
    sample_df: pd.DataFrame
    patient_df: pd.DataFrame
    release: str
    cbioportal_folders: dict
    syn: synapseclient.Synapse
    upload: bool = False

    def create_bpc_cbio_metafiles(self) -> List:
        """Create BPC cBioPortal meta* files.

        Returns:
            List: file paths of meta files.
        """

        mg_release_ent = self.syn.get(self.bpc_config.mg_release_synid)
        name = f"GENIE BPC {self.bpc_config.cohort} v{self.release}"
        description = (
            f"{self.bpc_config.cohort} cohort v{self.release} "
            f"(GENIE {date.today().year}) GENIE {mg_release_ent.name}. "
            f"Several hundred different variables are collected for each of "
            f'the BPC cohorts; consult the <a href="{self.bpc_config.url_bpc}">Documentation</a> '
            f"for further detail. To learn more about which variables are "
            f"visualized in cBioPortal and how, see the cBioPortal "
            f'<a href="{self.bpc_config.url_cbio}">ReadMe</a>. '
            '<font color="red">Although these data are de-identified, your analysis '
            "may require institutional review prior to publication.</font>"
        )
        short_name = f"{self.bpc_config.cohort} GENIE"
        study_identifier = f"{self.bpc_config.cohort.lower()}_genie_bpc"
        # Get list of files to create cBioPortal metadata files for
        to_create_meta = list(set(CBIO_FILEFORMATS_ALL) - set(self.bpc_config.exclude_files))

        meta_files = metafiles.create_cbio_metafiles(
            study_identifier=study_identifier,
            outdir=self.bpc_config.cohort,
            cbio_fileformats=to_create_meta,
        )
        meta_study = metafiles.create_study_meta_file(
            study_identifier=study_identifier,
            type_of_cancer="mixed",
            name=name,
            description=description,
            groups="GENIE",
            short_name=short_name,
        )
        study_file = metafiles.write_meta_file(
            meta_info=meta_study,
            filename="meta_study.txt",
            outdir=self.bpc_config.cohort,
        )
        meta_files.append(study_file)
        return meta_files

    def create_and_write_genematrix(
        self, clinicaldf: pd.DataFrame, cna_samples: list, used_ent: list = None
    ) -> None:
        """Create gene matrix dataframe.

        Args:
            clinicaldf (pd.DataFrame): sample information
            cna_samples (list): sample IDs with CNA data
            used_ent (list, optional): Synapse IDs used. Defaults to None.
        """

        data_gene_panel = clinicaldf[["SAMPLE_ID", "SEQ_ASSAY_ID"]]
        data_gene_panel = data_gene_panel.rename(columns={"SEQ_ASSAY_ID": "mutations"})
        data_gene_panel = data_gene_panel[data_gene_panel["SAMPLE_ID"] != ""]
        data_gene_panel.drop_duplicates("SAMPLE_ID", inplace=True)
        cna_seqids = data_gene_panel["mutations"][
            data_gene_panel["SAMPLE_ID"].isin(cna_samples)
        ].unique()
        data_gene_panel["cna"] = data_gene_panel["mutations"]
        data_gene_panel["cna"][~data_gene_panel["cna"].isin(cna_seqids)] = "NA"
        data_gene_panel.fillna("NA", inplace=True)
        gene_matrix_filepath = os.path.join(
            self.bpc_config.cohort, "data_gene_matrix.txt"
        )
        data_gene_panel.to_csv(gene_matrix_filepath, sep="\t", index=False)
        if self.upload:
            file_ent = File(
                gene_matrix_filepath, parent=self.cbioportal_folders["release"]
            )
            self.syn.store(file_ent, used=used_ent, executed=self.bpc_config.github_url)

    def write_and_storedf(
        self, df: pd.DataFrame, filepath: str, used_entities: list = []
    ):
        """Write and, if applicable, store data frame.

        Args:
            df (pd.DataFrame): data frame to store
            filepath (str): path to written file
            used_entities (list, optional): Synapse IDs used to generate the file. Defaults to [].
        """

        df_text = process_functions.removePandasDfFloat(df)
        with open(filepath, "w") as file_f:
            file_f.write(df_text)

        if self.upload:
            # Add the mapping file to the release file provenance
            used_entities.append(self.bpc_config.redcap_to_cbio_mapping_synid)
            ent = File(filepath, parent=self.cbioportal_folders["release"])
            self.syn.store(ent, executed=self.bpc_config.github_url, used=used_entities)

    def create_and_write_maf(self, keep_samples: list) -> str:
        """Create maf file from release maf

        Args:
            keep_samples (list): List of samples to keep

        Returns:
            str: file path to written data
        """
        file_name = "data_mutations_extended.txt"
        mafpath = os.path.join(self.bpc_config.cohort, file_name)
        maf_synid = self.extract.get_mg_synid(self.bpc_config.mg_release_synid, file_name)
        maf_ent = self.syn.get(maf_synid, followLink=True)
        maf_chunks = pd.read_table(maf_ent.path, chunksize=50000, low_memory=False)
        index = 0
        for maf_chunk in maf_chunks:
            mafdf = configure_mafdf(maf_chunk, keep_samples)
            # Skip to next chunk if empty
            if mafdf.empty:
                continue
            # If maf file has not been created
            if index == 0:
                maf_text = process_functions.removePandasDfFloat(mafdf)
                with open(mafpath, "w") as maf_f:
                    maf_f.write(maf_text)
            else:
                maf_text = mafdf.to_csv(sep="\t", header=None, index=False)
                maf_text = process_functions.removeStringFloat(maf_text)
                with open(mafpath, "a") as maf_f:
                    maf_f.write(maf_text)
            index += 1
        if self.upload:
            file_ent = File(mafpath, parent=self.cbioportal_folders["release"])
            self.syn.store(file_ent, used=[maf_synid], executed=self.bpc_config.github_url)

        return mafpath

    def create_and_write_cna(self, keep_samples: list) -> dict:
        """Create CNA file

        Args:
            keep_samples (list): List of samples to keep

        Returns:
            dict: "filepath" with path to written file
                    "cna_sample" list of CNA sample IDs
        """
        file_name = "data_CNA.txt"
        cna_synid = self.extract.get_mg_synid(self.bpc_config.mg_release_synid, file_name)
        cna_path = os.path.join(self.bpc_config.cohort, file_name)
        cna_ent = self.syn.get(cna_synid, followLink=True)
        cnadf = pd.read_table(cna_ent.path, low_memory=False)
        keep_cols = ["Hugo_Symbol"]
        keep_cols.extend(cnadf.columns[cnadf.columns.isin(keep_samples)].tolist())
        cnadf = cnadf[keep_cols]
        cna_text = process_functions.removePandasDfFloat(cnadf)
        # Must do this replace twice because \t\t\t ->
        # \tNA\t\t -> \tNA\tNA\t
        cna_text = (
            cna_text.replace("\t\t", "\tNA\t")
            .replace("\t\t", "\tNA\t")
            .replace("\t\n", "\tNA\n")
        )

        with open(cna_path, "w") as cna_file:
            cna_file.write(cna_text)

        if self.upload:
            file_ent = File(cna_path, parent=self.cbioportal_folders["release"])
            self.syn.store(file_ent, used=[cna_synid], executed=self.bpc_config.github_url)
        return {"filepath": cna_file, "cna_samples": cnadf.columns.tolist()}

    def create_and_write_fusion(self, keep_samples: list) -> str:
        """Create fusion file

        Args:
            keep_samples (list): List of samples to keep

        Returns:
            str: file path to written data
        """
        file_name = "data_fusions.txt"
        fusion_synid = self.extract.get_mg_synid(self.bpc_config.mg_release_synid, file_name)
        fusion_ent = self.syn.get(fusion_synid, followLink=True)
        fusiondf = pd.read_table(fusion_ent.path, low_memory=False)
        fusiondf = fusiondf[fusiondf["Tumor_Sample_Barcode"].isin(keep_samples)]
        # cBioPortal validation fails when Hugo Symbol is null
        fusiondf = fusiondf[~fusiondf["Hugo_Symbol"].isnull()]
        fusion_path = os.path.join(self.bpc_config.cohort, file_name)
        self.write_and_storedf(fusiondf, fusion_path, used_entities=[fusion_synid])
        return fusion_path

    def create_and_write_seg(self, keep_samples: list) -> str:
        """Create seg file

        Args:
            keep_samples (list): List of samples to keep

        Returns:
            str: file path to written data
        """
        # TODO: the seg filename will change 13.X release.
        file_name = "data_cna_hg19.seg"
        seg_synid = self.extract.get_mg_synid(self.bpc_config.mg_release_synid, file_name)
        seg_ent = self.syn.get(seg_synid, followLink=True)
        segdf = pd.read_table(seg_ent.path, low_memory=False)
        segdf = segdf[segdf["ID"].isin(keep_samples)]
        seg_path = os.path.join(self.bpc_config.cohort, "data_cna_hg19.seg")
        self.write_and_storedf(segdf, seg_path, used_entities=[seg_synid])
        return seg_path

    def create_and_write_sv(self, keep_samples):
        """Create sv file

        Args:
            keep_samples: List of samples to keep
        """
        file_name = "data_sv.txt"
        # TODO: remove try except after using main genie release >= 13
        try:
            sv_synid = self.extract.get_mg_synid(self.bpc_config.mg_release_synid, file_name)
        except ValueError:
            sv_synid = None
            logging.warning(
                f"data_sv.txt doesn't exist in main genie release: {self.bpc_config.mg_release_synid}"
            )
        if sv_synid is not None:
            sv_ent = self.syn.get(sv_synid, followLink=True)
            svdf = pd.read_table(sv_ent.path, low_memory=False)
            svdf = svdf[svdf["Sample_Id"].isin(keep_samples)]
            sv_path = os.path.join(self.bpc_config.cohort, "data_sv.txt")
            self.write_and_storedf(svdf, sv_path, used_entities=[sv_synid])

    def create_and_write_gene_panels(self, keep_seq_assay_ids: list) -> List:
        """Create gene panels.

        Args:
            keep_seq_assay_ids (list): list of sequence assay IDs

        Returns:
            list: file paths of written data
        """

        gene_panel_paths = []
        file_name = "genomic_information.txt"
        genomic_info_synid = self.extract.get_mg_synid(self.bpc_config.mg_release_synid, file_name)
        genomic_info_ent = self.syn.get(genomic_info_synid, followLink=True)
        genomic_infodf = pd.read_table(genomic_info_ent.path, low_memory=False)
        # Filter by SEQ_ASSAY_ID and only exonic regions
        genomic_infodf = genomic_infodf[
            (genomic_infodf["SEQ_ASSAY_ID"].isin(keep_seq_assay_ids))
            & (genomic_infodf["Feature_Type"] == "exon")
            & (~genomic_infodf["Hugo_Symbol"].isnull())
            & (genomic_infodf["includeInPanel"])
        ]

        seq_assay_groups = genomic_infodf.groupby("SEQ_ASSAY_ID")
        for seq_assay_id, seqdf in seq_assay_groups:
            unique_genes = seqdf.Hugo_Symbol.unique()
            gene_panel_text = (
                "stable_id: {seq_assay_id}\n"
                "description: {seq_assay_id}, "
                "Number of Genes - {num_genes}\n"
                "gene_list:\t{genelist}".format(
                    seq_assay_id=seq_assay_id,
                    num_genes=len(unique_genes),
                    genelist="\t".join(unique_genes),
                )
            )
            gene_panel_name = f"data_gene_panel_{seq_assay_id}.txt"
            gene_panel_path = os.path.join(self.bpc_config.cohort, gene_panel_name)
            gene_panel_paths.append(gene_panel_path)
            with open(gene_panel_path, "w+") as f:
                f.write(gene_panel_text)
            if self.upload:
                fileEnt = File(
                    gene_panel_path, parent=self.cbioportal_folders["release"]
                )
                self.syn.store(
                    fileEnt, used=[genomic_info_synid], executed=self.bpc_config.github_url
                )
        return gene_panel_paths

    def create_and_write_case_lists(
        self, subset_sampledf: pd.DataFrame, subset_patientdf: pd.DataFrame, used: list
    ) -> None:
        """Create, write, and, if applicable, store case list files.

        Args:
            subset_sampledf: Dataframe of sample information
            subset_patientdf: Dataframe of patient information
            used: Synapse IDs used to construct the case file data.
        """
        # Merged clinical file is created here because it is needed
        # for the case lists
        # Remove oncotree code here, because no longer need it
        merged_clinicaldf = subset_sampledf.merge(
            subset_patientdf, on="PATIENT_ID", how="outer"
        )
        missing_sample_idx = merged_clinicaldf["SAMPLE_ID"].isnull()
        # Make sure there are no missing sample ids
        if sum(missing_sample_idx) > 0:
            logging.info(
                "MISSING SAMPLE_ID for: {}".format(
                    ",".join(merged_clinicaldf["PATIENT_ID"][missing_sample_idx])
                )
            )
            merged_clinicaldf = merged_clinicaldf[~missing_sample_idx]

        # TODO Add this back in
        # # upload samples that are not part of the main GENIE cohort
        # if merged_clinicaldf.get("SAMPLE_ID") is not None:
        #     logging.info("Samples not in GENIE clinical databases (SP and normal)")
        #     not_found_samples = merged_clinicaldf["SAMPLE_ID"][
        #         ~merged_clinicaldf["SAMPLE_ID"].isin(self.genie_clinicaldf["SAMPLE_ID"])
        #     ]
        #     if not not_found_samples.empty:
        #         logging.info(not_found_samples[~not_found_samples.isnull()])
        #         not_found_samples.to_csv("notfoundsamples.csv")
        #         if self.upload:
        #             self.syn.store(
        #                 File(
        #                     "notfoundsamples.csv", parent=self.bpc_config.sp_redcap_exports_synid
        #                 )
        #             )
        # Hard coded most up to date oncotree version
        # oncotreelink = self.syn.get("syn13890902").externalURL
        # Use the old oncotree link for now
        # oncotreelink = (
        #     "http://oncotree.mskcc.org/api/tumorTypes/tree?version=oncotree_2018_06_01"
        # )
        oncotreelink = self.bpc_config.oncotreelink
        oncotree_dict = process_functions.get_oncotree_code_mappings(oncotreelink)
        # Map cancer type and cancer type detailed
        # This is to create case list files
        merged_clinicaldf["CANCER_TYPE"] = [
            oncotree_dict[code.upper()].get("CANCER_TYPE", float("nan"))
            for code in merged_clinicaldf["ONCOTREE_CODE"]
        ]
        merged_clinicaldf["CANCER_TYPE_DETAILED"] = [
            oncotree_dict[code.upper()].get("CANCER_TYPE_DETAILED", float("nan"))
            for code in merged_clinicaldf["ONCOTREE_CODE"]
        ]
        merged_clinicaldf["ONCOTREE_PRIMARY_NODE"] = [
            oncotree_dict[code.upper()].get("ONCOTREE_PRIMARY_NODE", float("nan"))
            for code in merged_clinicaldf["ONCOTREE_CODE"]
        ]
        merged_clinicaldf["ONCOTREE_SECONDARY_NODE"] = [
            oncotree_dict[code.upper()].get("ONCOTREE_SECONDARY_NODE", float("nan"))
            for code in merged_clinicaldf["ONCOTREE_CODE"]
        ]
        # Remove duplicated sample ids (there shouldn't be any)
        merged_clinicaldf = merged_clinicaldf.drop_duplicates("SAMPLE_ID")
        merged_clinicaldf.to_csv(
            os.path.join(self.bpc_config.cohort, "data_clinical.txt"),
            index=False,
            sep="\t",
        )

        # Create case lists
        case_list_path = os.path.join(self.bpc_config.cohort, "case_lists")

        if not os.path.exists(case_list_path):
            os.mkdir(case_list_path)
        else:
            caselists = os.listdir(case_list_path)
            for caselist in caselists:
                os.remove(os.path.join(case_list_path, caselist))

        # Write out cases sequenced so people can tell
        # which samples were sequenced
        assay_info = self.syn.tableQuery(
            f"select * from {self.bpc_config.mg_assay_synid}",
            includeRowIdAndRowVersion=False,
            separator="\t",
        )
        create_case_lists.main(
            os.path.join(self.bpc_config.cohort, "data_clinical.txt"),
            assay_info.filepath,
            case_list_path,
            f"{self.bpc_config.cohort.lower()}_genie_bpc",
        )

        case_list_files = os.listdir(case_list_path)
        for casepath in case_list_files:
            casepath = os.path.join(case_list_path, casepath)
            if self.upload:
                file_ent = File(casepath, parent=self.cbioportal_folders["case_lists"])
                self.syn.store(
                    file_ent,
                    used=used,
                    executed=self.bpc_config.github_url,
                )

    def run(self):
        """Runs the redcap export to export all files"""

        logging.info("writing CLINICAL-SAMPLE...")
        df_sample_final = self.sample_df
        df_patient_final = self.patient_df
        logging.info("writing genomic data files...")
        self.create_and_write_maf(df_sample_final["SAMPLE_ID"])
        dict_cna = self.create_and_write_cna(df_sample_final["SAMPLE_ID"])
        self.create_and_write_genematrix(df_sample_final, dict_cna["cna_samples"])
        # self.create_and_write_fusion(df_sample_final["SAMPLE_ID"])
        self.create_and_write_seg(df_sample_final["SAMPLE_ID"])
        self.create_and_write_sv(df_sample_final["SAMPLE_ID"])
        # TODO: probably add this in bpc_config or extract
        ids = get_synid_data(
            df_map=self.extract.mapping_df,
            df_file=self.extract.data_tables_df,
            sampletype=["PATIENT", "SAMPLE"],
            cohort=self.bpc_config.cohort,
        )
        self.create_and_write_case_lists(
            subset_sampledf=df_sample_final, subset_patientdf=df_patient_final, used=ids
        )
        self.create_and_write_gene_panels(df_sample_final["SEQ_ASSAY_ID"].unique())

        logging.info("writing metadata files...")
        metadata_files = self.create_bpc_cbio_metafiles()
        if self.upload:
            for metadata_file in metadata_files:
                file_ent = synapseclient.File(
                    metadata_file, parent=self.cbioportal_folders["release"]
                )
                self.syn.store(
                    file_ent,
                    executed=self.bpc_config.github_url
                )
