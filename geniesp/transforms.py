from abc import ABCMeta
from dataclasses import dataclass
import logging

import pandas as pd
import numpy as np

from geniesp.extract import Extract
from geniesp.config import BpcConfig
from geniesp.bpc_redcap_export_mapping import (
    _convert_to_int,
    fill_cancer_dx_start_date,
    create_regimens,
    get_drug_mapping,
    remap_os_values,
    remap_pfs_values,
    change_days_to_years,
    hack_remap_laterality
)


@dataclass
class Transforms(metaclass=ABCMeta):
    """Timeline data class."""
    # timeline_infodf: pd.DataFrame
    extract: Extract
    bpc_config: BpcConfig

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


@dataclass
class TimelinePerformanceTransform(Transforms):
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
        # TODO Fix this code chunk
        # Only patients and samples that exist in the
        # df_file = self.extract.data_tables_df
        # df_map = self.extract.timeline_infodf
        # # Patient and Sample mapping values
        # patient_sample_idx = df_map["sampleType"].isin(
        #     ["PATIENT", "SAMPLE", "SURVIVAL"]
        # )
        # infodf = df_map[patient_sample_idx].merge(df_file, on="dataset", how="left")
        # infodf.index = infodf["code"]

        # # Regimen mapping values
        # regimen_idx = df_map["sampleType"].isin(["REGIMEN"])
        # regimen_infodf = df_map[regimen_idx].merge(df_file, on="dataset", how="left")
        # regimen_infodf.index = regimen_infodf["code"]

        # # Create regimens data for patient file
        # drug_mapping = get_drug_mapping(
        #     syn=self.syn,
        #     cohort=self._SPONSORED_PROJECT,
        #     synid_table_prissmm=self._PRISSMM_SYNID,
        # )
        # regimens_data = create_regimens(
        #     self.syn,
        #     regimen_infodf,
        #     mapping=drug_mapping,
        #     top_x_regimens=20,
        #     cohort=self._SPONSORED_PROJECT,
        # )

        # survival_info = pd.concat([infodf, regimens_data["info"]])

        # TODO check the above
        # Only patients and samples that exist in the
        # sponsored project uploads are going to be pulled into the SP project
        subset_survivaldf = self.retract_samples_and_patients(df_final_survival)

        del subset_survivaldf["SP"]
        subset_survivaldf = remap_os_values(df=subset_survivaldf)
        subset_survivaldf = remap_pfs_values(df=subset_survivaldf)
        cols_to_order = ["PATIENT_ID"]
        cols_to_order.extend(subset_survivaldf.columns.drop(cols_to_order).tolist())

        return subset_survivaldf[cols_to_order]


class SurvivalTreatmentTransform(Transforms):

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
