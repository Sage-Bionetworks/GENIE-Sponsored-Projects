from abc import ABCMeta
from dataclasses import dataclass
import os

import pandas as pd
import synapseclient
from synapseclient import Synapse
from genie import process_functions

from extract import Extract
from config import BpcConfig


def _convert_to_int(value):
    """Convert object to integer or return nan"""
    try:
        return int(value)
    except ValueError:
        return float('nan')


def fill_cancer_dx_start_date(finaldf: pd.DataFrame) -> pd.DataFrame:
    """Fills in cancer dx start date for those missing start dates
    with zero.

    Args:
        finaldf (pd.DataFrame): Mapped cancer diagnosis timeline dataframe

    Returns:
        pd.DataFrame: dataframe with filled START_DATEs
    """
    # Get all time0 points for all records
    time0_dates_idx = finaldf["INDEX_CANCER"] == "Yes"
    # Make all time0 points 0
    finaldf["diagnosis_int"] = finaldf["START_DATE"]
    finaldf["START_DATE"][time0_dates_idx] = 0
    # Remove unused variable
    del finaldf["diagnosis_int"]
    return finaldf


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


class TimelineImagingTransform(Transforms):
    pass
