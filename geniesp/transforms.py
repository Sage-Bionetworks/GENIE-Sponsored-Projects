from abc import ABCMeta
from dataclasses import dataclass
import os

import pandas as pd
import synapseclient
from synapseclient import Synapse
from genie import process_functions

from extract import Extract


def _convert_to_int(value):
    """Convert object to integer or return nan"""
    try:
        return int(value)
    except ValueError:
        return float('nan')


@dataclass
class Transforms(metaclass=ABCMeta):
    """Timeline data class."""
    timeline_data: dict
    timeline_infodf: pd.DataFrame
    extract: Extract
    # config: BpcConfig

    def transforms(self, timelinedf, subset_infodf, filter_start):
        # Obtain portal value (EVENT_TYPE)
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
        timelinedf = self.transforms(timelinedf=self.timeline_data['df'], subset_infodf=self.timeline_infodf, filter_start=filter_start)
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
    