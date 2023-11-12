from abc import ABCMeta
from dataclasses import dataclass
from functools import cached_property

import pandas as pd
from synapseclient import Synapse

from config import BpcConfig


@dataclass
class Extract:
    """Timeline data class."""
    # Config = TimelineConfig
    bpc_config: BpcConfig
    syn: Synapse
    sample_type: str = None

    def get_mg_synid(self, synid_folder: str, file_name: str) -> str:
        """Get Synapse ID of main GENIE data file in release folder.

        Args:
            synid_folder (str): Synapse ID of main GENIE release folder
            file_name (str): File name for which to retrieve Synapse ID

        Returns:
            str: Synapse ID if file found; otherwise, None
        """
        synid_children = self.syn.getChildren(synid_folder)
        for synid_child in synid_children:
            if synid_child["name"] == file_name:
                return synid_child["id"]
        raise ValueError(f"file '{file_name}' not found in {synid_folder}")

    @cached_property
    def mapping_df(self) -> pd.DataFrame:
        """Extract relevant portions of the mapping table for the sponsored project
        variables and cBioPortal variable mappings and return as a data frame.

        Args:
            syn (Synapse): Synapse connection
            cohort (str): sponsored project label
            synid_table_cbio (str): Synapse ID of table containing variable to cbio mapping

        Returns:
            pd.DataFrame: data frame of all mapping columns for released variables
        """
        redcap_to_cbiomapping = self.syn.tableQuery(
            f"SELECT * FROM {self.bpc_config.redcap_to_cbio_mapping_synid} where "
            f"{self.bpc_config.cohort} is true AND sampleType <> 'TIMELINE-STATUS'"
        )
        redcap_to_cbiomappingdf = redcap_to_cbiomapping.asDataFrame()
        return redcap_to_cbiomappingdf

    @cached_property
    def data_tables_df(self) -> pd.DataFrame:
        """Get mapping of dataset files to Synapse IDs.

        Args:
            syn (Synapse): Synapse connection
            synid_table_files (str): Synapse ID of table containing data file to Synapse ID mapping

        Returns:
            pd.DataFrame: data frame with two columns representing the Synapse ID and dataset label
        """
        data_tables = self.syn.tableQuery(f"SELECT id, dataset FROM {self.bpc_config.data_tables_id}")
        data_tablesdf = data_tables.asDataFrame()
        return data_tablesdf

    @cached_property
    def timeline_infodf(self):
        timeline_infodf = self.mapping_df.query(f'sampleType == "{self.sample_type}"').merge(
            self.data_tables_df, on="dataset", how="left"
        )
        timeline_infodf.index = timeline_infodf["code"]
        timeline_infodf = pd.concat(
            [
                timeline_infodf,
                pd.DataFrame(
                    {
                        "code": "rt_rt_int",
                        "sampleType": "TIMELINE-TREATMENT-RT",
                        "dataset": "Cancer-Directed Radiation Therapy dataset",
                        "cbio": "TEMP",
                    },
                    index=["rt_rt_int"],
                ),
            ]
        )
        subset_infodf = timeline_infodf[timeline_infodf["sampleType"] == self.sample_type]

        return subset_infodf

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
        print(self.timeline_infodf)
        mappingdf = self.timeline_infodf[self.timeline_infodf["data_type"] != "portal_value"]
        # Group by dataset because different datasets could have the
        # same variable
        datasets = mappingdf.groupby("dataset")
        finaldf = pd.DataFrame()
        used_entities = []

        for _, df in datasets:
            # Get synapse id
            synid = df["id"].unique()[0]
            table = self.syn.get(synid)
            used_entities.append(f"{synid}.{table.versionNumber}")
            # obtain columns to subset df
            cols = df["code"][df["sampleType"] == self.sample_type]
            cols = cols.tolist()
            if "record_id" not in cols:
                cols.append("record_id")
            # Must add path_rep_number for sample and sample acquisition file
            if self.sample_type in ["SAMPLE", "TIMELINE-SAMPLE"]:
                cols.append("path_rep_number")
            # Must add path_proc_number to sample file
            if self.sample_type == "SAMPLE":
                cols.append("path_proc_number")
            # Only get specific cohort and subset cols
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

        return {"df": finaldf, "used": used_entities}


    @cached_property
    def to_keep_samples_and_patients(self) -> pd.DataFrame:
        """Get main GENIE clinical samples and perform retraction
        against BPC sample and patient database. It's important that I use the
        sample database, because if a sample from a patient with multiple samples
        is retracted from the main GENIE clinical samples, the patient
        will still exist. (Jira: GEN-260)

        Returns:
            pd.DataFrame: main GENIE information for sponsored project samples
        """
        # genie_clinicaldb = self.syn.tableQuery(
        #     f"select SAMPLE_ID, PATIENT_ID, ONCOTREE_CODE, SEQ_ASSAY_ID, "
        #     f"SAMPLE_TYPE, SEQ_YEAR from {self._CLINICAL_SYNID}"
        # )
        # genie_clinicaldf = genie_clinicaldb.asDataFrame()
        # DECISION 06/2023: use the release file to do retractions
        # This is due to the most recent releases potentially having
        # samples retracted. The consortium release matched with the
        # public release (14.7-consortium <-> 14.0-public) must be used
        # due to SEQ_YEAR being used through the code.
        sample_synid = self.get_mg_synid(
            self.bpc_config.mg_release_synid, "data_clinical_sample.txt"
        )
        genie_clinicaldf = pd.read_csv(
            self.syn.get(sample_synid, followLink=True).path, sep="\t", comment="#"
        )
        # Filter out cfDNA samples
        genie_clinicaldf = genie_clinicaldf[
            genie_clinicaldf['SAMPLE_CLASS'] != "cfDNA"
        ]
        # BPC retraction database
        bpc_sample_retraction_db = self.syn.tableQuery(
            f"select SAMPLE_ID from {self.bpc_config.sample_retraction_synid} where "
            f"{self.bpc_config.cohort} is true"
        )
        bpc_sample_retractiondf = bpc_sample_retraction_db.asDataFrame()

        bpc_patient_retraction_db = self.syn.tableQuery(
            f"select record_id from {self.bpc_config.patient_retraction_synid} where "
            f"{self.bpc_config.cohort} is true"
        )
        bpc_patient_retraction_df = bpc_patient_retraction_db.asDataFrame()

        bpc_temp_patient_retraction_db = self.syn.tableQuery(
            f"select record_id from {self.bpc_config.temporary_patient_retraction_synid} where "
            f"cohort = '{self.bpc_config.cohort}'"
        )
        bpc_temp_patient_retraction_df = bpc_temp_patient_retraction_db.asDataFrame()

        retraction_at_release = self.syn.tableQuery(
            f"select patient_id from {self.bpc_config.retraction_at_release_synid} where "
            f"cohort = '{self.bpc_config.cohort}'"
        )
        retraction_at_release_df = retraction_at_release.asDataFrame()
        # Retract samples from sample retraction db
        keep_clinicaldf = genie_clinicaldf[
            ~genie_clinicaldf["SAMPLE_ID"].isin(bpc_sample_retractiondf["SAMPLE_ID"])
        ]
        # Retract patients from patient retraction db
        keep_clinicaldf = keep_clinicaldf[
            ~keep_clinicaldf["PATIENT_ID"].isin(bpc_patient_retraction_df["record_id"])
        ]
        # Retract patients for at release retraction table
        keep_clinicaldf = keep_clinicaldf[
            ~keep_clinicaldf["PATIENT_ID"].isin(retraction_at_release_df["patient_id"])
        ]
        # Retract patients from temporary patient retraction db
        keep_clinicaldf = keep_clinicaldf[
            ~keep_clinicaldf["PATIENT_ID"].isin(
                bpc_temp_patient_retraction_df["record_id"]
            )
        ]
        return keep_clinicaldf