from abc import ABCMeta
from dataclasses import dataclass
from functools import cached_property

import pandas as pd
from synapseclient import Synapse

from geniesp.config import BpcConfig


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
        timeline_infodf = pd.concat(
            [
                timeline_infodf,
                pd.DataFrame(
                    [
                        {
                            "code": "rt_rt_int",
                            "sampleType": "TIMELINE-TREATMENT-RT",
                            "dataset": "Cancer-Directed Radiation Therapy dataset",
                            "cbio": "TEMP",
                        },
                        {
                            "code": "redcap_ca_index",
                            "sampleType": "TIMELINE-SEQUENCE",
                            "dataset": "Cancer-level dataset",
                            "cbio": "INDEX_CANCER",
                            "id": "syn22296816",  # HACK: hard coded synapse id
                        },
                        {
                            "code": "dob_ca_dx_days",
                            "sampleType": "TIMELINE-SEQUENCE",
                            "dataset": "Cancer-level dataset",
                            "cbio": "CA_DX_DAYS",
                            "id": "syn22296816",  # HACK: hard coded synapse id
                        },
                        {
                            "code": "dob_cpt_report_days",
                            "sampleType": "TIMELINE-SEQUENCE",
                            "dataset": "Cancer panel test level dataset",
                            "cbio": "DPT_REPORT_DAYS",
                            "id": "syn22296816",  # HACK: hard coded synapse id
                        },
                        {
                            "code": "first_index_ca_days",
                            "sampleType": "SURVIVAL",
                            "dataset": "Cancer-level index dataset",
                            "cbio": "CANCER_INDEX",
                        },
                        {
                            "code": "redcap_ca_index",
                            "sampleType": "PATIENT",
                            "dataset": "Cancer-level dataset",
                        }
                    ],
                    index=["rt_rt_int", "redcap_ca_index", "dob_ca_dx_days", "dob_cpt_report_days", "first_index_ca_days", "redcap_ca_index"],
                ),
            ]
        )
        subset_infodf = timeline_infodf[timeline_infodf["sampleType"] == self.sample_type]
        subset_infodf.index = subset_infodf["code"]

        return subset_infodf

    def get_derived_variable_files(
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
        used_entities = []
        derived_variable_entities = {}
        # Group by dataset because different datasets could have the
        # same variable
        datasets = self.timeline_infodf.groupby("dataset")
        for dataset, df in datasets:
            # Get synapse id
            synid = df["id"].unique()[0]
            file_entity = self.syn.get(synid)
            derived_variable_entities[dataset] = file_entity
            used_entities.append(f"{synid}.{file_entity.versionNumber}")

        return {"derived_variable_entities": derived_variable_entities, "used": used_entities}


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
