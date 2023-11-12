from dataclasses import dataclass
from functools import cached_property

import pandas as pd
import synapseclient
from synapseclient import Synapse


@dataclass
class BpcConfig():
    # json_path: str
    cohort = "BLADDER"
    # Redcap codes to cbioportal mapping synid and form key is in
    # version 38, 42 were last stable version(s)
    redcap_to_cbio_mapping_synid = "syn25712693.49"
    # Run `git rev-parse HEAD` in Genie_processing directory to obtain shadigest
    # github_repo = None
    # Mapping from Synapse Table to derived variables
    # TODO: Make versioned
    data_tables_id = "syn22296821"
    # Storage of not found samples
    sp_redcap_exports_synid = "syn21446571"
    # main GENIE release folder (14.8-consortium)
    # Must use consortium release, because SEQ_YEAR is used
    mg_release_synid = "syn52794528"
    # PRISSMM documentation table
    prissmm_synid = "syn22684834"
    # BPC sample retraction table
    sample_retraction_synid = "syn25779833"
    patient_retraction_synid = "syn25998970"
    retraction_at_release_synid = "syn52915299"
    temporary_patient_retraction_synid = "syn29266682"
    # main GENIE assay information table
    mg_assay_synid = "syn17009222"
    # exclude files to be created for cbioportal
    # TODO: need to support this feature in rest of code, for now
    # This is added for metadata files
    exclude_files = []
    # cohort-generic link to documentation for BPC datasets
    url_bpc = "https://aacr.box.com/s/en5dyu9zfw1krg2u58wlcz01jttc6y9h"
    # cohort-generic link to documentation for cBio files
    url_cbio = "https://docs.google.com/document/d/1IBVF-FLecUG8Od6mSEhYfWH3wATLNMnZcBw2_G0jSAo/edit"
    # syn: Synapse
    
    # def __post_init__(self):
    #     self.syn = synapseclient.login()

    # @cached_property
    # def mapping_df(self) -> pd.DataFrame:
    #     """Extract relevant portions of the mapping table for the sponsored project
    #     variables and cBioPortal variable mappings and return as a data frame.

    #     Args:
    #         syn (Synapse): Synapse connection
    #         cohort (str): sponsored project label
    #         synid_table_cbio (str): Synapse ID of table containing variable to cbio mapping

    #     Returns:
    #         pd.DataFrame: data frame of all mapping columns for released variables
    #     """
    #     redcap_to_cbiomapping = self.syn.tableQuery(
    #         f"SELECT * FROM {self.redcap_to_cbio_mapping_synid} where "
    #         f"{self.cohort} is true AND sampleType <> 'TIMELINE-STATUS'"
    #     )
    #     redcap_to_cbiomappingdf = redcap_to_cbiomapping.asDataFrame()
    #     return redcap_to_cbiomappingdf

    # @cached_property
    # def data_tables_df(self) -> pd.DataFrame:
    #     """Get mapping of dataset files to Synapse IDs.

    #     Args:
    #         syn (Synapse): Synapse connection
    #         synid_table_files (str): Synapse ID of table containing data file to Synapse ID mapping

    #     Returns:
    #         pd.DataFrame: data frame with two columns representing the Synapse ID and dataset label
    #     """
    #     data_tables = self.syn.tableQuery(f"SELECT id, dataset FROM {self.data_tables_id}")
    #     data_tablesdf = data_tables.asDataFrame()
    #     return data_tablesdf
