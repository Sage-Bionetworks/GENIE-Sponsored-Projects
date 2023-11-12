from dataclasses import dataclass
import os
import pandas as pd

import synapseclient
from genie import process_functions

from config import BpcConfig
from extract import Extract
from transforms import Transforms, TimelinePerformanceTransform



def write_and_storedf(
    df: pd.DataFrame, filepath: str
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

# @dataclass
# class ETL:
#     config: BpcConfig
#     extraction: Extract
#     transformations: Transforms
#     loading: Load

#     def workflow(self):
#         syn = synapseclient.login() 
#         cohort = "BLADDER"
#         _REDCAP_TO_CBIOMAPPING_SYNID = "syn25712693.49"
#         _DATA_TABLE_IDS = "syn22296821"
#         self.extraction.extract()
#         self.transformations.transform()
#         self.loading.write_and_storedf(
#             df=performance_data["df"],
#             filepath=os.path.join(
#                 cohort, "data_timeline_performance_status.txt"
#             )
#         )

syn = synapseclient.login()
cohort = "BLADDER"

config = BpcConfig()
temp_extract = Extract(
    bpc_config = config,
    sample_type = "TIMELINE-PERFORMANCE",
    syn = syn
)
temp_transform = TimelinePerformanceTransform(
    timeline_data= temp_extract.map_to_cbioportal_format(),
    timeline_infodf= temp_extract.timeline_infodf,
    extract = temp_extract
)

performance_data = temp_transform.create_timeline_file()

# Retraction...
write_and_storedf(
    df=performance_data,
    filepath=os.path.join(
        cohort, "data_timeline_performance_status.txt"
    )
)