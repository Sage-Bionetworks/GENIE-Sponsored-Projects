from dataclasses import dataclass
import os
import pandas as pd

import synapseclient
from genie import process_functions

from config import BpcConfig
from extract import Extract
from transforms import (
    TimelinePerformanceTransform,
    TimelineTreatmentRadTransform,
    TimelineTreatmentTransform,
    TimelineTransform,
    TimelineSampleTransform,
    TimelineSequenceTransform
)



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
timeline_files = {
    # "TIMELINE-PERFORMANCE": TimelinePerformanceTransform,
    # "TIMELINE-TREATMENT-RT":  TimelineTreatmentRadTransform,
    # "TIMELINE-DX": TimelineDxTransform,
    # "TIMELINE-IMAGING": TimelineTransform,
    # "TIMELINE-MEDONC": TimelineTransform,
    # "TIMELINE-PATHOLOGY": TimelineTransform,
    # "TIMELINE-SAMPLE": TimelineSampleTransform,
    # "TIMELINE-SEQUENCE": TimelineSequenceTransform,
    "TIMELINE-LAB": TimelineTransform
}
# Exception for timeline treatment file
temp_extract = Extract(
    bpc_config = config,
    sample_type = "TIMELINE-TREATMENT",
    syn = syn
)
temp_transform = TimelineTreatmentTransform(
    # timeline_infodf= temp_extract.timeline_infodf,
    extract = temp_extract,
    bpc_config = config
)

timeline_treatment_df = temp_transform.create_timeline_file()

for sample_type, transform_cls in timeline_files.items():
    if sample_type == "TIMELINE-LAB" and cohort in ["NSCLC", "BLADDER"]:
        continue
    temp_extract = Extract(
        bpc_config = config,
        sample_type = sample_type,
        syn = syn
    )
    temp_transform = transform_cls(
        # timeline_infodf= temp_extract.timeline_infodf,
        extract = temp_extract,
        bpc_config = config
    )
    if 'TIMELINE-DX':
        filter_start = False
    else:
        filter_start = True
    performance_data = temp_transform.create_timeline_file(filter_start=filter_start)
    if sample_type == "TIMELINE-TREATMENT-RT":
        performance_data = pd.concat(
            [timeline_treatment_df, performance_data]
        )
    # Retraction...
    write_and_storedf(
        df=performance_data,
        filepath=os.path.join(
            cohort, f"{sample_type}.txt"
        )
    )
