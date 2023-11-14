import os
import pandas as pd

import synapseclient
from genie import process_functions

# from geniesp.config import BpcConfig
# from geniesp.extract import Extract
# from geniesp.transforms import (
#     TimelinePerformanceTransform,
#     TimelineTreatmentRadTransform,
#     TimelineTreatmentTransform,
#     TimelineTransform,
#     TimelineSampleTransform,
#     TimelineSequenceTransform,
#     TimelineDxTransform,
#     SurvivalTransform,
#     SurvivalTreatmentTransform,
#     SampleTransform,
#     PatientTransform
# )


def write_and_storedf(
    syn, df: pd.DataFrame, filepath: str, used_entities: list = []
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
    ent = synapseclient.File(filepath, parent="syn52950402")
    syn.store(ent, executed="https://github.com/Sage-Bionetworks/GENIE-Sponsored-Projects", used=used_entities)


def write_clinical_file(
    clinicaldf: pd.DataFrame,
    redcap_to_cbiomappingdf: pd.DataFrame,
    clin_path: str
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
    redcap_to_cbiomappingdf.index = redcap_to_cbiomappingdf["cbio"]
    label_map = redcap_to_cbiomappingdf["labels"].to_dict()
    description_map = redcap_to_cbiomappingdf["description"].to_dict()
    coltype_map = redcap_to_cbiomappingdf["colType"].to_dict()
    # Priority column will determine which columns are shown on cBioPortal
    # Must columns should be shown on cBioPortal will have a 1
    # but survival columns should all be 0
    priority_map = redcap_to_cbiomappingdf["priority"].to_dict()

    labels = [str(label_map[col]) for col in clinicaldf]
    descriptions = [str(description_map[col]) for col in clinicaldf]
    coltype = [str(coltype_map[col]) for col in clinicaldf]
    priority = [str(int(priority_map[col])) for col in clinicaldf]

    # clin_path = os.path.join(
    #     self._SPONSORED_PROJECT, f"data_clinical_{filetype}.txt"
    # )

    with open(clin_path, "w+") as clin_file:
        clin_file.write("#{}\n".format("\t".join(labels)))
        clin_file.write("#{}\n".format("\t".join(descriptions)))
        clin_file.write("#{}\n".format("\t".join(coltype)))
        # TODO attributes in the supp file are PATIENT, so must
        # specify that.
        # if filetype.startswith("supp_survival"):
        #     clin_file.write("#{}\n".format("\t".join(["PATIENT"] * len(labels))))
        clin_file.write("#{}\n".format("\t".join(priority)))
        clin_file.write(
            process_functions.removeStringFloat(
                clinicaldf.to_csv(index=False, sep="\t")
            )
        )
    return clin_path
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

# syn = synapseclient.login()
# cohort = "BLADDER"

# config = BpcConfig(
#     cohort = cohort
# )
# timeline_files = {
#     # "TIMELINE-PERFORMANCE": TimelinePerformanceTransform,
#     # "TIMELINE-TREATMENT-RT":  TimelineTreatmentRadTransform,
#     # "TIMELINE-DX": TimelineDxTransform,
#     # "TIMELINE-IMAGING": TimelineTransform,
#     # "TIMELINE-MEDONC": TimelineTransform,
#     # "TIMELINE-PATHOLOGY": TimelineTransform,
#     # "TIMELINE-SAMPLE": TimelineSampleTransform,
#     # "TIMELINE-SEQUENCE": TimelineSequenceTransform,
#     # "TIMELINE-LAB": TimelineTransform,
#     # "SURVIVAL": SurvivalTransform,
#     "SAMPLE": SampleTransform,
#     "PATIENT": PatientTransform,

# }
# # Exception for timeline treatment file
# temp_extract = Extract(
#     bpc_config = config,
#     sample_type = "TIMELINE-TREATMENT",
#     syn = syn
# )
# temp_transform = TimelineTreatmentTransform(
#     # timeline_infodf= temp_extract.timeline_infodf,
#     extract = temp_extract,
#     bpc_config = config
# )

# timeline_treatment_df = temp_transform.create_timeline_file()

# for sample_type, transform_cls in timeline_files.items():
#     # Conditions to skip
#     if sample_type == "TIMELINE-LAB" and cohort in ["NSCLC", "BLADDER"]:
#         continue
#     if sample_type == "TIMELINE-PERFORMANCE" and cohort not in ['BLADDER']:
#         continue
#     if sample_type == "TIMELINE-TREATMENT-RT" and cohort in ["BrCa", "CRC"]:
#         continue

#     # Download all the files required for processing
#     temp_extract = Extract(
#         bpc_config = config,
#         sample_type = sample_type,
#         syn = syn
#     )
#     derived_variables = temp_extract.get_derived_variable_files()

#     # Leverage the cbioportal mapping and derived variables to create the timeline files
#     temp_transform = transform_cls(
#         # timeline_infodf= temp_extract.timeline_infodf,
#         extract = temp_extract,
#         bpc_config = config
#     )
#     if 'TIMELINE-DX':
#         filter_start = False
#     else:
#         filter_start = True

#     performance_data = temp_transform.create_timeline_file(filter_start=filter_start)
#     # This is specific to the timeline treatment file where it is concatenated with
#     # the timeline file
#     if sample_type == "TIMELINE-TREATMENT-RT":
#         performance_data = pd.concat(
#             [timeline_treatment_df, performance_data]
#         )
#     # store the files with provenance
#     # Generate provenance
#     used_entities = [
#         config.redcap_to_cbio_mapping_synid,
#         config.data_tables_id,
#         config.mg_release_synid,
#         config.prissmm_synid,
#         config.sample_retraction_synid,
#         config.patient_retraction_synid,
#         config.retraction_at_release_synid,
#         config.temporary_patient_retraction_synid,
#         config.mg_assay_synid,
#     ]
#     used_entities.extend(derived_variables['used'])
#     write_and_storedf(
#         syn=syn,
#         df=performance_data,
#         filepath=os.path.join(cohort, f"{sample_type}.txt"),
#         used_entities = used_entities
#     )
#     # TODO: need to have a write clinical method