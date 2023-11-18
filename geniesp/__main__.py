"""GENIE SP/BPC cBioPortal exporter CLI"""
import argparse
import logging
import os

import pandas as pd
import synapseclient

from geniesp.config import Brca, Crc, Nsclc, Panc, Prostate, Bladder
from geniesp.extract import Extract
from geniesp.transforms import (
    TimelinePerformanceTransform,
    TimelineTreatmentRadTransform,
    TimelineTreatmentTransform,
    TimelineTransform,
    TimelineSampleTransform,
    TimelineSequenceTransform,
    TimelineDxTransform,
    SurvivalTransform,
    SurvivalTreatmentTransform,
    SampleTransform,
    PatientTransform
)
from geniesp.loads import get_cbioportal_upload_folders


BPC_MAPPING = {
    "NSCLC": Nsclc,
    "CRC": Crc,
    "BrCa": Brca,
    "PANC": Panc,
    "Prostate": Prostate,
    "BLADDER": Bladder
}


def build_parser():
    parser = argparse.ArgumentParser(description="Run GENIE sponsored projects")
    parser.add_argument(
        "sp",
        type=str,
        help="Specify sponsored project to run",
        choices=BPC_MAPPING.keys(),
    )

    parser.add_argument("release", type=str, help="Specify bpc release")
    parser.add_argument(
        "--upload",
        action="store_true",
        help="Upload files into Synapse BPC staging directory. Default: false.",
    )
    parser.add_argument(
        "--log",
        "-l",
        type=str,
        choices=["debug", "info", "warning", "error"],
        default="info",
        help="Set logging output level " "(default: %(default)s)",
    )
    parser.add_argument(
        "--cbioportal",
        type=str,
        help="Optional parameter to specify cbioportal folder location",
    )
    return parser.parse_args()

def main():
    """Main"""
    args = build_parser()

    numeric_level = getattr(logging, args.log.upper(), None)
    if not isinstance(numeric_level, int):
        raise ValueError(f"Invalid log level: args.log")
    logging.basicConfig(level=numeric_level)

    syn = synapseclient.login()

    if args.cbioportal is None:
        cbiopath = "../cbioportal"
    else:
        cbiopath = args.cbioportal

    config = BPC_MAPPING[args.sp]

    timeline_files = {
        "TIMELINE-PERFORMANCE": TimelinePerformanceTransform,
        "TIMELINE-TREATMENT-RT":  TimelineTreatmentRadTransform,
        "TIMELINE-DX": TimelineDxTransform,
        "TIMELINE-IMAGING": TimelineTransform,
        "TIMELINE-MEDONC": TimelineTransform,
        "TIMELINE-PATHOLOGY": TimelineTransform,
        "TIMELINE-SAMPLE": TimelineSampleTransform,
        "TIMELINE-SEQUENCE": TimelineSequenceTransform,
        "TIMELINE-LAB": TimelineTransform,
        "SURVIVAL": SurvivalTransform,
        "REGIMEN": SurvivalTreatmentTransform,
        "SAMPLE": SampleTransform,
        "PATIENT": PatientTransform,

    }
    # Exception for timeline treatment file
    temp_extract = Extract(
        bpc_config = config,
        sample_type = "TIMELINE-TREATMENT",
        syn = syn
    )
    temp_transform = TimelineTreatmentTransform(
        extract = temp_extract,
        bpc_config = config

    )

    timeline_treatment_df = temp_transform.create_timeline_file()
    sample_type_dfs = {"TIMELINE-TREATMENT": timeline_treatment_df}

    for sample_type, transform_cls in timeline_files.items():
        # Conditions to skip
        if sample_type == "TIMELINE-LAB" and args.sp in ["NSCLC", "BLADDER"]:
            continue
        if sample_type == "TIMELINE-PERFORMANCE" and args.sp not in ['BLADDER']:
            continue
        if sample_type == "TIMELINE-TREATMENT-RT" and args.sp in ["BrCa", "CRC"]:
            continue

        # Download all the files required for processing
        temp_extract = Extract(
            bpc_config = config,
            sample_type = sample_type,
            syn = syn
        )
        derived_variables = temp_extract.get_derived_variable_files()
        # Leverage the cbioportal mapping and derived variables to create the timeline files
        filepath = f"{sample_type}.txt"

        temp_transform = transform_cls(
            extract = temp_extract,
            bpc_config = config,
            filepath = filepath
        )
        if sample_type == 'TIMELINE-DX':
            filter_start = False
        else:
            filter_start = True

        performance_data = temp_transform.create_timeline_file(filter_start=filter_start)
        # This is specific to the timeline treatment file where it is concatenated with
        # the timeline file
        if sample_type == "TIMELINE-TREATMENT-RT":
            performance_data = pd.concat(
                [timeline_treatment_df, performance_data]
            )
        sample_type_dfs[sample_type] = performance_data
        # write the dataframe
        temp_transform.write(df = performance_data)
        # store the files with provenance
        # Generate provenance
        used_entities = [
            config.redcap_to_cbio_mapping_synid,
            config.data_tables_id,
            config.mg_release_synid,
            config.prissmm_synid,
            config.sample_retraction_synid,
            config.patient_retraction_synid,
            config.retraction_at_release_synid,
            config.temporary_patient_retraction_synid,
            config.mg_assay_synid,
        ]
        used_entities.extend(derived_variables['used'])
        cbioportal_folders = get_cbioportal_upload_folders(
            syn=syn,
            staging_release_folder=config.staging_release_folder,
            cohort=args.cohort,
            release=args.release
        )
        ent = synapseclient.File(
            filepath, parent=cbioportal_folders['release']
        )
        if args.upload:
            ent = syn.store(
                ent, used=used_entities, executed=config.github_url
            )
    from geniesp.transforms import MainGenie
    MainGenie(
        bpc_config = config,
        extract = temp_extract,
        sample_df = sample_type_dfs['SAMPLE'],
        patient_df = sample_type_dfs['PATIENT'],
        release = args.release,
        cbioportal_folders = cbioportal_folders,
        syn = syn,
        upload = args.upload
    ).run()


if __name__ == "__main__":
    main()
