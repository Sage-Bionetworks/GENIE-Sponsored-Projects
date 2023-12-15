"""GENIE SP/BPC cBioPortal exporter CLI"""
import argparse
import logging
import json
import os
import subprocess

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
    PatientTransform,
    MainGenie
)
from geniesp.utils import create_release_folders
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

    create_release_folders(args.sp)
    bpc_conf = BPC_MAPPING[args.sp]
    bpc_config = bpc_conf(syn=syn, cohort=bpc_conf.cohort, exclude_files=bpc_conf.exclude_files)
    logging.info("using config")
    logging.info(json.dumps(bpc_config.to_dict(), indent=4))
    # Create a mapping between the timeline file types and the transform classes
    timeline_files = {
        "TIMELINE-PERFORMANCE": {
            'cls': TimelinePerformanceTransform,
            'filename': "data_timeline_performance_status.txt"
        },
        "TIMELINE-TREATMENT-RT": {
            'cls': TimelineTreatmentRadTransform,
            'filename': "data_timeline_treatment.txt"
        },
        "TIMELINE-DX": {'cls': TimelineDxTransform, 'filename': "data_timeline_cancer_diagnosis.txt"},
        "TIMELINE-IMAGING": {'cls': TimelineTransform, 'filename': "data_timeline_imaging.txt"},
        "TIMELINE-MEDONC": {'cls': TimelineTransform, 'filename': "data_timeline_medonc.txt"},
        "TIMELINE-PATHOLOGY": {'cls': TimelineTransform, 'filename': "data_timeline_pathology.txt"},
        "TIMELINE-SAMPLE": {'cls': TimelineSampleTransform, 'filename': "data_timeline_sample_acquisition.txt"},
        "TIMELINE-SEQUENCE": {'cls': TimelineSequenceTransform, 'filename': "data_timeline_sequencing.txt"},
        "TIMELINE-LAB": {'cls': TimelineTransform, 'filename': "data_timeline_labtest.txt"},
        "SURVIVAL": {'cls': SurvivalTransform, 'filename': "data_clinical_supp_survival.txt"},
        "REGIMEN": {'cls': SurvivalTreatmentTransform, 'filename': "data_clinical_supp_survival_treatment.txt"},
        "SAMPLE": {'cls': SampleTransform, 'filename': "data_clinical_sample.txt"},
        "PATIENT": {'cls': PatientTransform, 'filename': "data_clinical_patient.txt"},

    }
    # Exception for timeline treatment file
    extract_raw_treatment_data = Extract(
        bpc_config = bpc_config,
        sample_type = "TIMELINE-TREATMENT",
        syn = syn
    )
    treatment_transform = TimelineTreatmentTransform(
        extract = extract_raw_treatment_data,
        bpc_config = bpc_config,
        filepath = os.path.join(args.sp, "data_timeline_treatment.txt"),
        sample_type = "TIMELINE-TREATMENT"
    )
    timeline_treatment_df = treatment_transform.create_timeline_file()
    # TODO fix this later
    # sample_type_dfs = {"TIMELINE-TREATMENT": timeline_treatment_df}
    sample_type_dfs = {}
    # TODO don't create folders if upload is false
    cbioportal_folders = get_cbioportal_upload_folders(
        syn=syn,
        staging_release_folder=bpc_config.staging_release_folder,
        cohort=args.sp,
        release=args.release
    )
    for sample_type, transform_cls in timeline_files.items():
        # Conditions to skip
        if sample_type == "TIMELINE-LAB" and args.sp in ["NSCLC", "BLADDER"]:
            logging.info(f"skipping {sample_type}...")
            continue

        if sample_type == "TIMELINE-PERFORMANCE" and args.sp not in ['BLADDER']:
            logging.info(f"skipping {sample_type}...")
            continue

        logging.info(f"writing {sample_type}...")
        # Download all the files required for processing
        extract_for_sample_type = Extract(
            bpc_config = bpc_config,
            sample_type = sample_type,
            syn = syn
        )
        derived_variables = extract_for_sample_type.get_derived_variable_files()
        # Leverage the cbioportal mapping and derived variables to create the timeline files
        filepath = os.path.join(bpc_config.cohort, transform_cls['filename'])
        # Transformation class
        sample_type_transform_cls = transform_cls['cls'](
            extract = extract_for_sample_type,
            bpc_config = bpc_config,
            filepath = filepath,
            sample_type=sample_type
        )

        # HACK this is because timeline treatment RT isn't created for two cohorts
        if sample_type == "TIMELINE-TREATMENT-RT" and args.sp in ["BrCa", "CRC"]:
            sample_type_df = pd.DataFrame()
        else:
            sample_type_df = sample_type_transform_cls.create_timeline_file()
        # This is specific to the timeline treatment file where it is concatenated with
        # the timeline file
        if sample_type == "TIMELINE-TREATMENT-RT":
            sample_type_df = pd.concat(
                [timeline_treatment_df, sample_type_df]
            )
        sample_type_dfs[sample_type] = sample_type_df
        # write the dataframe
        sample_type_transform_cls.write(df = sample_type_df)
        # store the files with provenance
        # Generate provenance
        used_entities = [
            bpc_config.redcap_to_cbio_mapping_synid,
            bpc_config.data_tables_id,
            bpc_config.mg_release_synid,
            bpc_config.prissmm_synid,
            bpc_config.sample_retraction_synid,
            bpc_config.patient_retraction_synid,
            bpc_config.retraction_at_release_synid,
            bpc_config.temporary_patient_retraction_synid,
            bpc_config.mg_assay_synid,
        ]
        used_entities.extend(derived_variables['used'])

        if args.upload:
            ent = synapseclient.File(
                filepath, parent=cbioportal_folders['release']
            )
            ent = syn.store(
                ent, used=used_entities, executed=bpc_config.github_url
            )
    MainGenie(
        bpc_config = bpc_config,
        extract = extract_raw_treatment_data,
        sample_df = sample_type_dfs['SAMPLE'],
        patient_df = sample_type_dfs['PATIENT'],
        release = args.release,
        cbioportal_folders = cbioportal_folders,
        syn = syn,
        upload = args.upload
    ).run()
    logging.info("cBioPortal validation")
    cmd = [
        "python",
        os.path.join(
            cbiopath, "core/src/main/scripts/importer/validateData.py"
        ),
        "-s",
        args.sp,
        "-n",
    ]
    subprocess.run(cmd)

if __name__ == "__main__":
    main()
