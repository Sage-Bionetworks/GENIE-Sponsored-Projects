import yaml


def write_meta_file(meta_info: dict, filename: str) -> str:
    """Writes metadata file

    Args:
        meta_info (dict): metafile in a dict
        filename (str): cbioportal filename

    Returns:
        str: meta file path
    """
    meta_filename = filename.replace("data", "meta")
    with open(meta_filename, "w") as meta_f:
        yaml.dump(meta_info, meta_f)
    return meta_filename


def create_meta_file(
    study_identifier: str,
    alteration_type: str,
    datatype: str,
    filename: str
) -> dict:
    """Create cbio meta files

    Args:
        datatype (str): Can be SAMPLE_ATTRIBUTE or PATIENT_ATTRIBUTE
        study_identifier (str): cbioportal study identifier
        filename (str): cbioportal filename
        alteration_type (str): Alteration type

    Returns:
        dict: metafile in a dict
    """
    # datatype is SAMPLE_ATTRIBUTE or PATIENT_ATTRIBUTE
    clinical_meta = {
        "cancer_study_identifier": study_identifier,
        "genetic_alteration_type": alteration_type,
        "datatype": datatype,
        "data_filename": filename
    }
    write_meta_file(clinical_meta, filename)
    return clinical_meta


def create_cna_meta_file(
    datatype: str,
    study_identifier: str,
    filename: str
) -> dict:
    meta_cna = {
        "cancer_study_identifier": "crc_genie_bpc",
        "genetic_alteration_type": "COPY_NUMBER_ALTERATION",
        "datatype": "DISCRETE",
        "stable_id": "cna",
        "show_profile_in_analysis_tab": "true",
        "profile_name": "Copy-number alterations",
        "profile_description": "Copy-number alterations",
        "data_filename": "data_CNA.txt"
    }


def main():
    """Create metadata files
    """
    # Create CRC meta files
    create_meta_file(
        study_identifier="crc_genie_bpc",
        alteration_type="CLINICAL",
        datatype="SAMPLE_ATTRIBUTE",
        filename="data_clinical_sample.txt"
    )
    create_meta_file(
        study_identifier="crc_genie_bpc",
        alteration_type="CLINICAL",
        datatype="PATIENT_ATTRIBUTE",
        filename="data_clinical_patient.txt"
    )
    create_meta_file(
        study_identifier="crc_genie_bpc",
        alteration_type="GENE_PANEL_MATRIX",
        datatype="GENE_PANEL_MATRIX",
        filename="data_gene_matrix.txt"
    )
    timeline_meta_files = [
        "data_timeline_cancer_diagnosis.txt",
        "data_timeline_imaging.txt",
        "data_timeline_medonc.txt",
        "data_timeline_pathology.txt",
        "data_timeline_sample_acquisition.txt",
        "data_timeline_sequencing.txt",
        "data_timeline_treatment"
    ]
    for filenames in timeline_meta_files:
        create_meta_file(
            study_identifier="crc_genie_bpc",
            alteration_type="CLINICAL",
            datatype="TIMELINE",
            filename=filenames
        )
