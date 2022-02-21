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
        dict: cBioPortal meta file
    """
    # datatype is SAMPLE_ATTRIBUTE or PATIENT_ATTRIBUTE
    meta_info = {
        "cancer_study_identifier": study_identifier,
        "genetic_alteration_type": alteration_type,
        "datatype": datatype,
        "data_filename": filename
    }
    write_meta_file(meta_info, filename)
    return meta_info


def create_genomic_meta_file(
    study_identifier: str,
    alteration_type: str,
    datatype: str,
    stable_id: str,
    profile_name: str,
    profile_description: str,
    filename: str
) -> dict:
    """Create genomic metadata files

    Args:
        study_identifier (str): [description]
        alteration_type (str): [description]
        datatype (str): [description]
        stable_id (str): [description]
        profile_name (str): [description]
        profile_description (str): [description]
        filename (str): [description]

    Returns:
        dict: cBioPortal meta file
    """
    meta_info = {
        "cancer_study_identifier": study_identifier,
        "genetic_alteration_type": alteration_type,
        "datatype": datatype,
        "stable_id": stable_id,
        "show_profile_in_analysis_tab": "true",
        "profile_name": profile_name,
        "profile_description": profile_description,
        "data_filename": filename
    }
    write_meta_file(meta_info, filename)
    return meta_info


def main(
    study_identifier: str,
    cbio_fileformats: list = [
        "data_clinical_sample.txt",
        "data_clinical_patient.txt",
        "data_clinical_supp_survival.txt",
        "data_clinical_supp_survival_treatment.txt",
        "data_gene_matrix.txt",
        "data_timeline_cancer_diagnosis.txt",
        "data_timeline_imaging.txt",
        "data_timeline_medonc.txt",
        "data_timeline_pathology.txt",
        "data_timeline_sample_acquisition.txt",
        "data_timeline_sequencing.txt",
        "data_timeline_treatment.txt",
        "data_timeline_labtest.txt",
        "data_CNA.txt",
        "data_fusions.txt",
        "data_mutations_extended.txt"
    ]
):
    """Create metadata files

    Args:
        sp (str): cBioPortal study identifier
    """
    create_meta_file(
        study_identifier=study_identifier,
        alteration_type="CLINICAL",
        datatype="SAMPLE_ATTRIBUTE",
        filename="data_clinical_sample.txt"
    )
    create_meta_file(
        study_identifier=study_identifier,
        alteration_type="CLINICAL",
        datatype="PATIENT_ATTRIBUTE",
        filename="data_clinical_patient.txt"
    )
    # supp_survival* files don't need a meta file
    create_meta_file(
        study_identifier=study_identifier,
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
        "data_timeline_treatment.txt"
    ]
    for filenames in timeline_meta_files:
        create_meta_file(
            study_identifier=study_identifier,
            alteration_type="CLINICAL",
            datatype="TIMELINE",
            filename=filenames
        )

    create_genomic_meta_file(
        study_identifier=study_identifier,
        alteration_type="COPY_NUMBER_ALTERATION",
        datatype="DISCRETE",
        stable_id="cna",
        profile_name="Copy-number alterations",
        profile_description="Copy-number alterations",
        filename="data_CNA.txt"
    )
    create_genomic_meta_file(
        study_identifier=study_identifier,
        alteration_type="FUSION",
        datatype="FUSION",
        stable_id="fusion",
        profile_name="Fusions",
        profile_description="Fusions",
        filename="data_fusions.txt"
    )
    create_genomic_meta_file(
        study_identifier=study_identifier,
        alteration_type="MUTATIONS_EXTENDED",
        datatype="MAF",
        stable_id="mutations",
        profile_name="Mutations",
        profile_description="Mutation data from next-gen sequencing.",
        filename="data_mutations_extended.txt"
    )
