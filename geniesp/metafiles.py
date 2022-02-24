from typing import List
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


def create_clinical_meta_file(
    study_identifier: str,
    alteration_type: str,
    datatype: str,
    filename: str
) -> dict:
    """Create cbio meta files

    Args:
        datatype (str): Can be SAMPLE_ATTRIBUTE or PATIENT_ATTRIBUTE
        study_identifier (str): A string used to uniquely identify this cancer
                                study within the database
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
        study_identifier (str): A string used to uniquely identify this cancer study
                                within the database
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


def create_seg_meta_file(
    study_identifier: str,
    alteration_type: str,
    datatype: str,
    reference_genome_id: str,
    description: str,
    filename: str
) -> dict:
    """Create seg metadata files

    Args:
        study_identifier (str): A string used to uniquely identify this cancer study
                                within the database
        alteration_type (str): [description]
        datatype (str): [description]
        reference_genome_id (str): [description]
        description (str): [description]
        filename (str): [description]

    Returns:
        dict: cBioPortal meta file
    """
    meta_info = {
        "cancer_study_identifier": study_identifier,
        "genetic_alteration_type": alteration_type,
        "datatype": datatype,
        "reference_genome_id": reference_genome_id,
        "description": description,
        "data_filename": filename
    }
    write_meta_file(meta_info, filename)
    return meta_info


def create_meta_study(
    study_identifier: str,
    type_of_cancer: str,
    name: str,
    description: str,
    groups: str,
    short_name: str
) -> dict:
    """Create seg metadata files

    Args:
        study_identifier (str): A string used to uniquely identify this cancer study
                                within the database
        type_of_cancer (str): The cancer type abbreviation
        name (str): [description]
        description (str): [description]
        groups (str): [description]
        short_name (str): [description]

    Returns:
        dict: cBioPortal meta file
    """
    meta_info = {
        "type_of_cancer": type_of_cancer,
        "cancer_study_identifier": study_identifier,
        "name": name,
        "description": description,
        "groups": groups,
        "short_name": short_name
    }
    write_meta_file(meta_info, "meta_study.txt")
    return meta_info


def create_cbio_metafiles(
    study_identifier: str,
    cbio_fileformats: List[str] = [
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
    """Create cbioportal metadata files

    Args:
        study_identifier (str): A string used to uniquely identify this cancer study
                                within the database
        cbio_fileformats (list): List of cbioportal file names
    """
    # supp_survival* files don't need a meta file
    for cbio_file in cbio_fileformats:
        if cbio_file.startswith("data_clinical_sample"):
            create_clinical_meta_file(
                study_identifier=study_identifier,
                alteration_type="CLINICAL",
                datatype="SAMPLE_ATTRIBUTE",
                filename="data_clinical_sample.txt"
            )
        elif cbio_file.startswith("data_clinical_patient"):
            create_clinical_meta_file(
                study_identifier=study_identifier,
                alteration_type="CLINICAL",
                datatype="PATIENT_ATTRIBUTE",
                filename="data_clinical_patient.txt"
            )
        elif cbio_file.startswith("data_gene_matrix"):
            create_clinical_meta_file(
                study_identifier=study_identifier,
                alteration_type="GENE_PANEL_MATRIX",
                datatype="GENE_PANEL_MATRIX",
                filename="data_gene_matrix.txt"
            )
        elif cbio_file.startswith([
            "data_timeline_cancer_diagnosis.txt",
            "data_timeline_imaging.txt",
            "data_timeline_medonc.txt",
            "data_timeline_pathology.txt",
            "data_timeline_sample_acquisition.txt",
            "data_timeline_sequencing.txt",
            "data_timeline_treatment.txt"
        ]):
            create_clinical_meta_file(
                study_identifier=study_identifier,
                alteration_type="CLINICAL",
                datatype="TIMELINE",
                filename=cbio_file
            )
        elif cbio_file.startswith("data_CNA"):
            create_genomic_meta_file(
                study_identifier=study_identifier,
                alteration_type="COPY_NUMBER_ALTERATION",
                datatype="DISCRETE",
                stable_id="cna",
                profile_name="Copy-number alterations",
                profile_description="Copy-number alterations",
                filename="data_CNA.txt"
            )
        elif cbio_file.startswith("data_fusions"):
            create_genomic_meta_file(
                study_identifier=study_identifier,
                alteration_type="FUSION",
                datatype="FUSION",
                stable_id="fusion",
                profile_name="Fusions",
                profile_description="Fusions",
                filename="data_fusions.txt"
            )
        elif cbio_file.startswith("data_mutations_extended"):
            create_genomic_meta_file(
                study_identifier=study_identifier,
                alteration_type="MUTATIONS_EXTENDED",
                datatype="MAF",
                stable_id="mutations",
                profile_name="Mutations",
                profile_description="Mutation data from next-gen sequencing.",
                filename="data_mutations_extended.txt"
            )
        elif cbio_file.endswith("data_cna_hg19.seg"):
            create_seg_meta_file(
                study_identifier=study_identifier,
                alteration_type="COPY_NUMBER_ALTERATION",
                datatype="SEG",
                reference_genome_id="hg19",
                description="Segment data for the genie study",
                filename=cbio_file
            )
        else:
            print(f"{cbio_file} does not have associated metafile")
        create_meta_study(
            study_identifier=study_identifier,
            type_of_cancer="mixed",
            name="PLACEHOLDER",
            description="PLACEHOLDER",
            groups="GENIE",
            short_name="PLACEHOLDER"
        )
