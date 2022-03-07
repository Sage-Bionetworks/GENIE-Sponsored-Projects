from typing import List

import yaml


def write_meta_file(meta_info: dict, filename: str, outdir: str) -> str:
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
    study_identifier: str, alteration_type: str, datatype: str, filename: str
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
        "data_filename": filename,
    }
    return meta_info


def create_genomic_meta_file(
    study_identifier: str,
    alteration_type: str,
    datatype: str,
    stable_id: str,
    profile_name: str,
    profile_description: str,
    filename: str,
) -> dict:
    """Create genomic metadata files

    Args:
        study_identifier (str): A string used to uniquely identify this cancer study
                                within the database
        alteration_type (str): Allowed values - MUTATIONS_EXTENDED, FUSION, COPY_NUMBER_ALTERATION
        datatype (str): Allowed values - DISCRETE, FUSION, MAF
        stable_id (str): Allowed values - mutations, fusion, cna
        profile_name (str): A name for the data
        profile_description (str): A description of the data
        filename (str): cbioportal filename

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
        "data_filename": filename,
    }
    return meta_info


def create_seg_meta_file(
    study_identifier: str, reference_genome_id: str, description: str, filename: str
) -> dict:
    """Create seg metadata files

    Args:
        study_identifier (str): A string used to uniquely identify this cancer study
                                within the database
        reference_genome_id (str): Reference genome version.
        description (str): A description of the segmented data.
        filename (str): cbioportal filename

    Returns:
        dict: cBioPortal meta file
    """
    meta_info = {
        "cancer_study_identifier": study_identifier,
        "genetic_alteration_type": "COPY_NUMBER_ALTERATION",
        "datatype": "SEG",
        "reference_genome_id": reference_genome_id,
        "description": description,
        "data_filename": filename,
    }
    return meta_info


def create_meta_study(
    study_identifier: str,
    type_of_cancer: str,
    name: str,
    description: str,
    groups: str,
    short_name: str,
) -> dict:
    """Create study metadata file

    Args:
        study_identifier (str): A string used to uniquely identify this cancer study
                                within the database
        type_of_cancer (str): The cancer type abbreviation
        name (str): The name of the cancer study
        description (str): A description of the cancer study
        groups (str): When using an authenticating cBioPortal, lists the user-groups
                      that are allowed access to this study.
        short_name (str): The abbreviated name of the cancer study

    Returns:
        dict: cBioPortal meta file
    """
    meta_info = {
        "type_of_cancer": type_of_cancer,
        "cancer_study_identifier": study_identifier,
        "name": name,
        "description": description,
        "groups": groups,
        "short_name": short_name,
    }
    return meta_info


def get_cbio_file_metadata(study_identifier: str, cbio_filename: str) -> dict:
    """Get cBioPortal metadata in a dict

    Args:
        study_identifier (str): A string used to uniquely identify this cancer study
                                within the database
        cbio_filename (str): cbioportal file name

    Returns:
        dict: cBioPortal metadata
    """
    # supp_survival* files don't need a meta file
    if cbio_filename.startswith("data_clinical_sample"):
        meta_dict = create_clinical_meta_file(
            study_identifier=study_identifier,
            alteration_type="CLINICAL",
            datatype="SAMPLE_ATTRIBUTE",
            filename=cbio_filename,
        )
    elif cbio_filename.startswith("data_clinical_patient"):
        meta_dict = create_clinical_meta_file(
            study_identifier=study_identifier,
            alteration_type="CLINICAL",
            datatype="PATIENT_ATTRIBUTE",
            filename=cbio_filename,
        )
    elif cbio_filename.startswith("data_gene_matrix"):
        meta_dict = create_clinical_meta_file(
            study_identifier=study_identifier,
            alteration_type="GENE_PANEL_MATRIX",
            datatype="GENE_PANEL_MATRIX",
            filename=cbio_filename,
        )
    elif cbio_filename.startswith("data_timeline"):
        meta_dict = create_clinical_meta_file(
            study_identifier=study_identifier,
            alteration_type="CLINICAL",
            datatype="TIMELINE",
            filename=cbio_filename,
        )
    elif cbio_filename.startswith("data_CNA"):
        meta_dict = create_genomic_meta_file(
            study_identifier=study_identifier,
            alteration_type="COPY_NUMBER_ALTERATION",
            datatype="DISCRETE",
            stable_id="cna",
            profile_name="Copy-number alterations",
            profile_description="Copy-number alterations",
            filename=cbio_filename,
        )
    elif cbio_filename.startswith("data_fusions"):
        meta_dict = create_genomic_meta_file(
            study_identifier=study_identifier,
            alteration_type="FUSION",
            datatype="FUSION",
            stable_id="fusion",
            profile_name="Fusions",
            profile_description="Fusions",
            filename=cbio_filename,
        )
    elif cbio_filename.startswith("data_mutations_extended"):
        meta_dict = create_genomic_meta_file(
            study_identifier=study_identifier,
            alteration_type="MUTATIONS_EXTENDED",
            datatype="MAF",
            stable_id="mutations",
            profile_name="Mutations",
            profile_description="Mutation data from next-gen sequencing.",
            filename=cbio_filename,
        )
    elif cbio_filename.endswith("data_cna_hg19.seg"):
        meta_dict = create_seg_meta_file(
            study_identifier=study_identifier,
            reference_genome_id="hg19",
            description="Segment data for the genie study",
            filename=cbio_filename,
        )
    else:
        raise NotImplementedError(f"{cbio_filename} does not have associated metafile")
    return meta_dict


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
        "data_mutations_extended.txt",
    ],
    outdir: str = "./",
):
    """Create cbioportal metadata files

    Args:
        study_identifier (str): A string used to uniquely identify this cancer study
                                within the database
        cbio_fileformats (list): List of cbioportal file names.
        outdir (str): Directory to write metadata files to. Defaults to current dir.
    """
    # supp_survival* files don't need a meta file
    for cbio_file in cbio_fileformats:
        meta_dict = get_cbio_file_metadata(
            study_identifier=study_identifier, cbio_filename=cbio_file
        )
        write_meta_file(meta_info=meta_dict, filename=cbio_file, outdir=outdir)
