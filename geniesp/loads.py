import pandas as pd

import synapseclient
from genie import process_functions


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


# TODO: need to have a load class
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
    # HACK this is to remove the index cancer for the patient file
    # redcap_to_cbiomappingdf = redcap_to_cbiomappingdf[redcap_to_cbiomappingdf['cbio'] != "INDEX_CANCER"]
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
        if "SURVIVAL" in clin_path or "REGIMEN" in clin_path:
            clin_file.write("#{}\n".format("\t".join(["PATIENT"] * len(labels))))
        clin_file.write("#{}\n".format("\t".join(priority)))
        clin_file.write(
            process_functions.removeStringFloat(
                clinicaldf.to_csv(index=False, sep="\t")
            )
        )
    return clin_path
