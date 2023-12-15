"""BPC redcap export PROCESSES
- timeline file
- clinical file
  OS_MONTHS = death_date_int - date_first_met_int
  OS_MONTHS_PRIMARY = death_date_int - primary_dx_date_int
  All dates are converted from days to months (days/30.4)
  Add headers
  REMOVE PATIENTS/SAMPLES THAT DON'T HAVE GENIE SAMPLE IDS
"""

import math
import os
import subprocess
from typing import List

import pandas as pd
from synapseclient import Synapse

# All cbioportal file formats written in BPC
CBIO_FILEFORMATS_ALL = [
    "data_timeline_treatment.txt",
    "data_timeline_cancer_diagnosis.txt",
    "data_timeline_pathology.txt",
    "data_timeline_sample_acquisition.txt",
    "data_timeline_performance_status.txt",
    "data_timeline_medonc.txt",
    "data_timeline_imaging.txt",
    "data_timeline_sequencing.txt",
    "data_timeline_labtest.txt",
    "data_clinical_supp_survival.txt",
    "data_clinical_supp_survival_treatment.txt",
    "data_clinical_sample.txt",
    "data_clinical_patient.txt",
    "data_mutations_extended.txt",
    "data_gene_matrix.txt",
    "data_cna_hg19.seg",
    # "data_fusions.txt",
    "data_sv.txt",
    "data_CNA.txt",
]


def get_file_data(
    syn: Synapse, mappingdf: pd.DataFrame, sampletype: str, cohort: str = "NSCLC"
) -> dict:
    """Extracts the sample, patient and timeline data frame

    Args:
        syn (Synapse): Synapse connection
        mappingdf (pd.DataFrame): Mapping dataframe
        sampletype (str): sample type label
        cohort (str, optional): cohort label. Defaults to "NSCLC".

    Returns:
        dict: dictionary with two keys ('df' and 'used') corresponding to data frame
        of data for sample type and a list of Synapse IDs used
    """

    # Group by dataset because different datasets could have the
    # same variable
    datasets = mappingdf.groupby("dataset")
    finaldf = pd.DataFrame()
    used_entities = []

    for _, df in datasets:
        # Get synapse id
        synid = df["id"].unique()[0]
        table = syn.get(synid)
        used_entities.append(f"{synid}.{table.versionNumber}")
        # obtain columns to subset df
        cols = df["code"][df["sampleType"] == sampletype]
        cols = cols.tolist()
        if "record_id" not in cols:
            cols.append("record_id")
        # Must add path_rep_number for sample and sample acquisition file
        if sampletype in ["SAMPLE", "TIMELINE-SAMPLE"]:
            cols.append("path_rep_number")
        # Must add path_proc_number to sample file
        if sampletype == "SAMPLE":
            cols.append("path_proc_number")
        # Only get specific cohort and subset cols
        tabledf = pd.read_csv(table.path, low_memory=False)
        tabledf = tabledf[tabledf["cohort"] == cohort]
        tabledf = tabledf[cols]
        # Append to final dataframe if empty
        if finaldf.empty:
            finaldf = pd.concat([finaldf, tabledf])
        else:
            # Records missing pathology reports still have to be present
            # So a left merge has to happen.  This logic also assumes that
            # The pathology report dataset mapping info isn't first.
            # This also assumes that the TIMELINE-PATHOLOGY file only
            # uses columns from the pathology-report dataset
            if df["dataset"].iloc[0] == "Pathology-report level dataset":
                finaldf = finaldf.merge(
                    tabledf,
                    on=["record_id", "path_proc_number", "path_rep_number"],
                    how="left",
                )
                del finaldf["path_rep_number"]
            else:
                finaldf = finaldf.merge(tabledf, on="record_id")

    return {"df": finaldf, "used": used_entities}


def get_synid_data(
    df_map: pd.DataFrame,
    df_file: pd.DataFrame,
    sampletype: list,
    cohort: str,
) -> List:
    """Get Synapse IDs of data files used to create the corresponding
    sample type data frames.

    Args:
        df_map (pd.DataFrame): REDCap to cBioPortal mapping data frame
        df_file (pd.DataFrame): Synapse ID to dataset label data frame
        sampletype (list): list of sample type labels
        cohort (str): cohort label

    Returns:
        List: List of used synapse ids
    """
    datasets = df_map[(df_map["sampleType"].isin(sampletype)) & (df_map[cohort])][
        "dataset"
    ].unique()
    used = df_file[df_file["dataset"].isin(datasets)]["id"]
    return list(used)


def fill_cancer_dx_start_date(finaldf: pd.DataFrame) -> pd.DataFrame:
    """Fills in cancer dx start date for those missing start dates
    with zero.

    Args:
        finaldf (pd.DataFrame): Mapped cancer diagnosis timeline dataframe

    Returns:
        pd.DataFrame: dataframe with filled START_DATEs
    """
    # Get all time0 points for all records
    time0_dates_idx = finaldf["INDEX_CANCER"] == "Yes"
    # Make all time0 points 0
    finaldf["diagnosis_int"] = finaldf["START_DATE"]
    finaldf["START_DATE"][time0_dates_idx] = 0
    # Remove unused variable
    del finaldf["diagnosis_int"]
    return finaldf


def configure_mafdf(mafdf: pd.DataFrame, keep_samples: list) -> pd.DataFrame:
    """Configures a maf dataframe

    Args:
        mafdf (pd.DataFrame): Chunk of maf dataframe
        keep_samples (list):  Samples to keep in the maf file

    Returns:
        pd.DataFrame: Configured maf dataframe
    """
    keep_mafdf = mafdf[mafdf["Tumor_Sample_Barcode"].isin(keep_samples.tolist())]
    if not keep_mafdf.empty:
        fillnas = [
            "t_depth",
            "t_ref_count",
            "t_alt_count",
            "n_depth",
            "n_ref_count",
            "n_alt_count",
        ]
        for col in fillnas:
            keep_mafdf[col].loc[keep_mafdf[col] == "."] = ""
        keep_mafdf["Validation_Status"] = ""
    return keep_mafdf


def change_days_to_years(days: int) -> float:
    """Convert days into years.

    Args:
        days (int): number of days

    Returns:
        int: number of years or nan, if applicable
    """
    if math.isnan(days):
        return float("nan")
    else:
        return math.floor(days / 365.25)


def _get_synid_dd(syn: Synapse, cohort: str, synid_table_prissmm: str) -> str:
    """Get Synapse ID of the most current PRISSMM non-PHI data dictionary for the BPC cohort.

    Args:
        syn (Synapse): Synapse connection
        cohort (str): cohort label
        synid_table_prissmm (str): Synapse ID of PRISSMM documentation table

    Returns:
        str: Synapse ID of cohort PRISSMM data dictionary
    """

    query = f"SELECT id FROM {synid_table_prissmm} WHERE cohort = '{cohort}' ORDER BY name DESC LIMIT 1"
    query_results = syn.tableQuery(query)

    synid_folder_prissmm = query_results.asDataFrame()["id"][0]

    synid_prissmm_children = syn.getChildren(synid_folder_prissmm)

    for child in synid_prissmm_children:
        if child["name"] == "Data Dictionary non-PHI":
            return child["id"]
    return None


def get_drug_mapping(
    syn: Synapse, cohort: str, synid_table_prissmm: str
) -> dict:
    """Get a mapping between drug short names and NCIT code from BPC data dictionary
    and BPC global response set for a given BPC cohort.

    Args:
        syn (Synapse): Synapse connection
        cohort (str): cohort label
        synid_table_prissmm (str): Synapse ID of PRISSMM documentation table

    Returns:
        dict: map where keys are BPC drug short names and value is the
                corresponding NCIT drug code
    """

    mapping = {}
    var_names = []

    synid_file_dd = _get_synid_dd(syn, cohort, synid_table_prissmm)

    dd = pd.read_csv(
        syn.get(synid_file_dd).path, encoding="unicode_escape", low_memory=False
    )

    for i in ["1", "2", "3", "4", "5"]:
        var_names.append("drugs_drug_" + i)
        var_names.append("drugs_drug_oth" + i)

    for var_name in var_names:
        if var_name in dd["Variable / Field Name"].unique():
            choice_str = dd[dd["Variable / Field Name"] == var_name][
                "Choices, Calculations, OR Slider Labels"
            ].values[0]
            choice_str = choice_str.replace('"', "")

            for pair in choice_str.split("|"):
                if pair.strip() != "":
                    code = pair.split(",")[0].strip()
                    value = pair.split(",")[1].strip()
                    label = value.split("(")[0].strip()
                    mapping[label] = code
    return mapping


def get_regimen_abbr(regimen: str, mapping: dict) -> str:
    """Given a BPC regimen and mapping between drug names and NCIT codes,
    return the regimen abbreviation consisting of NCIT codes.

    Args:
        regimen (str): string representing a comma delimited list of drug names in the regimen
        mapping (dict): map where keys are BPC drug short names and value is the
                corresponding NCIT drug code

    Returns:
        str: regimen abbreviation with NCIT codes
    """

    abbr = ""
    drugs = regimen.split(",")
    for drug in drugs:
        if drug == drugs[0]:
            abbr = mapping[drug.strip()]
        else:
            abbr = abbr + "_" + mapping[drug.strip()]
    return abbr


def get_git_sha() -> str:
    """get git sha digest"""
    text = subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True, text=True)
    return text.stdout.rstrip("\n")


def create_regimens(
    syn: Synapse,
    regimen_infodf: pd.DataFrame,
    mapping: dict,
    top_x_regimens: int = 5,
    cohort: str = "NSCLC",
) -> dict:
    """Create regimens to merge into the patient file.

    Args:
        syn (Synapse): Synapse connection
        regimen_infodf (pd.DataFrame): data frame containing regimen info
        mapping (dict): map where keys are BPC drug short names and value is the
                corresponding NCIT drug code
        top_x_regimens (int, optional): number of regimens to catalog. Defaults to 5.
        cohort (str, optional): cohort label. Defaults to "NSCLC".

    Returns:
        dict: dictionary with three keys ('df', 'used', 'info')
    """

    regimen_synid = regimen_infodf["id"].unique()[0]
    regimens_to_exclude = ["Investigational Drug"]
    regimen_ent = syn.get(regimen_synid)
    regimendf = pd.read_csv(regimen_ent.path, low_memory=False)
    # Get only NSCLC cohort
    regimendf = regimendf[regimendf["cohort"] == cohort]
    # Use redcap_ca_index == Yes
    regimendf = regimendf[regimendf["redcap_ca_index"] == "Yes"]
    # Exclude regimens
    regimendf = regimendf[~regimendf["regimen_drugs"].isin(regimens_to_exclude)]
    regimendf = regimendf[
        ~regimendf["regimen_drugs"].str.contains("Investigational Drug")
    ]
    # Exclude all regimens with "Other"
    regimendf = regimendf[~regimendf["regimen_drugs"].str.contains("Other")]
    # sort file by regimen_number and drop rest of duplicates
    # (not all duplicates), if duplicated keep the first regimen
    regimendf.sort_values("regimen_number", inplace=True)
    regimendf.drop_duplicates(["record_id", "regimen_drugs"], inplace=True)

    count_of_regimens = regimendf["regimen_drugs"].value_counts()
    # Obtain top X number of regimens
    to_include_regimens = count_of_regimens[:top_x_regimens].index.tolist()

    subset_regimendf = regimendf[regimendf["regimen_drugs"].isin(to_include_regimens)]
    regimen_groups = subset_regimendf.groupby("regimen_drugs")
    new_regimen_info = pd.DataFrame()
    # Create regimen clinical headers
    final_regimendf = pd.DataFrame()
    for regimen, df in regimen_groups:
        regimen_drug_info = regimen_infodf.copy()
        # Create regimen drug abbreviations
        regimen_abbr = get_regimen_abbr(regimen, mapping)

        # Create correct column mappings for the clinical patient file
        regimen_drug_info["cbio"] = [
            value.format(regimen_abbr=regimen_abbr) for value in regimen_infodf["cbio"]
        ]
        regimen_drug_info["labels"] = [
            value.format(regimen=regimen) for value in regimen_infodf["labels"]
        ]
        regimen_drug_info["description"] = [
            value.format(regimen=regimen) for value in regimen_infodf["description"]
        ]
        regimen_drug_info["priority"] = [
            int(value) for value in regimen_infodf["priority"]
        ]
        new_regimen_info = pd.concat([new_regimen_info, regimen_drug_info])

        col_map = regimen_drug_info["cbio"].to_dict()
        col_map["record_id"] = "PATIENT_ID"
        regimen_patientdf = df[list(col_map.keys())].rename(columns=col_map)
        # Merge final regimen dataframe
        if final_regimendf.empty:
            final_regimendf = regimen_patientdf
        else:
            final_regimendf = final_regimendf.merge(
                regimen_patientdf, on="PATIENT_ID", how="outer"
            )
    return {"df": final_regimendf, "info": new_regimen_info, "used": regimen_synid}


def get_bpc_to_cbio_mapping_df(
    syn: Synapse, cohort: str, synid_table_cbio: str
) -> pd.DataFrame:
    """Extract relevant portions of the mapping table for the sponsored project
    variables and cBioPortal variable mappings and return as a data frame.

    Args:
        syn (Synapse): Synapse connection
        cohort (str): sponsored project label
        synid_table_cbio (str): Synapse ID of table containing variable to cbio mapping

    Returns:
        pd.DataFrame: data frame of all mapping columns for released variables
    """
    redcap_to_cbiomapping = syn.tableQuery(
        f"SELECT * FROM {synid_table_cbio} where "
        f"{cohort} is true AND sampleType <> 'TIMELINE-STATUS'"
    )
    redcap_to_cbiomappingdf = redcap_to_cbiomapping.asDataFrame()
    return redcap_to_cbiomappingdf


def get_data_file_synapse_id_df(syn: Synapse, synid_table_files: str) -> pd.DataFrame:
    """Get mapping of dataset files to Synapse IDs.

    Args:
        syn (Synapse): Synapse connection
        synid_table_files (str): Synapse ID of table containing data file to Synapse ID mapping

    Returns:
        pd.DataFrame: data frame with two columns representing the Synapse ID and dataset label
    """
    data_tables = syn.tableQuery(f"SELECT id, dataset FROM {synid_table_files} ")
    data_tablesdf = data_tables.asDataFrame()
    return data_tablesdf


def create_release_folders(cohort: str) -> None:
    """Create local folders for release folders.

    Args:
        cohort (str): sponsored project label
    """
    if not os.path.exists(cohort):
        os.mkdir(cohort)
    else:
        filelists = os.listdir(cohort)
        for each_file in filelists:
            if each_file != "case_lists":
                os.remove(os.path.join(cohort, each_file))


def hack_remap_laterality(df_patient_subset: pd.DataFrame) -> pd.DataFrame:
    """Temporary (?) remapping of specific values in a column.

    Args:
        df_patient_subset (pd.DataFrame): patient data

    Returns:
        pd.DataFrame: patient data wtih column remapped
    """

    laterality_mapping = {
        "0": "Not a paired site",
        "1": "Right: origin of primary",
        "2": "Left: origin of primary",
        "3": "Only one side involved, right or left origin unspecified",
        "4": "Bilateral involvement at time of diagnosis, lateral origin "
        "unknown for a single primary; or both ovaries involved "
        "simultaneously, single histology; bilateral retinoblastomas; "
        "bilateral Wilms' tumors",
        "5": "Paired site: midline tumor",
        "9": "Paired site, but no information concerning laterality",
        "Not paired": "Not a paired site",
    }
    if df_patient_subset.get("NAACCR_LATERALITY_CD") is not None:
        remapped_values = (
            df_patient_subset["NAACCR_LATERALITY_CD"]
            .astype(str)
            .map(laterality_mapping)
        )
        df_patient_subset["NAACCR_LATERALITY_CD"] = remapped_values

    return df_patient_subset


def remap_os_values(df: pd.DataFrame):
    """Remap OS numerical values to string values
    0 -> 0:LIVING
    1 -> 1:DECEASED
    """
    remap_values = {col: {0: "0:LIVING", 1: "1:DECEASED"}
                    for col in df.columns
                    if col.startswith('OS') and col.endswith("STATUS")}
    return df.replace(remap_values)


def remap_pfs_values(df: pd.DataFrame):
    """Remap PFS numerical values to string values
    0 -> 0:CENSORED
    1 -> 1:PROGRESSED
    """
    remap_values = {col: {0: "0:CENSORED", 1: "1:PROGRESSED"}
                    for col in df.columns
                    if col.startswith('PFS') and col.endswith("STATUS")}
    return df.replace(remap_values)


def _convert_to_int(value):
    """Convert object to integer or return nan"""
    try:
        return int(value)
    except ValueError:
        return float('nan')


def get_survival_info(syn, df_map, df_file, cohort, prissm_synid):
    # Patient and Sample mapping values

    patient_sample_idx = df_map["sampleType"].isin(
        ["PATIENT", "SAMPLE", "SURVIVAL"]
    )
    infodf = df_map[patient_sample_idx].merge(df_file, on="dataset", how="left")
    infodf.index = infodf["code"]

    # Regimen mapping values
    regimen_idx = df_map["sampleType"].isin(["REGIMEN"])
    regimen_infodf = df_map[regimen_idx].merge(df_file, on="dataset", how="left")
    regimen_infodf.index = regimen_infodf["code"]

    # Create regimens data for patient file
    drug_mapping = get_drug_mapping(
        syn=syn,
        cohort=cohort,
        synid_table_prissmm=prissm_synid
    )
    regimens_data = create_regimens(
        syn,
        regimen_infodf,
        mapping=drug_mapping,
        top_x_regimens=20,
        cohort=cohort
    )
    survival_info = pd.concat([infodf, regimens_data["info"]])
    return survival_info
