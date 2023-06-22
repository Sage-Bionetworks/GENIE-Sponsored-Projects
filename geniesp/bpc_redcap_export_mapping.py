"""BPC redcap export PROCESSES
- timeline file
- clinical file
  OS_MONTHS = death_date_int - date_first_met_int
  OS_MONTHS_PRIMARY = death_date_int - primary_dx_date_int
  All dates are converted from days to months (days/30.4)
  Add headers
  REMOVE PATIENTS/SAMPLES THAT DON'T HAVE GENIE SAMPLE IDS
"""
from abc import ABCMeta
from datetime import date
from functools import cached_property
import math
import os
import subprocess
import logging
from typing import List

from genie import create_case_lists, process_functions
import pandas as pd
from synapseclient import File, Folder, Synapse

from . import metafiles

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
    "data_fusions.txt",
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


def change_days_to_years(days: int) -> int:
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
    syn: Synapse, cohort: str, synid_file_grs: str, synid_table_prissmm: str
) -> dict:
    """Get a mapping between drug short names and NCIT code from BPC data dictionary
    and BPC global response set for a given BPC cohort.

    Args:
        syn (Synapse): Synapse connection
        cohort (str): cohort label
        synid_file_grs (str): Synapse ID of REDCap global response set file.
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
    grs = pd.read_csv(
        syn.get(synid_file_grs).path, encoding="unicode_escape", low_memory=False
    )
    grs.columns = ["Variable / Field Name", "Choices, Calculations, OR Slider Labels"]

    for i in ["1", "2", "3", "4", "5"]:
        var_names.append("drugs_drug_" + i)
        var_names.append("drugs_drug_oth" + i)

    for obj in dd, grs:

        for var_name in var_names:

            if var_name in obj["Variable / Field Name"].unique():
                choice_str = obj[obj["Variable / Field Name"] == var_name][
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


class BpcProjectRunner(metaclass=ABCMeta):
    """BPC redcap to cbioportal export"""

    # Sponsored project name
    _SPONSORED_PROJECT = ""
    # Redcap codes to cbioportal mapping synid and form key is in
    # version 38 was the last stable version
    # Use version 42 - but there is a bug in Synapse...
    _REDCAP_TO_CBIOMAPPING_SYNID = "syn25712693.42"
    # Run `git rev-parse HEAD` in Genie_processing directory to obtain shadigest
    _GITHUB_REPO = None
    # Mapping from Synapse Table to derived variables
    # TODO: Make versioned
    _DATA_TABLE_IDS = "syn22296821"
    # Storage of not found samples
    _SP_REDCAP_EXPORTS_SYNID = "syn21446571"
    # main GENIE release folder (12.0-public)
    _MG_RELEASE_SYNID = "syn32309524"
    # PRISSMM documentation table
    _PRISSMM_SYNID = "syn22684834"
    # REDCap global response set
    _GRS_SYNID = "syn24184523"
    # main GENIE sample clinical database
    _CLINICAL_SYNID = "syn7517674"
    # BPC sample retraction table
    _RETRACTION_SYNID = "syn25779833"
    # main GENIE assay information table
    _ASSAY_SYNID = "syn17009222"
    # exclude files to be created for cbioportal
    # TODO: need to support this feature in rest of code, for now
    # This is added for metadata files
    _exclude_files = []
    # cohort-generic link to documentation for BPC datasets
    _url_bpc = "https://aacr.box.com/s/en5dyu9zfw1krg2u58wlcz01jttc6y9h"
    # cohort-generic link to documentation for cBio files
    _url_cbio = "https://docs.google.com/document/d/1IBVF-FLecUG8Od6mSEhYfWH3wATLNMnZcBw2_G0jSAo/edit"

    def __init__(self, syn, cbiopath, release, upload=False):
        if not os.path.exists(cbiopath):
            raise ValueError("cbiopath doesn't exist")
        if self._SPONSORED_PROJECT == "":
            raise ValueError("Must configure _SPONSORED_PROJECT")
        self.syn = syn
        self.cbiopath = cbiopath
        self.upload = upload
        self.release = release
        self._GITHUB_REPO = f"https://github.com/Sage-Bionetworks/GENIE-Sponsored-Projects/tree/{get_git_sha()}"

    @cached_property
    def genie_clinicaldf(self) -> pd.DataFrame:
        """Get main GENIE clinical samples and perform retraction
        against BPC sample and patient database. It's important that I use the
        sample database, because if a sample from a patient with multiple samples
        is retracted from the main GENIE clinical samples, the patient
        will still exist. (Jira: GEN-260)

        Returns:
            pd.DataFrame: main GENIE information for sponsored project samples
        """
        genie_clinicaldb = self.syn.tableQuery(
            f"select SAMPLE_ID, PATIENT_ID, ONCOTREE_CODE, SEQ_ASSAY_ID, "
            f"SAMPLE_TYPE, SEQ_YEAR from {self._CLINICAL_SYNID}"
        )
        genie_clinicaldf = genie_clinicaldb.asDataFrame()
        # BPC retraction database
        bpc_retraction_db = self.syn.tableQuery(
            f"select SAMPLE_ID from {self._RETRACTION_SYNID} where "
            f"{self._SPONSORED_PROJECT} is true"
        )
        bpc_retractiondf = bpc_retraction_db.asDataFrame()

        # TODO: add patient retraction database...
        bpc_patient_retraction_db = self.syn.tableQuery(
            "select record_id from syn25998970 where "
            f"{self._SPONSORED_PROJECT} is true"
        )
        bpc_patient_retraction_df = bpc_patient_retraction_db.asDataFrame()

        bpc_temp_patient_retraction_db = self.syn.tableQuery(
            "select record_id from syn29266682 where "
            f"cohort = '{self._SPONSORED_PROJECT}'"
        )
        bpc_temp_patient_retraction_df = bpc_temp_patient_retraction_db.asDataFrame()
        # Retract samples from sample retraction db
        keep_clinicaldf = genie_clinicaldf[
            ~genie_clinicaldf["SAMPLE_ID"].isin(bpc_retractiondf["SAMPLE_ID"])
        ]
        # Retract patients from patient retraction db
        keep_clinicaldf = keep_clinicaldf[
            ~keep_clinicaldf["PATIENT_ID"].isin(bpc_patient_retraction_df["record_id"])
        ]
        # Retract patients from temporary patient retraction db
        keep_clinicaldf = keep_clinicaldf[
            ~keep_clinicaldf["PATIENT_ID"].isin(
                bpc_temp_patient_retraction_df["record_id"]
            )
        ]
        return keep_clinicaldf

    @cached_property
    def cbioportal_folders(self) -> dict:
        """Create case lists and release folder"""
        # parent_id = "syn50876969" if self.upload else "syn21241322"
        if self.upload:
            parent_id = "syn50876969"
            sp_data_folder = self.syn.store(
                Folder(self._SPONSORED_PROJECT, parentId=parent_id)
            )
            release_folder = self.syn.store(Folder(self.release, parent=sp_data_folder))
            # if not self.upload:
            #     release_folder = self.syn.store(
            #         Folder("cBioPortal_files", parent=release_folder)
            #     ).id
            case_lists = self.syn.store(Folder("case_lists", parent=release_folder))
            return {"release": release_folder, "case_lists": case_lists}
        return {}

    def create_bpc_cbio_metafiles(self) -> List:
        """Create BPC cBioPortal meta* files.

        Returns:
            List: file paths of meta files.
        """

        mg_release_ent = self.syn.get(self._MG_RELEASE_SYNID)
        name = f"GENIE BPC {self._SPONSORED_PROJECT} v{self.release}"
        description = (
            f"{self._SPONSORED_PROJECT} cohort v{self.release} "
            f"(GENIE {date.today().year}) GENIE {mg_release_ent.name}. "
            f"Several hundred different variables are collected for each of "
            f'the BPC cohorts; consult the <a href="{self._url_bpc}">Documentation</a> '
            f"for further detail. To learn more about which variables are "
            f"visualized in cBioPortal and how, see the cBioPortal "
            f'<a href="{self._url_cbio}">ReadMe</a>. '
            '<font color="red">Although these data are de-identified, your analysis '
            "may require institutional review prior to publication.</font>"
        )
        short_name = f"{self._SPONSORED_PROJECT} GENIE"
        study_identifier = f"{self._SPONSORED_PROJECT.lower()}_genie_bpc"
        # Get list of files to create cBioPortal metadata files for
        to_create_meta = list(set(CBIO_FILEFORMATS_ALL) - set(self._exclude_files))

        meta_files = metafiles.create_cbio_metafiles(
            study_identifier=study_identifier,
            outdir=self._SPONSORED_PROJECT,
            cbio_fileformats=to_create_meta,
        )
        meta_study = metafiles.create_study_meta_file(
            study_identifier=study_identifier,
            type_of_cancer="mixed",
            name=name,
            description=description,
            groups="GENIE",
            short_name=short_name,
        )
        study_file = metafiles.write_meta_file(
            meta_info=meta_study,
            filename="meta_study.txt",
            outdir=self._SPONSORED_PROJECT,
        )
        meta_files.append(study_file)
        return meta_files

    def create_and_write_genematrix(
        self, clinicaldf: pd.DataFrame, cna_samples: list, used_ent: list = None
    ) -> None:
        """Create gene matrix dataframe.

        Args:
            clinicaldf (pd.DataFrame): sample information
            cna_samples (list): sample IDs with CNA data
            used_ent (list, optional): Synapse IDs used. Defaults to None.
        """

        data_gene_panel = clinicaldf[["SAMPLE_ID", "SEQ_ASSAY_ID"]]
        data_gene_panel = data_gene_panel.rename(columns={"SEQ_ASSAY_ID": "mutations"})
        data_gene_panel = data_gene_panel[data_gene_panel["SAMPLE_ID"] != ""]
        data_gene_panel.drop_duplicates("SAMPLE_ID", inplace=True)
        cna_seqids = data_gene_panel["mutations"][
            data_gene_panel["SAMPLE_ID"].isin(cna_samples)
        ].unique()
        data_gene_panel["cna"] = data_gene_panel["mutations"]
        data_gene_panel["cna"][~data_gene_panel["cna"].isin(cna_seqids)] = "NA"
        data_gene_panel.fillna("NA", inplace=True)
        gene_matrix_filepath = os.path.join(
            self._SPONSORED_PROJECT, "data_gene_matrix.txt"
        )
        data_gene_panel.to_csv(gene_matrix_filepath, sep="\t", index=False)
        if self.upload:
            file_ent = File(
                gene_matrix_filepath, parent=self.cbioportal_folders["release"]
            )
            self.syn.store(file_ent, used=used_ent, executed=self._GITHUB_REPO)

    def configure_clinicaldf(
        self, clinicaldf: pd.DataFrame, redcap_to_cbiomappingdf: pd.DataFrame
    ) -> pd.DataFrame:
        """Create clinical file from sponsored project mapped dataframe

        Args:
            clinicaldf (pd.DataFrame): clinical information
            redcap_to_cbiomappingdf (pd.DataFrame): cBio mapping info

        Raises:
            ValueError: All column names must be in mapping dataframe
            ValueError: Must have no null patient ids
            ValueError: Must have no null sample ids

        Returns:
            pd.DataFrame: configured clinical information
        """

        if not clinicaldf.columns.isin(redcap_to_cbiomappingdf["code"]).all():
            raise ValueError("All column names must be in mapping dataframe")
        mapping = redcap_to_cbiomappingdf["cbio"].to_dict()
        clinicaldf = clinicaldf.rename(columns=mapping)
        clinicaldf = clinicaldf.drop_duplicates()
        if sum(clinicaldf["PATIENT_ID"].isnull()) > 0:
            raise ValueError("Must have no null patient ids")
        # Remove white spaces for PATIENT/SAMPLE ID
        clinicaldf["PATIENT_ID"] = [
            patient.strip() for patient in clinicaldf["PATIENT_ID"]
        ]

        if clinicaldf.get("SAMPLE_ID") is not None:
            # This line should not be here
            clinicaldf = clinicaldf[~clinicaldf["SAMPLE_ID"].isnull()]
            if sum(clinicaldf["SAMPLE_ID"].isnull()) > 0:
                raise ValueError("Must have no null sample ids")
            clinicaldf["SAMPLE_ID"] = [
                sample.strip() for sample in clinicaldf["SAMPLE_ID"]
            ]
        # JIRA: GEN-10- Must standardize SEQ_ASSAY_ID values to be uppercase
        if clinicaldf.get("SEQ_ASSAY_ID") is not None:
            clinicaldf["SEQ_ASSAY_ID"] = clinicaldf["SEQ_ASSAY_ID"].str.upper()

        clinicaldf["SP"] = self._SPONSORED_PROJECT

        for col in clinicaldf:
            num_missing = sum(clinicaldf[col].isnull())
            if num_missing > 0:
                logging.warning(f"Number of missing {col}: {num_missing}")

        return clinicaldf

    def write_clinical_file(
        self,
        clinicaldf: pd.DataFrame,
        redcap_to_cbiomappingdf: pd.DataFrame,
        filetype: str,
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

        if filetype not in [
            "patient",
            "sample",
            "supp_survival",
            "supp_survival_treatment",
        ]:
            raise ValueError(
                "sample type must be patient, sample, supp_survival or "
                "supp_survival_treatment"
            )
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

        clin_path = os.path.join(
            self._SPONSORED_PROJECT, f"data_clinical_{filetype}.txt"
        )

        with open(clin_path, "w+") as clin_file:
            clin_file.write("#{}\n".format("\t".join(labels)))
            clin_file.write("#{}\n".format("\t".join(descriptions)))
            clin_file.write("#{}\n".format("\t".join(coltype)))
            # attributes in the supp file are PATIENT, so must
            # specify that.
            if filetype.startswith("supp_survival"):
                clin_file.write("#{}\n".format("\t".join(["PATIENT"] * len(labels))))
            clin_file.write("#{}\n".format("\t".join(priority)))

            clin_file.write(
                process_functions.removeStringFloat(
                    clinicaldf.to_csv(index=False, sep="\t")
                )
            )
        return clin_path

    def make_timeline_treatmentdf(self, infodf: pd.DataFrame, sample_type: str) -> dict:
        """Make timeline treatment data frame.

        Args:
            infodf (pd.DataFrame): cBioPortal mapping merged with data file Synapse IDs
            sample_type (str): sample type label

        Returns:
            dict: dictionary with two keys ('df', 'used')
        """

        subset_infodf = infodf[infodf["sampleType"] == sample_type]
        # Exclude heme onc columns
        subset_infodf = subset_infodf[
            ~subset_infodf["data_type"].isin(["portal_value", "heme"])
        ]
        synid = subset_infodf["id"].unique()[0]
        ent = self.syn.get(synid)
        used_entity = f"{synid}.{ent.versionNumber}"
        timelinedf = pd.read_csv(ent.path, low_memory=False)
        # Only take lung cohort
        timelinedf = timelinedf[timelinedf["cohort"] == self._SPONSORED_PROJECT]
        # Only take samples where redcap_ca_index is Yes
        timelinedf = timelinedf[timelinedf["redcap_ca_index"] == "Yes"]
        # Flatten multiple columns values into multiple rows
        multiple_cols_idx = subset_infodf["code"].str.contains("[*]")
        final_timelinedf = pd.DataFrame()
        for _, row in subset_infodf[multiple_cols_idx].iterrows():
            code = row["code"]
            # Make sure to use regex to get colname + integer
            new_code = code.replace("*", "[\d]")
            cols = timelinedf.columns.str.contains(new_code)
            wanted_cols = timelinedf.columns[cols].tolist()
            wanted_cols.extend(["record_id", "regimen_drugs", "regimen_number"])
            # melt function creates multiple rows from multiple columns
            melted_df = pd.melt(
                timelinedf[wanted_cols],
                id_vars=["record_id", "regimen_drugs", "regimen_number"],
                value_name=row["cbio"],
            )
            del melted_df["variable"]
            if final_timelinedf.empty:
                final_timelinedf = pd.concat([final_timelinedf, melted_df])
            else:
                final_timelinedf[row["cbio"]] = melted_df[row["cbio"]]
        final_timelinedf["TREATMENT_TYPE"] = "Systemic Therapy"
        # Remove all START_DATE is NULL
        final_timelinedf = final_timelinedf[~final_timelinedf["START_DATE"].isnull()]
        non_multi_cols = subset_infodf[~multiple_cols_idx]["code"].tolist()
        non_multi_cols.append("record_id")

        # Merge final timeline
        final_timelinedf = final_timelinedf.merge(
            timelinedf[non_multi_cols],
            on=["record_id", "regimen_drugs", "regimen_number"],
        )

        # Make sure all events types are treatment
        final_timelinedf["EVENT_TYPE"] = "Treatment"

        # Make sure AGENT is not null and doesn't have parenthesis
        agents = []
        final_timelinedf = final_timelinedf[~final_timelinedf["AGENT"].isnull()]
        for index, agent in enumerate(final_timelinedf["AGENT"]):
            if "(" in agent:
                agents.append(agent.split("(")[0].strip())
            else:
                agents.append(agent.split(",")[0].strip())

        final_timelinedf["AGENT"] = agents
        # Map timeline treatment columns
        mapping = subset_infodf["cbio"].to_dict()
        # Must add in PATIENT_ID
        mapping["record_id"] = "PATIENT_ID"
        final_timelinedf = final_timelinedf.rename(columns=mapping)
        # timeline file must be in this order
        cols_to_order = [
            "PATIENT_ID",
            "START_DATE",
            "STOP_DATE",
            "EVENT_TYPE",
            "TREATMENT_TYPE",
            "AGENT",
        ]
        cols_to_order.extend(final_timelinedf.columns.drop(cols_to_order).tolist())
        final_timelinedf = self.filter_df(final_timelinedf)

        return {
            "df": final_timelinedf[cols_to_order].drop_duplicates(),
            "used": [used_entity],
        }

    def filter_df(self, df: pd.DataFrame) -> pd.DataFrame:
        """Make sure samples and patients exist in the main genie
        clinical samples with the exception of the removal of
        retraction database samples

        Args:
            df (pd.DataFrame): GENIE data with a patient or sample ID column

        Returns:
            pd.DataFrame: GENIE data with retracted patients or samples removed.
        """

        if df.get("SAMPLE_ID") is not None:
            to_keep_samples_idx = df["SAMPLE_ID"].isin(
                self.genie_clinicaldf["SAMPLE_ID"]
            )
            df = df[to_keep_samples_idx]
        elif df.get("PATIENT_ID") is not None:
            to_keep_patient_idx = df["PATIENT_ID"].isin(
                self.genie_clinicaldf["PATIENT_ID"]
            )
            df = df[to_keep_patient_idx]
        return df

    def write_and_storedf(
        self, df: pd.DataFrame, filepath: str, used_entities: list = []
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

        if self.upload:
            # Add the mapping file to the release file provenance
            used_entities.append(self._REDCAP_TO_CBIOMAPPING_SYNID)
            ent = File(filepath, parent=self.cbioportal_folders["release"])
            self.syn.store(ent, executed=self._GITHUB_REPO, used=used_entities)

    def create_fixed_timeline_files(
        self,
        timeline_infodf: pd.DataFrame,
        timeline_type: str,
        filter_start: bool = True,
    ) -> dict:
        """Create timeline files straight from derived variables.

        Args:
            timeline_infodf (pd.DataFrame): cBio mapping information relevant to the timeline
            timeline_type (str): timeline label
            filter_start (bool, optional): whether to filter out rows with null START_DATEs. Defaults to True.

        Returns:
            dict: 'df' mapped dataframe, 'used' list of entities
        """
        # Remove portal values
        subset_infodf = timeline_infodf[timeline_infodf["sampleType"] == timeline_type]
        # Obtain portal value (EVENT_TYPE)
        portal_value_idx = subset_infodf["data_type"] == "portal_value"
        # HACK: Passing in data mapping without portal_values should
        # still generate a file.
        # if sum(portal_value_idx) > 0:
        portal_value = subset_infodf["code"][portal_value_idx].values[0]
        # else:
        #    portal_value = ""
        subset_infodf = subset_infodf[subset_infodf["data_type"] != "portal_value"]
        timeline_data = get_file_data(
            self.syn, subset_infodf, timeline_type, cohort=self._SPONSORED_PROJECT
        )
        timelinedf = timeline_data["df"]
        used_entities = timeline_data["used"]

        timelinedf["EVENT_TYPE"] = portal_value
        mapping = subset_infodf["cbio"].to_dict()
        # Must add in PATIENT_ID
        mapping["record_id"] = "PATIENT_ID"
        timelinedf = timelinedf.rename(columns=mapping)
        timelinedf["STOP_DATE"] = ""
        # timeline file must be in this order
        cols_to_order = ["PATIENT_ID", "START_DATE", "STOP_DATE", "EVENT_TYPE"]
        cols_to_order.extend(timelinedf.columns.drop(cols_to_order).tolist())
        timelinedf = self.filter_df(timelinedf)
        # Remove all null START_DATE rows if requested
        if filter_start:
            timelinedf = timelinedf[~timelinedf["START_DATE"].isnull()]
        return {
            "df": timelinedf[cols_to_order].drop_duplicates(),
            "used": used_entities,
        }

    def get_mg_synid(self, synid_folder: str, file_name: str) -> str:
        """Get Synapse ID of main GENIE data file in release folder.

        Args:
            synid_folder (str): Synapse ID of main GENIE release folder
            file_name (str): File name for which to retrieve Synapse ID

        Returns:
            str: Synapse ID if file found; otherwise, None
        """
        synid_children = self.syn.getChildren(synid_folder)
        for synid_child in synid_children:
            if synid_child["name"] == file_name:
                return synid_child["id"]
        raise ValueError(f"file '{file_name}' not found in {synid_folder}")

    def create_and_write_maf(self, keep_samples: list) -> str:
        """Create maf file from release maf

        Args:
            keep_samples (list): List of samples to keep

        Returns:
            str: file path to written data
        """
        file_name = "data_mutations_extended.txt"
        mafpath = os.path.join(self._SPONSORED_PROJECT, file_name)
        maf_synid = self.get_mg_synid(self._MG_RELEASE_SYNID, file_name)
        maf_ent = self.syn.get(maf_synid)
        maf_chunks = pd.read_table(maf_ent.path, chunksize=50000, low_memory=False)
        index = 0
        for maf_chunk in maf_chunks:
            mafdf = configure_mafdf(maf_chunk, keep_samples)
            # Skip to next chunk if empty
            if mafdf.empty:
                continue
            # If maf file has not been created
            if index == 0:
                maf_text = process_functions.removePandasDfFloat(mafdf)
                with open(mafpath, "w") as maf_f:
                    maf_f.write(maf_text)
            else:
                maf_text = mafdf.to_csv(sep="\t", header=None, index=False)
                maf_text = process_functions.removeStringFloat(maf_text)
                with open(mafpath, "a") as maf_f:
                    maf_f.write(maf_text)
            index += 1
        if self.upload:
            file_ent = File(mafpath, parent=self.cbioportal_folders["release"])
            self.syn.store(file_ent, used=[maf_synid], executed=self._GITHUB_REPO)

        return mafpath

    def create_and_write_cna(self, keep_samples: list) -> dict:
        """Create CNA file

        Args:
            keep_samples (list): List of samples to keep

        Returns:
            dict: "filepath" with path to written file
                    "cna_sample" list of CNA sample IDs
        """
        file_name = "data_CNA.txt"
        cna_synid = self.get_mg_synid(self._MG_RELEASE_SYNID, file_name)
        cna_path = os.path.join(self._SPONSORED_PROJECT, file_name)
        cna_ent = self.syn.get(cna_synid)
        cnadf = pd.read_table(cna_ent.path, low_memory=False)
        keep_cols = ["Hugo_Symbol"]
        keep_cols.extend(cnadf.columns[cnadf.columns.isin(keep_samples)].tolist())
        cnadf = cnadf[keep_cols]
        cna_text = process_functions.removePandasDfFloat(cnadf)
        # Must do this replace twice because \t\t\t ->
        # \tNA\t\t -> \tNA\tNA\t
        cna_text = (
            cna_text.replace("\t\t", "\tNA\t")
            .replace("\t\t", "\tNA\t")
            .replace("\t\n", "\tNA\n")
        )

        with open(cna_path, "w") as cna_file:
            cna_file.write(cna_text)

        if self.upload:
            file_ent = File(cna_path, parent=self.cbioportal_folders["release"])
            self.syn.store(file_ent, used=[cna_synid], executed=self._GITHUB_REPO)
        return {"filepath": cna_file, "cna_samples": cnadf.columns.tolist()}

    def create_and_write_fusion(self, keep_samples: list) -> str:
        """Create fusion file

        Args:
            keep_samples (list): List of samples to keep

        Returns:
            str: file path to written data
        """
        file_name = "data_fusions.txt"
        fusion_synid = self.get_mg_synid(self._MG_RELEASE_SYNID, file_name)
        fusion_ent = self.syn.get(fusion_synid)
        fusiondf = pd.read_table(fusion_ent.path, low_memory=False)
        fusiondf = fusiondf[fusiondf["Tumor_Sample_Barcode"].isin(keep_samples)]
        # cBioPortal validation fails when Hugo Symbol is null
        fusiondf = fusiondf[~fusiondf["Hugo_Symbol"].isnull()]
        fusion_path = os.path.join(self._SPONSORED_PROJECT, file_name)
        self.write_and_storedf(fusiondf, fusion_path, used_entities=[fusion_synid])
        return fusion_path

    def create_and_write_seg(self, keep_samples: list) -> str:
        """Create seg file

        Args:
            keep_samples (list): List of samples to keep

        Returns:
            str: file path to written data
        """
        # TODO: the seg filename will change 13.X release.
        file_name = "genie_data_cna_hg19.seg"
        seg_synid = self.get_mg_synid(self._MG_RELEASE_SYNID, file_name)
        seg_ent = self.syn.get(seg_synid)
        segdf = pd.read_table(seg_ent.path, low_memory=False)
        segdf = segdf[segdf["ID"].isin(keep_samples)]
        seg_path = os.path.join(self._SPONSORED_PROJECT, "data_cna_hg19.seg")
        self.write_and_storedf(segdf, seg_path, used_entities=[seg_synid])
        return seg_path

    def create_and_write_sv(self, keep_samples):
        """Create sv file

        Args:
            keep_samples: List of samples to keep
        """
        file_name = "data_sv.txt"
        # TODO: remove try except after using main genie release >= 13
        try:
            sv_synid = self.get_mg_synid(self._MG_RELEASE_SYNID, file_name)
        except ValueError:
            sv_synid = None
            logging.warning(
                f"data_sv.txt doesn't exist in main genie release: {self._MG_RELEASE_SYNID}"
            )
        if sv_synid is not None:
            sv_ent = self.syn.get(sv_synid)
            svdf = pd.read_table(sv_ent.path, low_memory=False)
            svdf = svdf[svdf["Sample_ID"].isin(keep_samples)]
            sv_path = os.path.join(self._SPONSORED_PROJECT, "data_sv.txt")
            self.write_and_storedf(svdf, sv_path, used_entities=[sv_synid])

    def create_and_write_gene_panels(self, keep_seq_assay_ids: list) -> List:
        """Create gene panels.

        Args:
            keep_seq_assay_ids (list): list of sequence assay IDs

        Returns:
            list: file paths of written data
        """

        gene_panel_paths = []
        file_name = "genomic_information.txt"
        genomic_info_synid = self.get_mg_synid(self._MG_RELEASE_SYNID, file_name)
        genomic_info_ent = self.syn.get(genomic_info_synid)
        genomic_infodf = pd.read_table(genomic_info_ent.path, low_memory=False)
        # Filter by SEQ_ASSAY_ID and only exonic regions
        genomic_infodf = genomic_infodf[
            (genomic_infodf["SEQ_ASSAY_ID"].isin(keep_seq_assay_ids))
            & (genomic_infodf["Feature_Type"] == "exon")
            & (~genomic_infodf["Hugo_Symbol"].isnull())
            & (genomic_infodf["includeInPanel"])
        ]

        seq_assay_groups = genomic_infodf.groupby("SEQ_ASSAY_ID")
        for seq_assay_id, seqdf in seq_assay_groups:
            unique_genes = seqdf.Hugo_Symbol.unique()
            gene_panel_text = (
                "stable_id: {seq_assay_id}\n"
                "description: {seq_assay_id}, "
                "Number of Genes - {num_genes}\n"
                "gene_list:\t{genelist}".format(
                    seq_assay_id=seq_assay_id,
                    num_genes=len(unique_genes),
                    genelist="\t".join(unique_genes),
                )
            )
            gene_panel_name = f"data_gene_panel_{seq_assay_id}.txt"
            gene_panel_path = os.path.join(self._SPONSORED_PROJECT, gene_panel_name)
            gene_panel_paths.append(gene_panel_path)
            with open(gene_panel_path, "w+") as f:
                f.write(gene_panel_text)
            if self.upload:
                fileEnt = File(
                    gene_panel_path, parent=self.cbioportal_folders["release"]
                )
                self.syn.store(
                    fileEnt, used=[genomic_info_synid], executed=self._GITHUB_REPO
                )
        return gene_panel_paths

    def get_timeline_treatment(
        self, df_map: pd.DataFrame, df_file: pd.DataFrame
    ) -> dict:
        """Get TIMELINE-TREATMENT file data.

        Args:
            df_map (pd.DataFrame): variable to cBioPortal mapping info
            df_file (pd.DataFrame): data file to Synapse ID mapping

        Returns:
            dict: TIMELINE-TREATMENT data
        """

        timeline_infodf = df_map.query('sampleType == "TIMELINE-TREATMENT"').merge(
            df_file, on="dataset", how="left"
        )
        timeline_infodf.index = timeline_infodf["code"]
        treatment_data = self.make_timeline_treatmentdf(
            timeline_infodf, "TIMELINE-TREATMENT"
        )

        return treatment_data

    def get_timeline_performance(
        self, df_map: pd.DataFrame, df_file: pd.DataFrame
    ) -> dict:
        """Get TIMELINE-PERFORMANCE file data.

        Args:
            df_map (pd.DataFrame): variable to cBioPortal mapping info
            df_file (pd.DataFrame): data file to Synapse ID mapping

        Returns:
            dict: TIMELINE-PERFORMANCE data
        """

        timeline_infodf = df_map.query('sampleType == "TIMELINE-PERFORMANCE"').merge(
            df_file, on="dataset", how="left"
        )
        timeline_infodf.index = timeline_infodf["code"]
        data = self.create_fixed_timeline_files(timeline_infodf, "TIMELINE-PERFORMANCE")
        # HACK: Due to remapping logic, we will re-create RESULT column with correct
        # values
        data['MD_KARNOF'] = data['RESULT']
        data['RESULT'] = [val.split(":")[0] if pd.isnull(val) else val for val in data['RESULT']]
        return data

    def get_timeline_treatment_rad(
        self, df_map: pd.DataFrame, df_file: pd.DataFrame
    ) -> dict:
        """Get TIMELINE-TREATMENT-RAD file data.

        Args:
            df_map (pd.DataFrame): variable to cBioPortal mapping info
            df_file (pd.DataFrame): data file to Synapse ID mapping

        Returns:
            dict: TIMELINE-TREATMENT-RAD data
        """

        timeline_infodf = df_map.query('sampleType == "TIMELINE-TREATMENT-RT"').merge(
            df_file, on="dataset", how="left"
        )
        timeline_infodf = pd.concat(
            [
                timeline_infodf,
                pd.DataFrame(
                    {
                        "code": "rt_rt_int",
                        "sampleType": "TIMELINE-TREATMENT-RT",
                        "dataset": "Cancer-Directed Radiation Therapy dataset",
                        "cbio": "TEMP",
                    },
                    index=["rt_rt_int"],
                ),
            ]
        )
        timeline_infodf.index = timeline_infodf["code"]

        treatment_rad_data = self.create_fixed_timeline_files(
            timeline_infodf, "TIMELINE-TREATMENT-RT"
        )
        rad_df = treatment_rad_data["df"]
        rad_df["STOP_DATE"] = rad_df["START_DATE"] + rad_df["TEMP"]
        rad_df = rad_df[rad_df["INDEX_CANCER"] == "Yes"]
        rad_df["EVENT_TYPE"] = "Treatment"
        rad_df["TREATMENT_TYPE"] = "Radiation Therapy"
        del rad_df["INDEX_CANCER"]
        del rad_df["TEMP"]

        return rad_df

    def get_timeline_dx(self, df_map: pd.DataFrame, df_file: pd.DataFrame) -> dict:
        """Get TIMELINE-DX file data.

        Args:
            df_map (pd.DataFrame): variable to cBioPortal mapping info
            df_file (pd.DataFrame): data file to Synapse ID mapping

        Returns:
            dict: TIMELINE-DX data
        """
        timeline_infodf = df_map.query('sampleType == "TIMELINE-DX"').merge(
            df_file, on="dataset", how="left"
        )
        timeline_infodf.index = timeline_infodf["code"]
        cancerdx_data = self.create_fixed_timeline_files(
            timeline_infodf, "TIMELINE-DX", filter_start=False
        )
        cancerdx_data["df"] = fill_cancer_dx_start_date(cancerdx_data["df"])
        cancerdx_data["df"] = cancerdx_data["df"][
            ~cancerdx_data["df"]["START_DATE"].isnull()
        ]
        return cancerdx_data

    def get_timeline_pathology(
        self, df_map: pd.DataFrame, df_file: pd.DataFrame
    ) -> dict:
        """Get TIMELINE-PATHOLOGY file data.

        Args:
            df_map (pd.DataFrame): variable to cBioPortal mapping info
            df_file (pd.DataFrame): data file to Synapse ID mapping

        Returns:
            dict: TIMELINE-PATHOLOGY data
        """
        timeline_infodf = df_map.query('sampleType == "TIMELINE-PATHOLOGY"').merge(
            df_file, on="dataset", how="left"
        )
        timeline_infodf.index = timeline_infodf["code"]
        pathology_data = self.create_fixed_timeline_files(
            timeline_infodf, "TIMELINE-PATHOLOGY"
        )
        return pathology_data

    def get_timeline_sample(self, df_map: pd.DataFrame, df_file: pd.DataFrame) -> dict:
        """Get TIMELINE-SAMPLE file data.

        Args:
            df_map (pd.DataFrame): variable to cBioPortal mapping info
            df_file (pd.DataFrame): data file to Synapse ID mapping

        Returns:
            dict: TIMELINE-SAMPLE data
        """
        timeline_infodf = df_map.query('sampleType == "TIMELINE-SAMPLE"').merge(
            df_file, on="dataset", how="left"
        )
        timeline_infodf.index = timeline_infodf["code"]
        acquisition_data = self.create_fixed_timeline_files(
            timeline_infodf, "TIMELINE-SAMPLE"
        )

        # acquisition_path = os.path.join(
        #     self._SPONSORED_PROJECT, "data_timeline_sample_acquisition.txt"
        # )
        # acquisitiondf = acquisition_data['df']
        # keep_idx = acquisitiondf['TEMP'] == acquisitiondf['PATH_PROC_NUMBER']
        # acquisitiondf.drop(columns="TEMP", inplace=True)
        # acquisition_data['df'] = acquisitiondf[keep_idx]

        # TODO: Can add getting of samples with NULL start dates in
        # self.create_fixed_timeline_files
        null_dates_idx = acquisition_data["df"]["START_DATE"].isnull()
        if null_dates_idx.any():
            logging.warning(
                "timeline sample with null START_DATE: {}".format(
                    ", ".join(acquisition_data["df"]["SAMPLE_ID"][null_dates_idx])
                )
            )
            acquisition_data["df"] = acquisition_data["df"][~null_dates_idx]
        return acquisition_data

    def get_timeline_medonc(self, df_map: pd.DataFrame, df_file: pd.DataFrame) -> dict:
        """Get TIMELINE-MEDONC file data.

        Args:
            df_map (pd.DataFrame): variable to cBioPortal mapping info
            df_file (pd.DataFrame): data file to Synapse ID mapping

        Returns:
            dict: TIMELINE-MEDONC data
        """
        timeline_infodf = df_map.query('sampleType == "TIMELINE-MEDONC"').merge(
            df_file, on="dataset", how="left"
        )
        timeline_infodf.index = timeline_infodf["code"]
        dict_medonc = self.create_fixed_timeline_files(
            timeline_infodf, "TIMELINE-MEDONC"
        )
        return dict_medonc

    def get_timeline_imaging(self, df_map: pd.DataFrame, df_file: pd.DataFrame) -> dict:
        """Get TIMELINE-IMAGING file data.

        Args:
            df_map (pd.DataFrame): variable to cBioPortal mapping info
            df_file (pd.DataFrame): data file to Synapse ID mapping

        Returns:
            dict: TIMELINE-IMAGING data
        """
        timeline_infodf = df_map.query('sampleType == "TIMELINE-IMAGING"').merge(
            df_file, on="dataset", how="left"
        )
        timeline_infodf.index = timeline_infodf["code"]
        dict_data = self.create_fixed_timeline_files(
            timeline_infodf, "TIMELINE-IMAGING"
        )
        return dict_data

    def get_timeline_sequence(
        self, df_map: pd.DataFrame, df_file: pd.DataFrame
    ) -> dict:
        """Get TIMELINE-SEQUENCE file data.

        Args:
            df_map (pd.DataFrame): variable to cBioPortal mapping info
            df_file (pd.DataFrame): data file to Synapse ID mapping

        Returns:
            dict: TIMELINE-SEQUENCE data
        """
        timeline_infodf = df_map.query('sampleType == "TIMELINE-SEQUENCE"').merge(
            df_file, on="dataset", how="left"
        )
        timeline_infodf.index = timeline_infodf["code"]
        # GEN-94: Add these cancer level dataset fields to get the
        # index cancer dob_ca_dx_days value to calculate START_DATE
        # and dob_cpt_report_days value
        df_timeline = pd.concat(
            [
                timeline_infodf,
                pd.DataFrame(
                    [
                        {
                            "code": "redcap_ca_index",
                            "sampleType": "TIMELINE-SEQUENCE",
                            "dataset": "Cancer-level dataset",
                            "cbio": "INDEX_CANCER",
                            "id": "syn22296816",  # HACK: hard coded synapse id
                        },
                        {
                            "code": "dob_ca_dx_days",
                            "sampleType": "TIMELINE-SEQUENCE",
                            "dataset": "Cancer-level dataset",
                            "cbio": "CA_DX_DAYS",
                            "id": "syn22296816",  # HACK: hard coded synapse id
                        },
                        {
                            "code": "dob_cpt_report_days",
                            "sampleType": "TIMELINE-SEQUENCE",
                            "dataset": "Cancer panel test level dataset",
                            "cbio": "DPT_REPORT_DAYS",
                            "id": "syn22296816",  # HACK: hard coded synapse id
                        },
                    ],
                    index=["redcap_ca_index", "dob_ca_dx_days", "dob_cpt_report_days"],
                ),
            ]
        )
        sequence_data = self.create_fixed_timeline_files(
            df_timeline, "TIMELINE-SEQUENCE"
        )
        # HACK: Manually calculate the START_DATE based on criteria defined
        # in GEN-94
        # sequence_data["df"]["START_DATE"]
        seq_df = sequence_data["df"]
        index_seq_df = seq_df[seq_df["INDEX_CANCER"] == "Yes"]
        index_seq_df["START_DATE"] = (
            index_seq_df["DPT_REPORT_DAYS"] - index_seq_df["CA_DX_DAYS"]
        )
        del seq_df["START_DATE"]
        seq_df = seq_df.merge(index_seq_df[["SAMPLE_ID", "START_DATE"]], on="SAMPLE_ID")
        # Delete unneeded columns
        del seq_df["DPT_REPORT_DAYS"]
        del seq_df["INDEX_CANCER"]
        del seq_df["CA_DX_DAYS"]
        # reorder Columns to match requirement
        cols_to_order = ["PATIENT_ID", "START_DATE", "STOP_DATE", "EVENT_TYPE"]
        cols_to_order.extend(seq_df.columns.drop(cols_to_order).tolist())
        seq_df = seq_df[cols_to_order]
        seq_df.drop_duplicates(inplace=True)
        sequence_data["df"] = seq_df
        # sequence_path = os.path.join(
        #     self._SPONSORED_PROJECT, "data_timeline_sequencing.txt"
        # )
        # self.write_and_storedf(
        #     seq_df, sequence_path, used_entities=sequence_data["used"]
        # )

        return sequence_data

    def get_timeline_lab(
        self, df_map: pd.DataFrame, df_file: pd.DataFrame
    ) -> pd.DataFrame:
        """Get TIMELINE-LAB file data.

        Args:
            df_map (pd.DataFrame): variable to cBioPortal mapping info
            df_file (pd.DataFrame): data file to Synapse ID mapping

        Returns:
            pd.DataFrame: TIMELINE-LAB data
        """
        timeline_infodf = df_map.query('sampleType == "TIMELINE-LAB"').merge(
            df_file, on="dataset", how="left"
        )
        timeline_infodf.index = timeline_infodf["code"]
        dict_data = self.create_fixed_timeline_files(timeline_infodf, "TIMELINE-LAB")
        return dict_data

    def get_survival(self, df_map: pd.DataFrame, df_file: pd.DataFrame) -> dict:
        """Get SURVIVAL file data.

        Args:
            df_map (pd.DataFrame): variable to cBioPortal mapping info
            df_file (pd.DataFrame): data file to Synapse ID mapping

        Returns:
            dict: SURVIVAL data
        """
        df_info = df_map.query('sampleType == "SURVIVAL"').merge(
            df_file, on="dataset", how="left"
        )
        df_info.index = df_info["code"]
        df_info = pd.concat(
            [
                df_info,
                pd.DataFrame(
                    {
                        "code": "first_index_ca_days",
                        "sampleType": "SURVIVAL",
                        "dataset": "Cancer-level index dataset",
                        "cbio": "CANCER_INDEX",
                    },
                    index=["first_index_ca_days"],
                ),
            ]
        )
        df_info.index = df_info["code"]

        df_info_survial = df_info[df_info["sampleType"] == "SURVIVAL"]
        dict_data_survial = get_file_data(
            self.syn, df_info_survial, "SURVIVAL", cohort=self._SPONSORED_PROJECT
        )
        df_raw_survival = dict_data_survial["df"]
        df_final_survival = self.configure_clinicaldf(df_raw_survival, df_info_survial)
        # Only take rows where cancer index is null
        df_final_survival = df_final_survival[
            df_final_survival["CANCER_INDEX"].isnull()
        ]
        # Remove cancer index column
        del df_final_survival["CANCER_INDEX"]
        # remove a row if patient ID is duplicated and PFS_I_ADV_STATUS is null or empty
        # tested on current survival data file and produces unique patient list
        if "PFS_I_ADV_STATUS" in df_final_survival.columns:
            pfs_not_null_idx = ~df_final_survival["PFS_I_ADV_STATUS"].isnull()
            pfs_not_blank_idx = df_final_survival["PFS_I_ADV_STATUS"] != ""
            nondup_patients_idx = ~df_final_survival["PATIENT_ID"].duplicated(
                keep=False
            )
            df_final_survival = df_final_survival[
                (pfs_not_null_idx & pfs_not_blank_idx) | (nondup_patients_idx)
            ]

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
            syn=self.syn,
            cohort=self._SPONSORED_PROJECT,
            synid_file_grs=self._GRS_SYNID,
            synid_table_prissmm=self._PRISSMM_SYNID,
        )
        regimens_data = create_regimens(
            self.syn,
            regimen_infodf,
            mapping=drug_mapping,
            top_x_regimens=20,
            cohort=self._SPONSORED_PROJECT,
        )

        survival_info = pd.concat([infodf, regimens_data["info"]])
        # Only patients and samples that exist in the
        # sponsored project uploads are going to be pulled into the SP project
        subset_survivaldf = df_final_survival[
            df_final_survival["PATIENT_ID"].isin(self.genie_clinicaldf["PATIENT_ID"])
        ]
        del subset_survivaldf["SP"]
        subset_survivaldf = remap_os_values(df=subset_survivaldf)
        subset_survivaldf = remap_pfs_values(df=subset_survivaldf)
        cols_to_order = ["PATIENT_ID"]
        cols_to_order.extend(subset_survivaldf.columns.drop(cols_to_order).tolist())

        return {"df": subset_survivaldf[cols_to_order], "survival_info": survival_info}

    def get_survival_treatment(
        self, df_map: pd.DataFrame, df_file: pd.DataFrame
    ) -> pd.DataFrame:
        """Get SURVIVAL and REGIMEN file data.

        Args:
            df_map (pd.DataFrame): variable to cBioPortal mapping info
            df_file (pd.DataFrame): data file to Synapse ID mapping

        Returns:
            pd.DataFrame: SURVIVAL and REGIMEN data
        """
        # Regimen mapping values
        regimen_infodf = df_map.query('sampleType == "REGIMEN"').merge(
            df_file, on="dataset", how="left"
        )
        regimen_infodf.index = regimen_infodf["code"]

        drug_mapping = get_drug_mapping(
            syn=self.syn,
            cohort=self._SPONSORED_PROJECT,
            synid_file_grs=self._GRS_SYNID,
            synid_table_prissmm=self._PRISSMM_SYNID,
        )
        regimens_data = create_regimens(
            self.syn,
            regimen_infodf,
            mapping=drug_mapping,
            top_x_regimens=20,
            cohort=self._SPONSORED_PROJECT,
        )

        df_survival_treatment = regimens_data["df"]
        df_survival_treatment = remap_os_values(df=df_survival_treatment)
        df_survival_treatment = remap_pfs_values(df=df_survival_treatment)
        cols_to_order = ["PATIENT_ID"]
        cols_to_order.extend(df_survival_treatment.columns.drop(cols_to_order).tolist())
        # Retract patients from survival treatment file
        df_survival_treatment = df_survival_treatment[
            df_survival_treatment["PATIENT_ID"].isin(
                self.genie_clinicaldf["PATIENT_ID"]
            )
        ]
        return df_survival_treatment[cols_to_order]

    def get_patient(self, df_map: pd.DataFrame, df_file: pd.DataFrame) -> pd.DataFrame:
        """Patient data file.

        Custom rule: Fix patient duplicated values due to cancer index DOB

        Args:
            df_map (pd.DataFrame): variable to cBioPortal mapping info
            df_file (pd.DataFrame): data file to Synapse ID mapping

        Returns:
            pd.DataFrame: PATIENT data
        """
        df_info_patient = df_map.query('sampleType == "PATIENT"').merge(
            df_file, on="dataset", how="left"
        )
        df_info_patient.index = df_info_patient["code"]
        df_info_patient = pd.concat(
            [
                df_info_patient,
                pd.DataFrame(
                    {
                        "code": "redcap_ca_index",
                        "sampleType": "PATIENT",
                        "dataset": "Cancer-level dataset",
                    },
                    index=["redcap_ca_index"],
                ),
            ],
            ignore_index=True,
        )
        df_info_patient.index = df_info_patient["code"]

        dict_patient = get_file_data(
            self.syn, df_info_patient, "PATIENT", cohort=self._SPONSORED_PROJECT
        )

        df_patient = dict_patient["df"]
        df_patient = df_patient[df_patient["redcap_ca_index"] == "Yes"]
        df_patient.drop(columns="redcap_ca_index", inplace=True)
        df_patient_final = self.configure_clinicaldf(df_patient, df_info_patient)

        df_patient_subset = df_patient_final[
            df_patient_final["PATIENT_ID"].isin(self.genie_clinicaldf["PATIENT_ID"])
        ]

        # Fix patient duplicated values due to cancer index DOB
        # Take the larger DX_LASTALIVE_INT_MOS value for all records
        # subset_patientdf.sort_values("DX_LASTALIVE_INT_MOS", inplace=True,
        #                              ascending=False)
        df_patient_subset.drop_duplicates("PATIENT_ID", inplace=True)
        duplicated = df_patient_subset.PATIENT_ID.duplicated()
        if duplicated.any():
            logging.warning(
                "DUPLICATED PATIENT_IDs: {}".format(
                    ",".join(df_patient_subset["PATIENT_ID"][duplicated])
                )
            )

        del df_patient_subset["SP"]
        cols_to_order = ["PATIENT_ID"]
        cols_to_order.extend(df_patient_subset.columns.drop(cols_to_order).tolist())

        df_patient_subset = hack_remap_laterality(df_patient_subset=df_patient_subset)

        return df_patient_subset[cols_to_order]

    def get_sample(self, df_map: pd.DataFrame, df_file: pd.DataFrame) -> pd.DataFrame:
        """SAMPLE data file

        Args:
            df_map (pd.DataFrame): variable to cBioPortal mapping info
            df_file (pd.DataFrame): data file to Synapse ID mapping

        Returns:
            pd.DataFrame: SAMPLE data
        """
        df_info_sample = df_map.query('sampleType == "SAMPLE"').merge(
            df_file, on="dataset", how="left"
        )
        df_info_sample.index = df_info_sample["code"]
        dict_sample = get_file_data(
            self.syn, df_info_sample, "SAMPLE", cohort=self._SPONSORED_PROJECT
        )

        df_sample = dict_sample["df"]
        del df_sample["path_proc_number"]
        df_sample_final = self.configure_clinicaldf(df_sample, df_info_sample)

        df_sample_subset = df_sample_final[
            df_sample_final["SAMPLE_ID"].isin(self.genie_clinicaldf["SAMPLE_ID"])
        ]
        del df_sample_subset["SP"]
        days_to_years_col = [
            "AGE_AT_SEQ_REPORT_YEARS",
            "CPT_ORDER_INT",
            "CPT_REPORT_INT",
        ]
        for col in days_to_years_col:
            # not all columns could exist, so check if column exists
            if col in df_sample_subset:
                years = df_sample_subset[col].apply(change_days_to_years)
                df_sample_subset[col] = years
        df_sample_subset["AGE_AT_SEQUENCING"] = df_sample_subset[
            "AGE_AT_SEQUENCING"
        ].apply(math.floor)
        # Remove SAMPLE_TYPE and CPT_SEQ_DATE because the values are incorrect
        del df_sample_subset["CPT_SEQ_DATE"]
        # Obtain this information from the main GENIE cohort
        df_sample_subset = df_sample_subset.merge(
            self.genie_clinicaldf[["SAMPLE_ID", "SEQ_YEAR"]],
            on="SAMPLE_ID",
            how="left",
        )
        df_sample_subset.rename(columns={"SEQ_YEAR": "CPT_SEQ_DATE"}, inplace=True)
        df_sample_subset.sort_values("PDL1_POSITIVE_ANY", ascending=False, inplace=True)
        df_sample_subset.drop_duplicates("SAMPLE_ID", inplace=True)

        return df_sample_subset

    def create_and_write_case_lists(
        self, subset_sampledf: pd.DataFrame, subset_patientdf: pd.DataFrame, used: list
    ) -> None:
        """Create, write, and, if applicable, store case list files.

        Args:
            subset_sampledf: Dataframe of sample information
            subset_patientdf: Dataframe of patient information
            used: Synapse IDs used to construct the case file data.
        """
        # Merged clinical file is created here because it is needed
        # for the case lists
        # Remove oncotree code here, because no longer need it
        merged_clinicaldf = subset_sampledf.merge(
            subset_patientdf, on="PATIENT_ID", how="outer"
        )
        missing_sample_idx = merged_clinicaldf["SAMPLE_ID"].isnull()
        # Make sure there are no missing sample ids
        if sum(missing_sample_idx) > 0:
            logging.info(
                "MISSING SAMPLE_ID for: {}".format(
                    ",".join(merged_clinicaldf["PATIENT_ID"][missing_sample_idx])
                )
            )
            merged_clinicaldf = merged_clinicaldf[~missing_sample_idx]

        # upload samples that are not part of the main GENIE cohort
        if merged_clinicaldf.get("SAMPLE_ID") is not None:
            logging.info("Samples not in GENIE clinical databases (SP and normal)")
            not_found_samples = merged_clinicaldf["SAMPLE_ID"][
                ~merged_clinicaldf["SAMPLE_ID"].isin(self.genie_clinicaldf["SAMPLE_ID"])
            ]
            if not not_found_samples.empty:
                logging.info(not_found_samples[~not_found_samples.isnull()])
                not_found_samples.to_csv("notfoundsamples.csv")
                if self.upload:
                    self.syn.store(
                        File(
                            "notfoundsamples.csv", parent=self._SP_REDCAP_EXPORTS_SYNID
                        )
                    )
        # Hard coded most up to date oncotree version
        oncotreelink = self.syn.get("syn13890902").externalURL
        # Use the old oncotree link for now
        # TODO: need to update oncotree link for 11.0 public
        oncotreelink = (
            "http://oncotree.mskcc.org/api/tumorTypes/tree?version=oncotree_2018_06_01"
        )
        oncotree_dict = process_functions.get_oncotree_code_mappings(oncotreelink)
        # Map cancer type and cancer type detailed
        # This is to create case list files
        merged_clinicaldf["CANCER_TYPE"] = [
            oncotree_dict[code.upper()].get("CANCER_TYPE", float("nan"))
            for code in merged_clinicaldf["ONCOTREE_CODE"]
        ]
        merged_clinicaldf["CANCER_TYPE_DETAILED"] = [
            oncotree_dict[code.upper()].get("CANCER_TYPE_DETAILED", float("nan"))
            for code in merged_clinicaldf["ONCOTREE_CODE"]
        ]
        merged_clinicaldf["ONCOTREE_PRIMARY_NODE"] = [
            oncotree_dict[code.upper()].get("ONCOTREE_PRIMARY_NODE", float("nan"))
            for code in merged_clinicaldf["ONCOTREE_CODE"]
        ]
        merged_clinicaldf["ONCOTREE_SECONDARY_NODE"] = [
            oncotree_dict[code.upper()].get("ONCOTREE_SECONDARY_NODE", float("nan"))
            for code in merged_clinicaldf["ONCOTREE_CODE"]
        ]
        # Remove duplicated sample ids (there shouldn't be any)
        merged_clinicaldf = merged_clinicaldf.drop_duplicates("SAMPLE_ID")
        merged_clinicaldf.to_csv(
            os.path.join(self._SPONSORED_PROJECT, "data_clinical.txt"),
            index=False,
            sep="\t",
        )

        # Create case lists
        case_list_path = os.path.join(self._SPONSORED_PROJECT, "case_lists")

        if not os.path.exists(case_list_path):
            os.mkdir(case_list_path)
        else:
            caselists = os.listdir(case_list_path)
            for caselist in caselists:
                os.remove(os.path.join(case_list_path, caselist))

        # Write out cases sequenced so people can tell
        # which samples were sequenced
        assay_info = self.syn.tableQuery(
            f"select * from {self._ASSAY_SYNID}",
            includeRowIdAndRowVersion=False,
            separator="\t",
        )
        create_case_lists.main(
            os.path.join(self._SPONSORED_PROJECT, "data_clinical.txt"),
            assay_info.filepath,
            case_list_path,
            f"{self._SPONSORED_PROJECT.lower()}_genie_bpc",
        )

        case_list_files = os.listdir(case_list_path)
        for casepath in case_list_files:
            casepath = os.path.join(case_list_path, casepath)
            if self.upload:
                file_ent = File(casepath, parent=self.cbioportal_folders["case_lists"])
                self.syn.store(
                    file_ent,
                    used=used,
                    executed=self._GITHUB_REPO,
                )

    def run(self):
        """Runs the redcap export to export all files"""

        logging.info("creating release folders...")
        create_release_folders(cohort=self._SPONSORED_PROJECT)

        logging.info("reading variable to cBio mapping...")
        redcap_to_cbiomappingdf = get_bpc_to_cbio_mapping_df(
            self.syn,
            cohort=self._SPONSORED_PROJECT,
            synid_table_cbio=self._REDCAP_TO_CBIOMAPPING_SYNID,
        )

        logging.info("reading and dataset label mappings...")
        data_tablesdf = get_data_file_synapse_id_df(
            syn=self.syn, synid_table_files=self._DATA_TABLE_IDS
        )

        logging.info("writing TIMELINE-TREATMENT...")
        treatment_data = self.get_timeline_treatment(
            df_map=redcap_to_cbiomappingdf, df_file=data_tablesdf
        )
        if self._SPONSORED_PROJECT not in ["BrCa", "CRC", "NSCLC"]:
            logging.info("writing TIMELINE-TREATMENT-RT...")
            rad_df = self.get_timeline_treatment_rad(
                df_map=redcap_to_cbiomappingdf, df_file=data_tablesdf
            )
            treatment_data["df"] = pd.concat([treatment_data["df"], rad_df])
        else:
            logging.info("skipping TIMELINE-TREATMENT-RT")
        self.write_and_storedf(
            treatment_data["df"],
            os.path.join(self._SPONSORED_PROJECT, "data_timeline_treatment.txt"),
            used_entities=treatment_data["used"],
        )

        logging.info("writing TIMELINE-DX...")
        cancerdx_data = self.get_timeline_dx(
            df_map=redcap_to_cbiomappingdf, df_file=data_tablesdf
        )
        self.write_and_storedf(
            cancerdx_data["df"],
            os.path.join(self._SPONSORED_PROJECT, "data_timeline_cancer_diagnosis.txt"),
            used_entities=cancerdx_data["used"],
        )

        logging.info("writing TIMELINE-PATHOLOGY...")
        pathology_data = self.get_timeline_pathology(
            df_map=redcap_to_cbiomappingdf, df_file=data_tablesdf
        )
        self.write_and_storedf(
            pathology_data["df"],
            os.path.join(self._SPONSORED_PROJECT, "data_timeline_pathology.txt"),
            used_entities=pathology_data["used"],
        )

        logging.info("writing TIMELINE-SAMPLE-ACQUISITION...")
        acquisition_data = self.get_timeline_sample(
            df_map=redcap_to_cbiomappingdf, df_file=data_tablesdf
        )
        self.write_and_storedf(
            acquisition_data["df"],
            os.path.join(
                self._SPONSORED_PROJECT, "data_timeline_sample_acquisition.txt"
            ),
            used_entities=acquisition_data["used"],
        )

        logging.info("writing TIMELINE-MEDONC...")
        medonc_data = self.get_timeline_medonc(
            df_map=redcap_to_cbiomappingdf, df_file=data_tablesdf
        )
        self.write_and_storedf(
            medonc_data["df"],
            os.path.join(self._SPONSORED_PROJECT, "data_timeline_medonc.txt"),
            used_entities=medonc_data["used"],
        )

        logging.info("writing TIMELINE-IMAGING...")
        imaging_data = self.get_timeline_imaging(
            df_map=redcap_to_cbiomappingdf, df_file=data_tablesdf
        )
        self.write_and_storedf(
            imaging_data["df"],
            os.path.join(self._SPONSORED_PROJECT, "data_timeline_imaging.txt"),
            used_entities=imaging_data["used"],
        )

        logging.info("writing TIMELINE-SEQUENCE...")
        sequence_data = self.get_timeline_sequence(
            df_map=redcap_to_cbiomappingdf, df_file=data_tablesdf
        )
        self.write_and_storedf(
            df=sequence_data["df"],
            filepath=os.path.join(
                self._SPONSORED_PROJECT, "data_timeline_sequencing.txt"
            ),
            used_entities=sequence_data["used"],
        )

        if self._SPONSORED_PROJECT not in ["NSCLC", "BLADDER"]:
            logging.info("writing TIMELINE-LABTEST...")
            lab_data = self.get_timeline_lab(
                df_map=redcap_to_cbiomappingdf, df_file=data_tablesdf
            )
            self.write_and_storedf(
                df=lab_data["df"],
                filepath=os.path.join(
                    self._SPONSORED_PROJECT, "data_timeline_labtest.txt"
                ),
                used_entities=lab_data["used"],
            )
        else:
            logging.info("skipping TIMELINE-LABTEST...")

        logging.info("writing TIMELINE-PERFORMANCE...")
        if self._SPONSORED_PROJECT in ["BLADDER"]:
            performance_data = self.get_timeline_performance(
                df_map=redcap_to_cbiomappingdf, df_file=data_tablesdf
            )
            self.write_and_storedf(
                df=performance_data["df"],
                filepath=os.path.join(
                    self._SPONSORED_PROJECT, "data_timeline_performance_status.txt"
                ),
                used_entities=performance_data["used"],
            )
        else:
            logging.info("skipping TIMELINE-PERFORMANCE...")

        logging.info("writing CLINICAL-SURVIVAL...")
        final_survival_data = self.get_survival(
            df_map=redcap_to_cbiomappingdf, df_file=data_tablesdf
        )
        survival_path = self.write_clinical_file(
            final_survival_data["df"],
            final_survival_data["survival_info"],
            "supp_survival",
        )
        survival_used = get_synid_data(
            df_map=redcap_to_cbiomappingdf,
            df_file=data_tablesdf,
            sampletype=["SURVIVAL", "REGIMEN"],
            cohort=self._SPONSORED_PROJECT,
        )
        if self.upload:
            survival_fileent = File(
                survival_path, parent=self.cbioportal_folders["release"]
            )
            survival_ent = self.syn.store(
                survival_fileent, used=survival_used, executed=self._GITHUB_REPO
            )

        logging.info("writing CLINICAL-SURVIVAL-TREATMENT...")
        df_survival_treatment = self.get_survival_treatment(
            df_map=redcap_to_cbiomappingdf, df_file=data_tablesdf
        )
        surv_treatment_path = self.write_clinical_file(
            df_survival_treatment,
            final_survival_data["survival_info"],
            "supp_survival_treatment",
        )
        if self.upload:
            survival_treatment_fileent = File(
                surv_treatment_path, parent=self.cbioportal_folders["release"]
            )
            survival_treatment_fileent = self.syn.store(
                survival_treatment_fileent,
                used=survival_used,
                executed=self._GITHUB_REPO,
            )

        logging.info("writing CLINICAL-SAMPLE...")
        df_sample_final = self.get_sample(
            df_map=redcap_to_cbiomappingdf, df_file=data_tablesdf
        )
        sample_path = self.write_clinical_file(
            df_sample_final, final_survival_data["survival_info"], "sample"
        )
        if self.upload:
            sample_fileent = File(
                sample_path, parent=self.cbioportal_folders["release"]
            )
            sample_used = get_synid_data(
                df_map=redcap_to_cbiomappingdf,
                df_file=data_tablesdf,
                sampletype=["SAMPLE"],
                cohort=self._SPONSORED_PROJECT,
            )
            sample_ent = self.syn.store(
                sample_fileent, used=sample_used, executed=self._GITHUB_REPO
            )

        logging.info("writing CLINICAL-PATIENT...")
        df_patient_final = self.get_patient(
            df_map=redcap_to_cbiomappingdf, df_file=data_tablesdf
        )
        patient_path = self.write_clinical_file(
            df_patient_final, final_survival_data["survival_info"], "patient"
        )
        if self.upload:
            logging.info("uploading clinical data files to Synapse...")
            patient_fileent = File(
                patient_path, parent=self.cbioportal_folders["release"]
            )
            patient_used = get_synid_data(
                df_map=redcap_to_cbiomappingdf,
                df_file=data_tablesdf,
                sampletype=["PATIENT"],
                cohort=self._SPONSORED_PROJECT,
            )
            patient_ent = self.syn.store(
                patient_fileent, used=patient_used, executed=self._GITHUB_REPO
            )

        logging.info("writing genomic data files...")
        self.create_and_write_maf(df_sample_final["SAMPLE_ID"])
        dict_cna = self.create_and_write_cna(df_sample_final["SAMPLE_ID"])
        self.create_and_write_genematrix(df_sample_final, dict_cna["cna_samples"])
        self.create_and_write_fusion(df_sample_final["SAMPLE_ID"])
        self.create_and_write_seg(df_sample_final["SAMPLE_ID"])
        self.create_and_write_sv(df_sample_final["SAMPLE_ID"])

        ids = get_synid_data(
            df_map=redcap_to_cbiomappingdf,
            df_file=data_tablesdf,
            sampletype=["PATIENT", "SAMPLE"],
            cohort=self._SPONSORED_PROJECT,
        )
        self.create_and_write_case_lists(
            subset_sampledf=df_sample_final, subset_patientdf=df_patient_final, used=ids
        )
        self.create_and_write_gene_panels(df_sample_final["SEQ_ASSAY_ID"].unique())

        logging.info("writing metadata files...")
        metadata_files = self.create_bpc_cbio_metafiles()
        if self.upload:
            for metadata_file in metadata_files:
                file_ent = File(
                    metadata_file, parent=self.cbioportal_folders["release"]
                )
                self.syn.store(
                    file_ent,
                    executed=self._GITHUB_REPO,
                )

        logging.info("cBioPortal validation")
        cmd = [
            "python",
            os.path.join(
                self.cbiopath, "core/src/main/scripts/importer/validateData.py"
            ),
            "-s",
            self._SPONSORED_PROJECT,
            "-n",
        ]
        subprocess.run(cmd)
