from abc import ABCMeta
from dataclasses import dataclass
from functools import cached_property
from typing import List

import pandas as pd
from synapseclient import Synapse

from geniesp.config import BpcConfig


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


@dataclass
class Extract:
    """Timeline data class."""
    # Config = TimelineConfig
    bpc_config: BpcConfig
    syn: Synapse
    sample_type: str = None

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

    @cached_property
    def mapping_df(self) -> pd.DataFrame:
        """Extract relevant portions of the mapping table for the sponsored project
        variables and cBioPortal variable mappings and return as a data frame.

        Args:
            syn (Synapse): Synapse connection
            cohort (str): sponsored project label
            synid_table_cbio (str): Synapse ID of table containing variable to cbio mapping

        Returns:
            pd.DataFrame: data frame of all mapping columns for released variables
        """
        redcap_to_cbiomapping = self.syn.tableQuery(
            f"SELECT * FROM {self.bpc_config.redcap_to_cbio_mapping_synid} where "
            f"{self.bpc_config.cohort} is true AND sampleType <> 'TIMELINE-STATUS'"
        )
        redcap_to_cbiomappingdf = redcap_to_cbiomapping.asDataFrame()
        return redcap_to_cbiomappingdf

    @cached_property
    def data_tables_df(self) -> pd.DataFrame:
        """Get mapping of dataset files to Synapse IDs.

        Args:
            syn (Synapse): Synapse connection
            synid_table_files (str): Synapse ID of table containing data file to Synapse ID mapping

        Returns:
            pd.DataFrame: data frame with two columns representing the Synapse ID and dataset label
        """
        data_tables = self.syn.tableQuery(f"SELECT id, dataset FROM {self.bpc_config.data_tables_id}")
        data_tablesdf = data_tables.asDataFrame()
        return data_tablesdf

    @cached_property
    def survival_info_df(self):
        # Patient and Sample mapping values
        df_map = self.mapping_df
        df_file = self.data_tables_df
        cohort = self.bpc_config.cohort
        prissm_synid = self.bpc_config.prissmm_synid

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
            cohort=cohort,
            synid_table_prissmm=prissm_synid
        )
        regimens_data = create_regimens(
            self.syn,
            regimen_infodf,
            mapping=drug_mapping,
            top_x_regimens=20,
            cohort=cohort
        )
        survival_info = pd.concat([infodf, regimens_data["info"]])
        return survival_info

    @cached_property
    def timeline_infodf(self):
        timeline_infodf = self.mapping_df.query(f'sampleType == "{self.sample_type}"').merge(
            self.data_tables_df, on="dataset", how="left"
        )
        timeline_infodf = pd.concat(
            [
                timeline_infodf,
                pd.DataFrame(
                    [
                        {
                            "code": "rt_rt_int",
                            "sampleType": "TIMELINE-TREATMENT-RT",
                            "dataset": "Cancer-Directed Radiation Therapy dataset",
                            "cbio": "TEMP",
                        },
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
                        {
                            "code": "first_index_ca_days",
                            "sampleType": "SURVIVAL",
                            "dataset": "Cancer-level index dataset",
                            "cbio": "CANCER_INDEX",
                        },
                        {
                            "code": "redcap_ca_index",
                            "sampleType": "PATIENT",
                            "dataset": "Cancer-level dataset",
                            "cbio": "INDEX_CANCER",
                        }
                    ],
                    index=["rt_rt_int", "redcap_ca_index", "dob_ca_dx_days", "dob_cpt_report_days", "first_index_ca_days", "redcap_ca_index"],
                ),
            ]
        )
        subset_infodf = timeline_infodf[timeline_infodf["sampleType"] == self.sample_type]
        subset_infodf.index = subset_infodf["code"]

        return subset_infodf

    def get_derived_variable_files(
        self
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
        used_entities = []
        derived_variable_entities = {}
        # Group by dataset because different datasets could have the
        # same variable
        datasets = self.timeline_infodf.groupby("dataset")
        for dataset, df in datasets:
            # Get synapse id
            synid = df["id"].unique()[0]
            file_entity = self.syn.get(synid)
            derived_variable_entities[dataset] = file_entity
            used_entities.append(f"{synid}.{file_entity.versionNumber}")

        return {"derived_variable_entities": derived_variable_entities, "used": used_entities}


    @cached_property
    def to_keep_samples_and_patients(self) -> pd.DataFrame:
        """Get main GENIE clinical samples and perform retraction
        against BPC sample and patient database. It's important that I use the
        sample database, because if a sample from a patient with multiple samples
        is retracted from the main GENIE clinical samples, the patient
        will still exist. (Jira: GEN-260)

        Returns:
            pd.DataFrame: main GENIE information for sponsored project samples
        """
        # genie_clinicaldb = self.syn.tableQuery(
        #     f"select SAMPLE_ID, PATIENT_ID, ONCOTREE_CODE, SEQ_ASSAY_ID, "
        #     f"SAMPLE_TYPE, SEQ_YEAR from {self._CLINICAL_SYNID}"
        # )
        # genie_clinicaldf = genie_clinicaldb.asDataFrame()
        # DECISION 06/2023: use the release file to do retractions
        # This is due to the most recent releases potentially having
        # samples retracted. The consortium release matched with the
        # public release (14.7-consortium <-> 14.0-public) must be used
        # due to SEQ_YEAR being used through the code.
        sample_synid = self.get_mg_synid(
            self.bpc_config.mg_release_synid, "data_clinical_sample.txt"
        )
        genie_clinicaldf = pd.read_csv(
            self.syn.get(sample_synid, followLink=True).path, sep="\t", comment="#"
        )
        # Filter out cfDNA samples
        genie_clinicaldf = genie_clinicaldf[
            genie_clinicaldf['SAMPLE_CLASS'] != "cfDNA"
        ]
        # BPC retraction database
        bpc_sample_retraction_db = self.syn.tableQuery(
            f"select SAMPLE_ID from {self.bpc_config.sample_retraction_synid} where "
            f"{self.bpc_config.cohort} is true"
        )
        bpc_sample_retractiondf = bpc_sample_retraction_db.asDataFrame()

        bpc_patient_retraction_db = self.syn.tableQuery(
            f"select record_id from {self.bpc_config.patient_retraction_synid} where "
            f"{self.bpc_config.cohort} is true"
        )
        bpc_patient_retraction_df = bpc_patient_retraction_db.asDataFrame()

        bpc_temp_patient_retraction_db = self.syn.tableQuery(
            f"select record_id from {self.bpc_config.temporary_patient_retraction_synid} where "
            f"cohort = '{self.bpc_config.cohort}'"
        )
        bpc_temp_patient_retraction_df = bpc_temp_patient_retraction_db.asDataFrame()

        retraction_at_release = self.syn.tableQuery(
            f"select patient_id from {self.bpc_config.retraction_at_release_synid} where "
            f"cohort = '{self.bpc_config.cohort}'"
        )
        retraction_at_release_df = retraction_at_release.asDataFrame()
        # Retract samples from sample retraction db
        keep_clinicaldf = genie_clinicaldf[
            ~genie_clinicaldf["SAMPLE_ID"].isin(bpc_sample_retractiondf["SAMPLE_ID"])
        ]
        # Retract patients from patient retraction db
        keep_clinicaldf = keep_clinicaldf[
            ~keep_clinicaldf["PATIENT_ID"].isin(bpc_patient_retraction_df["record_id"])
        ]
        # Retract patients for at release retraction table
        keep_clinicaldf = keep_clinicaldf[
            ~keep_clinicaldf["PATIENT_ID"].isin(retraction_at_release_df["patient_id"])
        ]
        # Retract patients from temporary patient retraction db
        keep_clinicaldf = keep_clinicaldf[
            ~keep_clinicaldf["PATIENT_ID"].isin(
                bpc_temp_patient_retraction_df["record_id"]
            )
        ]
        return keep_clinicaldf
