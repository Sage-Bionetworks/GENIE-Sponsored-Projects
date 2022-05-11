"""
Description: Validate the BPC to cBioPortal mapping file. 
Author: Haley Hunter-Zinck
Date: 2022-01-27
"""

import argparse
import logging
import re
from typing import Dict, List

import pandas as pd

import synapseclient
from synapseclient import Synapse
from synapseclient.core.exceptions import (
    SynapseAuthenticationError,
    SynapseNoCredentialsError,
)
import yaml


def get_sor_column_name(syn: Synapse, synid_table_rel: str, cohort: str, release: str) -> str:
    """Get Scope of Release column name for requested cohort and release pair.

    Args:
        syn (Synapse): Synapse object
        synid_table_rel (str): Synapse ID of BPC release info table
        cohort (str): cohort label
        release (str): release label

    Returns:
        str: SOR column for cohort and release
    """

    release_version = release.split("-")[0]
    release_type = release.split("-")[1]
    query = f"SELECT sor_column FROM {synid_table_rel} WHERE cohort = '{cohort}' AND release_version = '{release_version}' AND release_type = '{release_type}'"
    df = syn.tableQuery(query).asDataFrame()
    return df["sor_column"][0]


def get_codes_to_remove(codes: List) -> List:
    """Remove codes with wildcards or nan.

    Args:
        codes (list): list of codes

    Returns:
        list: codes that are nan or that contain wildcards
    """
    code_remove = []
    for code in codes:
        if bool(re.match(r"^.+[*]$", str(code).strip())):
            code_remove.append(code)
        elif pd.isna(code):
            code_remove.append(code)
    return code_remove


def check_code_name_empty(
    df: pd.DataFrame, syn: Synapse, config: Dict, cohort: str, release: str
) -> List:
    """Check for any code that is empty.

     Args:
        df (pd.DataFrame): dataframe representing map
        syn (Synapse): Synapse object
        config (dict): configuration parameters
        cohort (str): cohort label to check
        release (str): release label to check

    Returns:
        list: missing code names
    """
    empty = df.loc[pd.isna(df["code"])]["code"]
    return list(empty)


def check_code_name_absent(
    df: pd.DataFrame, syn: Synapse, config: Dict, cohort: str, release: str
) -> List:
    """Check for any code that is not code name that
    does not appear in its associated data file.

    Args:
        df (pd.DataFrame): dataframe representing map
        syn (Synapse): Synapse object
        config (dict): configuration parameters
        cohort (str): cohort label to check
        release (str): release label to check

    Returns:
        list: code names not found in datasets
    """
    absent = []

    query = f'SELECT id, dataset FROM {config["synapse"]["dataset"]["id"]} WHERE dataset IS NOT NULL'
    res = syn.tableQuery(query)

    # only examine released codes
    df = df[df[cohort].str.lower() == "y"]

    for row in res:

        synapse_id = row[0]
        dataset = row[1]

        data = pd.read_csv(syn.get(synapse_id)["path"], low_memory=False)
        code_data = data.columns

        # get codes associated with the dataset and of types derived or curated
        code_map = list(
            df.loc[
                (
                    (df["dataset"] == dataset)
                    & (
                        (df["data_type"].str.lower() == "derived")
                        | (df["data_type"].str.lower() == "curated")
                    )
                )
            ]["code"]
        )

        # do not check wildcard code names or NA code names
        code_remove = get_codes_to_remove(code_map)
        for code in code_remove:
            code_map.remove(code)

        absent.extend(list(set(code_map) - set(code_data)))
    return absent


def check_dataset_names(
    df: pd.DataFrame, syn: Synapse, config: Dict, cohort: str, release: str
) -> List:
    """Check for any dataset name that is not associated with a dataset.

    Args:
        df (pd.DataFrame): dataframe representing map
        syn (Synapse): Synapse object
        config (dict): configuration parameters
        cohort (str): cohort label to check
        release (str): release label to check

    Returns:
        list: dataset names in map but not on synapse
    """

    query = f'SELECT DISTINCT dataset FROM {config["synapse"]["dataset"]["id"]} WHERE dataset IS NOT NULL'
    table_ds = syn.tableQuery(query).asDataFrame()
    map_ds = df["dataset"].unique()
    res = set([x for x in map_ds if pd.isnull(x) == False]) - set(table_ds["dataset"])
    return list(res)


def check_release_status_ambiguous(
    df: pd.DataFrame, syn: Synapse, config: Dict, cohort: str, release: str
) -> List:
    """Check for any codes with release status that is not y or n.

    Args:
        df (pd.DataFrame): dataframe representing map
        syn (Synapse): Synapse object
        config (dict): configuration parameters
        cohort (str): cohort label to check
        release (str): release label to check

    Returns:
        list: codes with ambiguous release status
    """
    status = df[cohort].str.lower()
    codes = df.loc[[item not in ["y", "n"] for item in status]]["code"]
    return codes


def check_release_status_map_yes_sor_not(self,
    df: pd.DataFrame, syn: Synapse, config: Dict, cohort: str, release: str
) -> List:
    """Check for codes where release status in mapping file is yes
    but relase status in scope of release is not yes.

    Args:
        df (pd.DataFrame): dataframe representing map
        syn (Synapse): Synapse object
        config (dict): configuration parameters
        cohort (str): cohort label to check
        release (str): release label to check

    Returns:
        list: codes with lenient release status
    """
    map_status = df[cohort].str.lower()
    map_type = df["data_type"].str.lower()
    map_codes = df.loc[
        ((map_status == "y") & ((map_type == "derived") | (map_type == "curated")))
    ]["code"]

    column_name = self.get_sor_column_name(syn, config["synapse"]["release"]["id"], cohort, release)
    file_sor = syn.get(config["synapse"]["sor"]["id"])["path"]
    sor = pd.read_excel(file_sor, engine="openpyxl", sheet_name=1)
    sor_status = sor[column_name].str.lower()
    sor_codes = sor.loc[sor_status.isin(["yes", "always"])]["VARNAME"]

    map_not_sor = list(set(map_codes) - set(sor_codes))
    code_remove = get_codes_to_remove(map_not_sor)
    for code in code_remove:
        map_not_sor.remove(code)

    return map_not_sor


def check_release_status_sor_yes_map_not(self,
    df: pd.DataFrame, syn: Synapse, config: Dict, cohort: str, release: str
) -> List:
    """Check for codes where release status in scope of release is yes
    but relase status in mapping file is not yes.

    Args:
        df (pd.DataFrame): dataframe representing map
        syn (Synapse): Synapse object
        config (dict): configuration parameters
        cohort (str): cohort label to check
        release (str): release label to check

    Returns:
        list: codes with lenient release status
    """
    map_status = df[cohort].str.lower()
    map_type = df["data_type"].str.lower()
    map_codes = df.loc[
        ((map_status == "y") & ((map_type == "derived") | (map_type == "curated")))
    ]["code"]

    column_name = self.get_sor_column_name(syn, cohort, release)
    file_sor = syn.get(config["synapse"]["sor"]["id"])["path"]
    sor = pd.read_excel(file_sor, engine="openpyxl", sheet_name=1)
    sor_status = sor[column_name].str.lower()
    sor_codes = sor.loc[sor_status.isin(["yes", "always"])]["VARNAME"]

    inter = set(sor_codes).intersection(set(df["code"]))
    sor_not_map = list(inter - set(map_codes))
    code_remove = get_codes_to_remove(sor_not_map)
    for code in code_remove:
        sor_not_map.remove(code)

    return sor_not_map


def check_code_name_catalog(
    df: pd.DataFrame, syn: Synapse, config: Dict, cohort: str, release: str
) -> List:
    """Check for any variable name not in the Sage data element catalog.

    Args:
        df (pd.DataFrame): dataframe representing map
        syn (Synapse): Synapse object
        config (dict): configuration parameters
        cohort (str): cohort label to check
        release (str): release label to check

    Returns:
        list: dataset names in map but not on catalog
    """

    query = f'SELECT DISTINCT variable FROM {config["synapse"]["catalog"]["id"]}'
    table_var = syn.tableQuery(query).asDataFrame()

    map_type = df["data_type"].str.lower()
    map_nonwild = ["*" not in code for code in df["code"]]
    map_codes = df.loc[
        ((map_nonwild) & ((map_type == "derived") | (map_type == "curated")))
    ]["code"]

    res = set(map_codes) - set(table_var["variable"])
    return list(res)


def format_result(codes: List, config: Dict, check_no: int) -> pd.DataFrame:
    """Format output for interpretable log file.

    Args:
        codes (list): problematic values
        config (dict): configuration parameters
        check_no (int): check number for which to format results

    Returns:
        pd.DataFrame: [description]
    """
    formatted = pd.DataFrame()
    formatted["code"] = codes
    formatted["check_no"] = str(check_no)
    formatted["description"] = config["check"][check_no]["description"]
    formatted["action"] = config["check"][check_no]["request"]
    return formatted


def create_function_map() -> Dict:
    """Create a dictionary that maps function to string representation.

    Returns:
        dict: map from function to name as string
    """
    fxns = {
        "check_code_name_absent": check_code_name_absent,
        "check_code_name_empty": check_code_name_empty,
        "check_dataset_names": check_dataset_names,
        "check_release_status_ambiguous": check_release_status_ambiguous,
        "check_release_status_map_yes_sor_not": check_release_status_map_yes_sor_not,
        "check_release_status_sor_yes_map_not": check_release_status_sor_yes_map_not,
        "check_code_name_catalog": check_code_name_catalog,
    }
    return fxns


def validate_map(
    synapse_id: str,
    file: str,
    syn: Synapse,
    config: Dict,
    version: int,
    cohort: str,
    release: str,
) -> pd.DataFrame:
    """Run all implemented checks on mapping file.

    Args:
        synapse_id (str): Synapse ID of mapping file
        file (str): local path to mapping file
        syn (Synapse): Synapse object
        config (dict): configuration parameters
        version (int):  Version number of Synapse ID
        cohort (str): cohort label to check
        release (str): release label to check

    Returns:
        pd.DataFrame: [description]
    """

    errors = pd.DataFrame()
    df = pd.DataFrame()
    fxns = create_function_map()
    if synapse_id is not None:
        if version == "None":
            df = pd.read_csv(syn.get(synapse_id)["path"])
        else:
            df = pd.read_csv(syn.get(synapse_id, version=version)["path"])
    else:
        df = pd.read_csv(file)

    for check_no in config["check"]:

        logging.info(f"Check {check_no}...")

        if (
            config["check"][check_no]["implemented"]
            and not config["check"][check_no]["deprecated"]
        ):
            fxn_name = config["check"][check_no]["function"]
            result = fxns[fxn_name](df, syn, config, cohort, release)
            errors = errors.append(format_result(result, config, check_no))
            logging.info(f"  Found {len(result)} error(s).")
        else:
            logging.info("  Check deprecated or not implemented.")

    errors.insert(0, "issue", range(1, errors.shape[0] + 1, 1))

    return errors


def build_parser(cohorts: List, releases: List):
    """Build command line parser.

    Args:
        cohorts (List): cohort label options
        releases (List): release label options

    Returns:
        [type]: [description]
    """
    parser = argparse.ArgumentParser(
        description="Checks validity of BPC to cBioPortal mapping file "
    )
    group = parser.add_mutually_exclusive_group()
    group.add_argument(
        "--synapse_id",
        "-s",
        metavar="SYNAPSE_ID",
        type=str,
        help="Synapse ID of mapping file",
    )
    group.add_argument(
        "--file",
        "-f",
        metavar="FILE",
        type=str,
        help="Local path to mapping file",
    )
    parser.add_argument(
        "--version",
        "-v",
        metavar="VERSION",
        type=str,
        default="None",
        help="Synapse entity version number " "(default: current)",
    )
    parser.add_argument(
        "--cohort",
        "-c",
        choices=cohorts,
        default=cohorts[0],
        help="BPC cohort label " "(default: %(default)s)",
    )
    parser.add_argument(
        "--release",
        "-r",
        choices=releases,
        default=releases[0],
        help="Release label " "(default: %(default)s)",
    )
    parser.add_argument(
        "--outfile",
        "-o",
        metavar="OUTFILE",
        type=str,
        default="output.csv",
        help="Name of output file " "(default: %(default)s)",
    )
    parser.add_argument(
        "--log",
        "-l",
        type=str,
        choices=["debug", "info", "warning", "error"],
        default="error",
        help="Set logging output level " "(default: %(default)s)",
    )
    return parser


def get_cohorts(syn: Synapse, config: Dict) -> List:
    """Get sorted list of cohort options.

    Args:
        syn (Synapse): Synapse object
        config (dict): configuration file contents

    Returns:
        list: sorted list of possible cohort labels
    """
    synid_table_rel = config["synapse"]["release"]["id"]
    df = syn.tableQuery(f"SELECT DISTINCT cohort AS cohort FROM {synid_table_rel} ORDER BY cohort").asDataFrame()
    return df["cohort"].values.tolist()


def get_releases(syn: Synapse, config: Dict) -> List:
    """Get sorted list of release options.

    Args:
        syn (Synapse): Synapse object
        config (dict): configuration file contents

    Returns:
        list: sorted list of possible release labels
    """
    synid_table_rel = config["synapse"]["release"]["id"]
    df = syn.tableQuery(f"SELECT DISTINCT CONCAT(release_version, '-', release_type) AS rel FROM {synid_table_rel} ORDER BY rel").asDataFrame()
    return df["rel"].values.tolist()


def read_config(file: str) -> Dict:
    config = None
    with open(file, "r") as stream:
        try:
            config = yaml.safe_load(stream)
        except yaml.YAMLError as exc:
            print(exc)
    return config


def synapse_login(synapse_config=synapseclient.client.CONFIG_FILE):
    """Login to Synapse

    Args:
        synapse_config ([type], optional): Path to synapse configuration file.. Defaults to synapseclient.client.CONFIG_FILE.

    Raises:
        ValueError: error logging into Synapse with given credentials

    Returns:
        Synapse: Synapse object
    """
    try:
        syn = synapseclient.Synapse(skip_checks=True, configPath=synapse_config)
        syn.login(silent=True)
    except (SynapseNoCredentialsError, SynapseAuthenticationError):
        raise ValueError(
            "Login error: please make sure you have correctly "
            "configured your client.  Instructions here: "
            "https://help.synapse.org/docs/Client-Configuration.1985446156.html.  "
            "You can also create a Synapse Personal Access Token and set it "
            "as an environmental variable: "
            "SYNAPSE_AUTH_TOKEN='<my_personal_access_token>'"
        )
    return syn


def main():

    config = read_config("config.yaml")
    args = build_parser(
        cohorts=get_cohorts(config), releases=get_releases(config)
    ).parse_args()
    syn = synapse_login()

    numeric_level = getattr(logging, args.log.upper(), None)
    if not isinstance(numeric_level, int):
        raise ValueError("Invalid log level: %s" % args.log)
    logging.basicConfig(level=numeric_level)

    res = validate_map(
        args.synapse_id, args.file, syn, config, args.version, args.cohort, args.release
    )
    res.to_csv(args.outfile, index=False)

    logging.info(f"Output written to '{args.outfile}'")


if __name__ == "__main__":
    main()
