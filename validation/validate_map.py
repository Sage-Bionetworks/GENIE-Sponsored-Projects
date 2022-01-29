"""
Description: Validate the BPC to cBioPortal mapping file. 
Author: Haley Hunter-Zinck
Date: 2022-01-27
"""

import argparse
import logging
import re

import pandas as pd
import synapseclient
from synapseclient import Synapse
from synapseclient.core.exceptions import (
    SynapseAuthenticationError,
    SynapseNoCredentialsError,
)
import yaml


def _check_code_name_empty(df: pd.DataFrame, syn: Synapse, config: dict) -> list:
    """Check for any code that is empty.
     Args:
      df: dataframe representing map
      syn: Synapse object
      config: configuration parameters
    Returns:
        dataframe with metadata on any empty codes.
    """
    empty = df.loc[pd.isna(df["code"])]["code"]
    return list(empty)


def _check_code_name_absent(df: pd.DataFrame, syn: Synapse, config: dict) -> list:
    """Check for any code that is not code name that
    does not appear in its associated data file.
    Args:
        df: dataframe representing map
        syn: Synapse object
        config: configuration parameters
    Returns:
        dataframe with metadata on any missing codes.
    """
    absent = []
    for dataset in config["dataset"]:
        data = pd.read_csv(
            syn.get(config["dataset"][dataset]["id"])["path"], low_memory=False
        )
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
        code_remove = []
        for code in code_map:
            if bool(re.match(r"^.+[*]$", str(code).strip())):
                code_remove.append(code)
            elif pd.isna(code):
                code_remove.append(code)
        for code in code_remove:
            code_map.remove(code)

        absent.extend(list(set(code_map) - set(code_data)))
    return absent


def _format_result(codes: list, config: dict, check_no: int) -> pd.DataFrame:
    """Format output for interpretable log file.
    Args:
        df: dataframe representing map
        config: configuration parameters
        check_no: check number for which to format results
    Returns:
        dataframe with additional metadata on any errors.
    """
    formatted = pd.DataFrame()
    formatted["code"] = codes
    formatted["check_no"] = str(check_no)
    formatted["description"] = config["check"][check_no]["description"]
    formatted["action"] = config["check"][check_no]["request"]
    return formatted


def _create_function_map() -> dict:
  fxns = {
  "_check_code_name_absent": _check_code_name_absent,
  "_check_code_name_empty": _check_code_name_empty
  }
  return fxns


def validate_map(
    synapse_id: str, syn: Synapse, config: dict, version: int
) -> pd.DataFrame:
    """Run all implemented checks on mapping file.
    Args:
        synapse_id: Synapse ID of mapping file
        syn: Synapse object
        config: configuration parameters
        version: Version number of Synapse ID
    Returns:
        dataframe with additional metadata on any errors.
    """

    errors = pd.DataFrame()
    df = pd.DataFrame()
    fxns = _create_function_map()
    if version == "None":
        df = pd.read_csv(syn.get(synapse_id)["path"])
    else:
        df = pd.read_csv(syn.get(synapse_id, version=version)["path"])

    for check_no in config["check"]:

        logging.info(f"Check {check_no}...")

        if (
            config["check"][check_no]["implemented"]
            and not config["check"][check_no]["deprecated"]
        ):
            fxn_name = config["check"][check_no]["function"]
            result = fxns[fxn_name](df, syn, config)
            errors = errors.append(_format_result(result, config, check_no))
            logging.info(f"  Found {errors.shape[0]} error(s).")
        else:
            logging.info("  Check deprecated or not implemented.")

    errors.insert(0, "issue", range(1, errors.shape[0] + 1, 1))

    return errors


def build_parser():
    parser = argparse.ArgumentParser(
        description="Checks validity of BPC to cBioPortal mapping file "
    )
    parser.add_argument(
        "synapse_id",
        metavar="SYNAPSE_ID",
        type=str,
        help="Synapse ID of mapping file",
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


def read_config(file: str) -> dict:
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
        synapse_config: Path to synapse configuration file.
                        Defaults to ~/.synapseConfig
    Returns:
        Synapse connection
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

    args = build_parser().parse_args()
    config = read_config("config.yaml")
    syn = synapse_login()

    numeric_level = getattr(logging, args.log.upper(), None)
    if not isinstance(numeric_level, int):
        raise ValueError("Invalid log level: %s" % args.log)
    logging.basicConfig(level=numeric_level)

    res = validate_map(args.synapse_id, syn, config, args.version)
    res.to_csv(args.outfile, index=False)


if __name__ == "__main__":
    main()
