"""GENIE SP/BPC cBioPortal exporter CLI"""
import argparse
import logging

import synapseclient

from .bpc_config import Brca, Crc, Nsclc, Panc, Prostate, Bladder
from .sp_config import Akt1, Erbb2, Fgfr4

BPC_MAPPING = {
    "NSCLC": Nsclc,
    "CRC": Crc,
    "BrCa": Brca,
    "PANC": Panc,
    "Prostate": Prostate,
    "BLADDER": Bladder,
    "AKT1": Akt1,
    "ERRB2": Erbb2,
    "FGFR4": Fgfr4,
}


def main():
    """Main"""
    parser = argparse.ArgumentParser(description="Run GENIE sponsored projects")
    parser.add_argument(
        "sp",
        type=str,
        help="Specify sponsored project to run",
        choices=BPC_MAPPING.keys(),
    )

    parser.add_argument("release", type=str, help="Specify bpc release")
    parser.add_argument(
        "--upload",
        action="store_true",
        help="Upload files into Synapse BPC staging directory. Default: false.",
    )
    parser.add_argument(
        "--log",
        "-l",
        type=str,
        choices=["debug", "info", "warning", "error"],
        default="info",
        help="Set logging output level " "(default: %(default)s)",
    )
    parser.add_argument(
        "--cbioportal",
        type=str,
        help="Optional parameter to specify cbioportal folder location",
    )
    args = parser.parse_args()

    numeric_level = getattr(logging, args.log.upper(), None)
    if not isinstance(numeric_level, int):
        raise ValueError("Invalid log level: %s" % args.log)
    logging.basicConfig(level=numeric_level)

    syn = synapseclient.login()

    if args.cbioportal is None:
        cbiopath = "../cbioportal"
    else:
        cbiopath = args.cbioportal

    BPC_MAPPING[args.sp](
        syn, cbiopath, release=args.release, upload=args.upload
    ).run()


if __name__ == "__main__":
    main()
