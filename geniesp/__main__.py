"""GENIE SP/BPC cBioPortal exporter CLI"""
import argparse

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
        "--staging",
        action="store_true",
        help="If true, files aren't uploaded onto synapse",
    )
    parser.add_argument(
        "--cbioportal",
        type=str,
        help="Optional parameter to specify cbioportal folder location",
    )
    args = parser.parse_args()

    syn = synapseclient.login()

    if args.cbioportal is None:
        cbiopath = "../cbioportal"
    else:
        cbiopath = args.cbioportal

    BPC_MAPPING[args.sp](
        syn, cbiopath, release=args.release, staging=args.staging
    ).run()


if __name__ == "__main__":
    main()
