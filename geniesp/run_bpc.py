"""Run BPC projects

>>> python run_bpc.py -h
"""
import argparse
import os

import synapseclient

from bpc_config import Brca, Crc, Nsclc

BPC_MAPPING = {"NSCLC": Nsclc,
               'CRC': Crc,
               'BrCa': Brca}


def main():
    """Main"""
    parser = argparse.ArgumentParser(
        description='Run GENIE sponsored projects'
    )
    parser.add_argument("sp", type=str,
                        help='Specify sponsored project to run',
                        choices=BPC_MAPPING.keys())
    parser.add_argument(
        "cBioPath", type=str,
        help='Specify path to cbio: must do '
             '`git clone https://github.com/cBioPortal/cbioportal.git`'
    )
    parser.add_argument("release", type=str,
                        help='Specify bpc release')
    parser.add_argument("--staging", action='store_true',
                        help="If true, files aren't uploaded onto synapse")
    args = parser.parse_args()
    syn = synapseclient.login()

    BPC_MAPPING[args.sp](syn, args.cBioPath, release=args.release,
                         staging=args.staging).run()


if __name__ == '__main__':
    main()
