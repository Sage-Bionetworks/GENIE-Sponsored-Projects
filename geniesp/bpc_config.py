"""
BPC configuration classes
>>> git clone https://github.com/cBioPortal/cbioportal.git
>>> python run_bpc.py NSCLC ../../cbioportal 1.1-consortium --staging
"""
from bpc_redcap_export_mapping import BpcProjectRunner


class Brca(BpcProjectRunner):
    """NSCLC BPC sponsored project"""
    # Sponsorted project name
    _SPONSORED_PROJECT = 'BrCa'
    # Redcap codes to cbioportal mapping synid and form key is in
    _REDCAP_TO_CBIOMAPPING_SYNID = "syn25712693.3"
    # Mapping from Synapse Table to form (derived files)
    _DATA_TABLE_IDS = "syn22296821"
    # Storage of not found samples
    _SP_REDCAP_EXPORTS_SYNID = "syn21446571"
    # Run `git rev-parse HEAD` in Genie_processing directory to
    # obtain shadigest
    _GITHUB_REPO = (
        "https://github.com/Sage-Bionetworks/GENIE-Sponsored-Projects/"
        # "tree/a672a70ad5195e6e8359325f5cea10bef384b2ff/"
        # "geniesp/bpc_config.py"
    )


class Crc(BpcProjectRunner):
    """NSCLC BPC sponsored project"""
    # Sponsorted project name
    _SPONSORED_PROJECT = 'CRC'
    # Redcap codes to cbioportal mapping synid and form key is in
    _REDCAP_TO_CBIOMAPPING_SYNID = "syn25712693.3"
    # Mapping from Synapse Table to form (derived files)
    # TODO: Make versioned
    _DATA_TABLE_IDS = "syn22296821"
    # Storage of not found samples
    _SP_REDCAP_EXPORTS_SYNID = "syn21446571"
    # Run `git rev-parse HEAD` in Genie_processing directory to
    # obtain shadigest
    _GITHUB_REPO = (
        "https://github.com/Sage-Bionetworks/GENIE-Sponsored-Projects/"
        # "tree/765a209402a0e4c8517ec826ddad1f05d842f54a/"
        # "geniesp/bpc_config.py"
    )


class Nsclc(BpcProjectRunner):
    """NSCLC BPC sponsored project"""
    # Sponsorted project name
    _SPONSORED_PROJECT = 'NSCLC'
    # Redcap codes to cbioportal mapping synid and form key is in
    _REDCAP_TO_CBIOMAPPING_SYNID = "syn25712693.3"
    # Mapping from Synapse Table to form (derived files)
    _DATA_TABLE_IDS = "syn22296821"
    # Storage of not found samples
    _SP_REDCAP_EXPORTS_SYNID = "syn21446571"
    # Run `git rev-parse HEAD` in Genie_processing directory to
    # obtain shadigest
    _GITHUB_REPO = (
        "https://github.com/Sage-Bionetworks/GENIE-Sponsored-Projects/"
        # "tree/765a209402a0e4c8517ec826ddad1f05d842f54a/"
        # "geniesp/bpc_config.py"
    )
