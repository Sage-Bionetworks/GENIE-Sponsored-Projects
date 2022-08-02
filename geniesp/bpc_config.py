"""
BPC configuration classes
>>> git clone https://github.com/cBioPortal/cbioportal.git
>>> python run_bpc.py NSCLC ../../cbioportal 1.1-consortium --staging
"""
from .bpc_redcap_export_mapping import BpcProjectRunner


class Brca(BpcProjectRunner):
    """NSCLC BPC sponsored project"""

    # Sponsored project name
    _SPONSORED_PROJECT = "BrCa"
    # Redcap codes to cbioportal mapping synid and form key is in
    _REDCAP_TO_CBIOMAPPING_SYNID = "syn25712693.33"
    # Mapping from Synapse Table to form (derived files)
    _DATA_TABLE_IDS = "syn22296821"
    # Storage of not found samples
    _SP_REDCAP_EXPORTS_SYNID = "syn21446571"
    # main GENIE release folder (11.0-public)
    _MG_RELEASE_SYNID = "syn26706564"


class Crc(BpcProjectRunner):
    """NSCLC BPC sponsored project"""

    # Sponsored project name
    _SPONSORED_PROJECT = "CRC"
    # Redcap codes to cbioportal mapping synid and form key is in
    _REDCAP_TO_CBIOMAPPING_SYNID = "syn25712693.33"
    # Mapping from Synapse Table to form (derived files)
    # TODO: Make versioned
    _DATA_TABLE_IDS = "syn22296821"
    # Storage of not found samples
    _SP_REDCAP_EXPORTS_SYNID = "syn21446571"
    # main GENIE release folder (11.0-public)
    _MG_RELEASE_SYNID = "syn26706564"


class Nsclc(BpcProjectRunner):
    """NSCLC BPC sponsored project"""

    # Sponsored project name
    _SPONSORED_PROJECT = "NSCLC"
    # Redcap codes to cbioportal mapping synid and form key is in
    _REDCAP_TO_CBIOMAPPING_SYNID = "syn25712693.33"
    # Mapping from Synapse Table to form (derived files)
    _DATA_TABLE_IDS = "syn22296821"
    # Storage of not found samples
    _SP_REDCAP_EXPORTS_SYNID = "syn21446571"
    # main GENIE release folder (11.0-public)
    _MG_RELEASE_SYNID = "syn26706564"
    _exclude_files = ["data_timeline_labtest.txt"]


class Panc(BpcProjectRunner):
    """PANC BPC sponsored project"""

    # Sponsored project name
    _SPONSORED_PROJECT = "PANC"
    # Redcap codes to cbioportal mapping synid and form key is in
    _REDCAP_TO_CBIOMAPPING_SYNID = "syn25712693.33"
    # Mapping from Synapse Table to form (derived files)
    _DATA_TABLE_IDS = "syn22296821"
    # Storage of not found samples
    _SP_REDCAP_EXPORTS_SYNID = "syn21446571"
    # main GENIE release folder (11.0-public)
    _MG_RELEASE_SYNID = "syn26706564"


class Prostate(BpcProjectRunner):
    """Prostate BPC sponsored project"""

    # Sponsored project name
    _SPONSORED_PROJECT = "Prostate"
    # Redcap codes to cbioportal mapping synid and form key is in
    _REDCAP_TO_CBIOMAPPING_SYNID = "syn25712693.33"
    # Mapping from Synapse Table to form (derived files)
    _DATA_TABLE_IDS = "syn22296821"
    # Storage of not found samples
    _SP_REDCAP_EXPORTS_SYNID = "syn21446571"
    # main GENIE release folder (11.0-public)
    _MG_RELEASE_SYNID = "syn26706564"


class Bladder(BpcProjectRunner):
    """BLADDER BPC sponsored project"""

    # Sponsored project name
    _SPONSORED_PROJECT = "BLADDER"
    # Redcap codes to cbioportal mapping synid and form key is in
    _REDCAP_TO_CBIOMAPPING_SYNID = "syn25712693.33"
    # Mapping from Synapse Table to form (derived files)
    _DATA_TABLE_IDS = "syn22296821"
    # Storage of not found samples
    _SP_REDCAP_EXPORTS_SYNID = "syn21446571"
    # main GENIE release folder (11.0-public)
    _MG_RELEASE_SYNID = "syn26706564"
    _exclude_files = ["data_timeline_labtest.txt"]
