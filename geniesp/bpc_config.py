"""BPC configuration classes"""
from .bpc_redcap_export_mapping import BpcProjectRunner


class Brca(BpcProjectRunner):
    """NSCLC BPC sponsored project"""

    # Sponsored project name
    _SPONSORED_PROJECT = "BrCa"


class Crc(BpcProjectRunner):
    """NSCLC BPC sponsored project"""

    # Sponsored project name
    _SPONSORED_PROJECT = "CRC"


class Nsclc(BpcProjectRunner):
    """NSCLC BPC sponsored project"""

    # Sponsored project name
    _SPONSORED_PROJECT = "NSCLC"
    _exclude_files = ["data_timeline_labtest.txt"]


class Panc(BpcProjectRunner):
    """PANC BPC sponsored project"""

    # Sponsored project name
    _SPONSORED_PROJECT = "PANC"


class Prostate(BpcProjectRunner):
    """Prostate BPC sponsored project"""

    # Sponsored project name
    _SPONSORED_PROJECT = "Prostate"


class Bladder(BpcProjectRunner):
    """BLADDER BPC sponsored project"""

    # Sponsored project name
    _SPONSORED_PROJECT = "BLADDER"
    _exclude_files = ["data_timeline_labtest.txt"]
