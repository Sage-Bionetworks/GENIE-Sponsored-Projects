"""BPC configuration classes"""
from .bpc_redcap_export_mapping import BpcProjectRunner


class Brca(BpcProjectRunner):
    """BrCa BPC sponsored project"""

    # Sponsored project name
    _SPONSORED_PROJECT = "BrCa"
    _exclude_files = ["data_timeline_performance_status.txt"]


class Crc(BpcProjectRunner):
    """CRC BPC sponsored project"""

    # Sponsored project name
    _SPONSORED_PROJECT = "CRC"
    _exclude_files = ["data_timeline_performance_status.txt"]


class Nsclc(BpcProjectRunner):
    """NSCLC BPC sponsored project"""

    # Sponsored project name
    _SPONSORED_PROJECT = "NSCLC"
    _exclude_files = ["data_timeline_labtest.txt", "data_timeline_performance_status.txt"]


class Panc(BpcProjectRunner):
    """PANC BPC sponsored project"""

    # Sponsored project name
    _SPONSORED_PROJECT = "PANC"
    _exclude_files = ["data_timeline_performance_status.txt"]


class Prostate(BpcProjectRunner):
    """Prostate BPC sponsored project"""

    # Sponsored project name
    _SPONSORED_PROJECT = "Prostate"
    _exclude_files = ["data_timeline_performance_status.txt"]


class Bladder(BpcProjectRunner):
    """BLADDER BPC sponsored project"""

    # Sponsored project name
    _SPONSORED_PROJECT = "BLADDER"
    _exclude_files = ["data_timeline_labtest.txt"]
