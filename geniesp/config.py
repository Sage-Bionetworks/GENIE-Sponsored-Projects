from dataclasses import dataclass, field
import subprocess
import os
from typing import List


def get_git_sha() -> str:
    """get git sha digest"""
    text = subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True, text=True)
    return text.stdout.rstrip("\n")


@dataclass
class BpcConfig:
    # json_path: str
    cohort: str
    # Redcap codes to cbioportal mapping synid and form key is in
    # version 38, 42 were last stable version(s)
    redcap_to_cbio_mapping_synid: str = "syn25712693.49"
    # Run `git rev-parse HEAD` in Genie_processing directory to obtain shadigest
    # github_repo = None
    # Mapping from Synapse Table to derived variables
    # TODO: Make versioned
    data_tables_id: str = "syn22296821"
    # Storage of not found samples
    sp_redcap_exports_synid: str = "syn21446571"
    # main GENIE release folder (14.8-consortium)
    # Must use consortium release, because SEQ_YEAR is used
    mg_release_synid: str = "syn52794528"
    # PRISSMM documentation table
    prissmm_synid: str = "syn22684834"
    # BPC sample retraction table
    sample_retraction_synid: str = "syn25779833"
    patient_retraction_synid: str = "syn25998970"
    retraction_at_release_synid: str = "syn52915299"
    temporary_patient_retraction_synid: str = "syn29266682"
    # main GENIE assay information table
    mg_assay_synid: str = "syn17009222"
    # staging release folder to upload cbioportal exports
    staging_release_folder = "syn50876969"
    # exclude files to be created for cbioportal
    # TODO: need to support this feature in rest of code, for now
    # This is added for metadata files
    exclude_files: List[str] = field(default_factory=list)
    # cohort-generic link to documentation for BPC datasets
    url_bpc: str = "https://aacr.box.com/s/en5dyu9zfw1krg2u58wlcz01jttc6y9h"
    # cohort-generic link to documentation for cBio files
    url_cbio: str = "https://docs.google.com/document/d/1IBVF-FLecUG8Od6mSEhYfWH3wATLNMnZcBw2_G0jSAo/edit"
    # syn: Synapse
    oncotreelink: str = "https://oncotree.info/api/tumorTypes/tree?version=oncotree_2021_11_02"
    github_url: str = f"https://github.com/Sage-Bionetworks/GENIE-Sponsored-Projects/tree/{get_git_sha()}"

    def __post_init__(self):
        if self.cohort is None:
            raise ValueError("cohort must be specified")
        if not os.path.exists(self.cohort):
            os.mkdir(self.cohort)
        else:
            import glob
            filelists = glob.glob("**/*", recursive=True)
            for each_file in filelists:
                os.remove(os.path.join(self.cohort, each_file))


class Brca(BpcConfig):
    """BrCa BPC sponsored project"""

    # Sponsored project name
    cohort = "BrCa"
    exclude_files = ["data_timeline_performance_status.txt"]


class Crc(BpcConfig):
    """CRC BPC sponsored project"""

    # Sponsored project name
    cohort = "CRC"
    exclude_files = ["data_timeline_performance_status.txt"]


class Nsclc(BpcConfig):
    """NSCLC BPC sponsored project"""

    # Sponsored project name
    cohort = "NSCLC"
    exclude_files = ["data_timeline_labtest.txt", "data_timeline_performance_status.txt"]


class Panc(BpcConfig):
    """PANC BPC sponsored project"""

    # Sponsored project name
    cohort = "PANC"
    exclude_files = ["data_timeline_performance_status.txt"]


class Prostate(BpcConfig):
    """Prostate BPC sponsored project"""

    # Sponsored project name
    cohort = "Prostate"
    exclude_files = ["data_timeline_performance_status.txt"]


class Bladder(BpcConfig):
    """BLADDER BPC sponsored project"""

    # Sponsored project name
    cohort = "BLADDER"
    exclude_files = ["data_timeline_labtest.txt"]
