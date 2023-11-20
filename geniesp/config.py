from dataclasses import dataclass, field
import subprocess
from typing import List

from synapseclient import Synapse


def get_git_sha() -> str:
    """get git sha digest"""
    text = subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True, text=True)
    return text.stdout.rstrip("\n")


def get_latest_entity_version(syn: Synapse, synid: str) -> str:
    """Get latest version of an entity

    Args:
        syn (Synapse): Synapse connection
        synid (str): Synapse Id

    Returns:
        dict:
            id: Synapse Id
            versionNumber: version of entity
    """
    versions = syn.restGET(f"/entity/{synid}/version")
    return versions['results'][0]


@dataclass
class BpcConfig:
    # json_path: str
    cohort: str
    syn: Synapse
    # Redcap codes to cbioportal mapping synid and form key is in
    # version 49 were last stable version(s)
    redcap_to_cbio_mapping_synid: str = "syn25712693"
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
    oncotreelink: str = "https://oncotree.info/api/tumorTypes/tree?version=oncotree_2021_11_02"
    github_url: str = f"https://github.com/Sage-Bionetworks/GENIE-Sponsored-Projects/tree/{get_git_sha()}"

    def __post_init__(self):
        redcap_to_cbio_mapping = get_latest_entity_version(
            syn=self.syn,
            synid=self.redcap_to_cbio_mapping_synid
        )
        self.redcap_to_cbio_mapping_synid = '{}.{}'.format(
            redcap_to_cbio_mapping['id'],
            redcap_to_cbio_mapping['versionNumber']
        )
        mg_assay = get_latest_entity_version(
            syn=self.syn,
            synid=self.mg_assay_synid
        )
        self.mg_assay_synid = '{}.{}'.format(
            mg_assay['id'],
            mg_assay['versionNumber']
        )

    def to_dict(self) -> dict:
        """Return configuration used

        Returns:
            dict: Configuration
        """
        return {
            "cohort": self.cohort,
            "redcap_to_cbio_mapping_synid": self.redcap_to_cbio_mapping_synid,
            "data_tables_id": self.data_tables_id,
            "sp_redcap_exports_synid": self.sp_redcap_exports_synid,
            "mg_release_synid": self.mg_release_synid,
            "prissmm_synid": self.prissmm_synid,
            "sample_retraction_synid": self.sample_retraction_synid,
            "patient_retraction_synid": self.patient_retraction_synid,
            "retraction_at_release_synid": self.retraction_at_release_synid,
            "temporary_patient_retraction_synid": self.temporary_patient_retraction_synid,
            "mg_assay_synid": self.mg_assay_synid,
            "staging_release_folder": self.staging_release_folder,
            "exclude_files": self.exclude_files,
            "url_bpc": self.url_bpc,
            "url_cbio": self.url_cbio,
            "oncotreelink": self.oncotreelink,
            "github_url": self.github_url
        }


@dataclass
class Brca(BpcConfig):
    """BrCa BPC sponsored project"""

    # Sponsored project name
    cohort = "BrCa"
    exclude_files = ["data_timeline_performance_status.txt"]


@dataclass
class Crc(BpcConfig):
    """CRC BPC sponsored project"""

    # Sponsored project name
    cohort = "CRC"
    exclude_files = ["data_timeline_performance_status.txt"]


@dataclass
class Nsclc(BpcConfig):
    """NSCLC BPC sponsored project"""

    # Sponsored project name
    cohort = "NSCLC"
    exclude_files = ["data_timeline_labtest.txt", "data_timeline_performance_status.txt"]


@dataclass
class Panc(BpcConfig):
    """PANC BPC sponsored project"""

    # Sponsored project name
    cohort = "PANC"
    exclude_files = ["data_timeline_performance_status.txt"]


@dataclass
class Prostate(BpcConfig):
    """Prostate BPC sponsored project"""

    # Sponsored project name
    cohort = "Prostate"
    exclude_files = ["data_timeline_performance_status.txt"]


@dataclass
class Bladder(BpcConfig):
    """BLADDER BPC sponsored project"""

    # Sponsored project name
    cohort = "BLADDER"
    exclude_files = ["data_timeline_labtest.txt"]
