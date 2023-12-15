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
    """Configuration class for BPC sponsored projects
    that will be leveraged throughout the code. All Synapse Ids
    and configurations should be stored in this class.

    Attributes:
        cohort: BPC cohort (e.g. BLADDER, NSCLC, etc)
        syn: Synapse connection
        redcap_to_cbio_mapping_synid: Redcap to cbioportal mapping. Defaults to syn25712693.
                                      The latest version will be used.
        data_tables_id: Synapse file view containing derived variable files. Defaults to syn22296821.
        mg_release_synid: Synapse folder containing main GENIE release. Defaults to syn52794528
                          (14.8-consortium). Must use consortium release, because SEQ_YEAR is used.
        prissmm_synid: PRISSMM documentation table. Defaults to syn22684834.
        sample_retraction_synid: BPC sample retraction table. Defaults to syn25779833.
        patient_retraction_synid: BPC patient retraction table. Defaults to syn25998970.
        retraction_at_release_synid: BPC retraction at release table. Defaults to syn52915299.
        temporary_patient_retraction_synid: BPC temporary patient retraction table. Defaults to syn29266682.
        mg_assay_synid: Main GENIE assay information table. Defaults to syn17009222.
                        The latest version will be used.
        staging_release_folder: Staging release folder to upload cbioportal exports. Defaults to syn50876969.
        exclude_files: Exclude files to be created for cbioportal.
        url_bpc: Cohort-generic link to documentation for BPC datasets. Defaults to
                 https://aacr.box.com/s/en5dyu9zfw1krg2u58wlcz01jttc6y9h.
        url_cbio: Cohort-generic link to documentation for cBio files. Defaults to
                  https://docs.google.com/document/d/1IBVF-FLecUG8Od6mSEhYfWH3wATLNMnZcBw2_G0jSAo/edit.
        oncotreelink: Link to oncotree api. Defaults to https://oncotree.info/api/tumorTypes/tree?version=oncotree_2021_11_02.
        github_url: Link to github repo. Defaults to https://github.com/Sage-Bionetworks/GENIE-Sponsored-Projects.
    """
    cohort: str
    """BPC cohort (e.g. BLADDER, NSCLC, etc)"""
    syn: Synapse
    """Synapse connection"""
    redcap_to_cbio_mapping_synid: str = "syn25712693"
    """Redcap codes to cbioportal mapping synid and form key is in
    version 49 were last stable version(s)"""
    # TODO: Make versioned
    data_tables_id: str = "syn22296821"
    """Mapping from Synapse Table to derived variables"""
    mg_release_synid: str = "syn52794528"
    """Main GENIE release folder (14.8-consortium)
    Must use consortium release, because SEQ_YEAR is used"""
    prissmm_synid: str = "syn22684834"
    """PRISSMM documentation table"""
    sample_retraction_synid: str = "syn25779833"
    """BPC sample retraction table"""
    patient_retraction_synid: str = "syn25998970"
    """BPC patient retraction table"""
    retraction_at_release_synid: str = "syn52915299"
    """Retract at release patients"""
    temporary_patient_retraction_synid: str = "syn29266682"
    """Temporary patient retraction table"""
    mg_assay_synid: str = "syn17009222"
    """main GENIE assay information table"""
    staging_release_folder = "syn50876969"
    """staging release folder to upload cbioportal exports"""
    exclude_files: List[str] = field(default_factory=list)
    """exclude files to be created for cbioportal. This is added for metadata files"""
    url_bpc: str = "https://aacr.box.com/s/en5dyu9zfw1krg2u58wlcz01jttc6y9h"
    """cohort-generic link to documentation for BPC datasets"""
    url_cbio: str = "https://docs.google.com/document/d/1IBVF-FLecUG8Od6mSEhYfWH3wATLNMnZcBw2_G0jSAo/edit"
    """cohort-generic link to documentation for cBio files"""
    oncotreelink: str = "https://oncotree.info/api/tumorTypes/tree?version=oncotree_2021_11_02"
    """Oncotree Link"""
    github_url: str = "https://github.com/Sage-Bionetworks/GENIE-Sponsored-Projects"
    """GitHub URL to sponsored project repo"""

    def __post_init__(self):
        """Post initialization"""
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
        self.github_url = f"https://github.com/Sage-Bionetworks/GENIE-Sponsored-Projects/tree/{get_git_sha()}"


    def to_dict(self) -> dict:
        """Return configuration used

        Returns:
            dict: Configuration
        """
        return {
            "cohort": self.cohort,
            "redcap_to_cbio_mapping_synid": self.redcap_to_cbio_mapping_synid,
            "data_tables_id": self.data_tables_id,
            # "sp_redcap_exports_synid": self.sp_redcap_exports_synid,
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
    cohort = "BrCa"
    exclude_files = ["data_timeline_performance_status.txt"]


@dataclass
class Crc(BpcConfig):
    """CRC BPC sponsored project"""
    cohort = "CRC"
    exclude_files = ["data_timeline_performance_status.txt"]


@dataclass
class Nsclc(BpcConfig):
    """NSCLC BPC sponsored project"""
    cohort = "NSCLC"
    exclude_files = ["data_timeline_labtest.txt", "data_timeline_performance_status.txt"]


@dataclass
class Panc(BpcConfig):
    """PANC BPC sponsored project"""
    cohort = "PANC"
    exclude_files = ["data_timeline_performance_status.txt"]


@dataclass
class Prostate(BpcConfig):
    """Prostate BPC sponsored project"""
    cohort = "Prostate"
    exclude_files = ["data_timeline_performance_status.txt"]


@dataclass
class Bladder(BpcConfig):
    """BLADDER BPC sponsored project"""
    cohort = "BLADDER"
    exclude_files = ["data_timeline_labtest.txt"]
