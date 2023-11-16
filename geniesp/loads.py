from dataclasses import dataclass
from functools import cached_property
from typing import Union

import pandas as pd

import synapseclient
from synapseclient import Folder, Synapse
from genie import process_functions

from geniesp.config import BpcConfig


@dataclass
class Loads():
    bpc_config: BpcConfig
    upload: bool = False
    syn: Synapse

    @cached_property
    def cbioportal_folders(self) -> dict:
        """Create case lists and release folder"""
        if self.upload:
            return {}
        sp_data_folder = self.syn.store(
            Folder(self.bpc_config.cohort, parentId=self.bpc_config.staging_release_folder)
        )
        release_folder = self.syn.store(Folder(self.release, parent=sp_data_folder))
        # Store in cBioPortal files because there is another folder for
        # clinical files
        release_folder = self.syn.store(
            Folder("cBioPortal_files", parent=release_folder)
        )
        case_lists = self.syn.store(Folder("case_lists", parent=release_folder))
        return {"release": release_folder, "case_lists": case_lists}

    def store(self, parent: str, filepath: str, used_entities: Union[list, None] = None):
        ent = synapseclient.File(filepath, parent=parent)
        self.syn.store(
            ent,
            executed="https://github.com/Sage-Bionetworks/GENIE-Sponsored-Projects",
            used=used_entities
        )
