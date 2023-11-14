from dataclasses import dataclass


@dataclass
class BpcConfig():
    # json_path: str
    cohort = "BLADDER"
    # Redcap codes to cbioportal mapping synid and form key is in
    # version 38, 42 were last stable version(s)
    redcap_to_cbio_mapping_synid = "syn25712693.49"
    # Run `git rev-parse HEAD` in Genie_processing directory to obtain shadigest
    # github_repo = None
    # Mapping from Synapse Table to derived variables
    # TODO: Make versioned
    data_tables_id = "syn22296821"
    # Storage of not found samples
    sp_redcap_exports_synid = "syn21446571"
    # main GENIE release folder (14.8-consortium)
    # Must use consortium release, because SEQ_YEAR is used
    mg_release_synid = "syn52794528"
    # PRISSMM documentation table
    prissmm_synid = "syn22684834"
    # BPC sample retraction table
    sample_retraction_synid = "syn25779833"
    patient_retraction_synid = "syn25998970"
    retraction_at_release_synid = "syn52915299"
    temporary_patient_retraction_synid = "syn29266682"
    # main GENIE assay information table
    mg_assay_synid = "syn17009222"
    # exclude files to be created for cbioportal
    # TODO: need to support this feature in rest of code, for now
    # This is added for metadata files
    exclude_files = []
    # cohort-generic link to documentation for BPC datasets
    url_bpc = "https://aacr.box.com/s/en5dyu9zfw1krg2u58wlcz01jttc6y9h"
    # cohort-generic link to documentation for cBio files
    url_cbio = "https://docs.google.com/document/d/1IBVF-FLecUG8Od6mSEhYfWH3wATLNMnZcBw2_G0jSAo/edit"
    # syn: Synapse
