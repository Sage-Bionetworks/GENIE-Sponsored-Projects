"""BPC redcap export PROCESSES
- timeline file
- clinical file
  OS_MONTHS = death_date_int - date_first_met_int
  OS_MONTHS_PRIMARY = death_date_int - primary_dx_date_int
  All dates are converted from days to months (days/30.4)
  Add headers
  REMOVE PATIENTS/SAMPLES THAT DON'T HAVE GENIE SAMPLE IDS

USAGE:
>>> python runSP.py NSCLC ../cbioportal/ --staging
"""
from abc import ABCMeta, abstractmethod
import math
import os
import subprocess

import synapseclient
from synapseclient import File, Folder
import pandas as pd

import genie
from genie import create_case_lists, process_functions


def get_data(syn, mappingdf, sampletype):
    """Extracts the sample, patient and timeline df

    Args:
        syn: Synapse connection
        mappingdf: Mapping dataframe
        sampletype: sample type

    Returns:
        df
    """
    synids = mappingdf['id'][mappingdf['sampleType'] == sampletype].unique()
    finaldf = pd.DataFrame()
    for synid in synids:
        table = syn.tableQuery(f"select * from {synid}")
        tabledf = table.asDataFrame()
        if finaldf.empty:
            finaldf = finaldf.append(tabledf)
        else:
            # must remove this or else columns will be duplicated
            del tabledf['redcap_data_access_group']
            finaldf = finaldf.merge(tabledf, on="record_id")
    # Subset needed columns
    cols = mappingdf['code'][mappingdf['sampleType'] == sampletype]
    cols = cols.tolist()
    cols.extend(['redcap_data_access_group', 'record_id'])
    finaldf = finaldf[cols]
    return finaldf


def get_file_data(syn, mappingdf, sampletype, cohort='NSCLC'):
    """Extracts the sample, patient and timeline df

    Args:
        syn: Synapse connection
        mappingdf: Mapping dataframe
        sampletype: sample type

    Returns:
        dict: mapped dataframe,
              list of used entities
    """
    # Group by dataset because different datasets could have the
    # same variable
    datasets = mappingdf.groupby("dataset")
    finaldf = pd.DataFrame()
    used_entities = []

    for _, df in datasets:
        # Get synapse id
        synid = df['id'].unique()[0]
        table = syn.get(synid)
        used_entities.append(f'{synid}.{table.versionNumber}')
        # obtain columns to subset df
        cols = df['code'][df['sampleType'] == sampletype]
        cols = cols.tolist()
        if "record_id" not in cols:
            cols.append("record_id")
        # Must add path_rep_number for sample and sample acquisition file
        if sampletype in ["SAMPLE", "TIMELINE-SAMPLE"]:
            cols.append("path_rep_number")
        # Must add path_proc_number to sample file
        if sampletype == "SAMPLE":
            cols.append("path_proc_number")
        # Only get specific cohort and subset cols
        tabledf = pd.read_csv(table.path)
        tabledf = tabledf[tabledf['cohort'] == cohort]
        tabledf = tabledf[cols]
        # Append to final dataframe if empty
        if finaldf.empty:
            finaldf = finaldf.append(tabledf)
        else:
            # Records missing pathology reports still have to be present
            # So a left merge has to happen.  This logic also assumes that
            # The pathology report dataset mapping info isn't first.
            # This also assumes that the TIMELINE-PATHOLOGY file only
            # uses columns from the pathology-report dataset
            if df['dataset'][0] == "Pathology-report level dataset":
                finaldf = finaldf.merge(
                    tabledf,
                    on=["record_id", "path_proc_number", "path_rep_number"],
                    how="left"
                )
                del finaldf['path_rep_number']
            else:
                finaldf = finaldf.merge(tabledf, on="record_id")
    # This is to map values (sample acquisition, samples)
    # duplicated_codes_idx = mappingdf['code'].duplicated()
    # if duplicated_codes_idx.any():
    #     duplicated = mappingdf['code'][duplicated_codes_idx].values[0]
    #     finaldf = finaldf[finaldf[f'{duplicated}_x'] == finaldf[f'{duplicated}_y']]
    #     del finaldf[f'{duplicated}_x']
    #     finaldf[duplicated] = finaldf[f'{duplicated}_y']
    #     del finaldf[f'{duplicated}_y']

    return {"df": finaldf,
            "used": used_entities}


# def _remap_os_pfs_values(clinicaldf):
#     """Remap numerical values to string values
#     0 -> 0:LIVING
#     1 -> 1:DECEASED
#     """
#     os_pfs_cols = [col for col in clinicaldf.columns
#                    if col.startswith(('OS', 'PFS')) and
#                    col.endswith("STATUS")]
#     remap_os_values = {col: {0: "0:LIVING", 1: "1:DECEASED"}
#                        for col in os_pfs_cols}
#     return clinicaldf.replace(remap_os_values)


def fill_cancer_dx_start_date(finaldf):
    """Fills in cancer dx start date for those missing pathology dates

    Args:
        finaldf: Mapped cancer diagnosis timeline dataframe

    Returns:
        dataframe with filled START_DATEs
    """
    # Get all time0 points for all records
    time0_dates_idx = finaldf['INDEX_CANCER'] == "Yes"
    # Make all time0 points 0
    finaldf['diagnosis_int'] = finaldf['START_DATE']
    finaldf['START_DATE'][time0_dates_idx] = 0
    # finaldf['START_DATE'][~time0_dates_idx] = 1

    # records = finaldf['PATIENT_ID'].unique()
    # for record in records:
    #     df = finaldf[finaldf['PATIENT_ID'] == record]
    #     # Grab the time 0 START_DATE and calculate the other start times
    #     time0_idx = df['START_DATE'] == 0
    #     values = df['CA_DOB_DX_INT'][time0_idx].unique()
    #     # Sometimes there are more than one index cancers, so must take
    #     # the first value which is the minimum value
    #     time0 = min(values)
    #     # time0 = df['CA_DOB_DX_INT'][time0_idx].values[0]
    #     finaldf.loc[df.index[~time0_idx],
    #                 'START_DATE'] = df['CA_DOB_DX_INT'][~time0_idx]-time0
    #     missing_idx = df['CA_DOB_DX_INT'].isnull()
    #     finaldf.loc[df.index[missing_idx],
    #                 'START_DATE'] = df['diagnosis_int'][missing_idx]-time0
    # Remove unused variable
    del finaldf['diagnosis_int']
    return finaldf


def configure_mafdf(mafdf, keep_samples):
    """Configures a maf dataframe

    Args:
        mafdf: Chunk of maf dataframe
        keep_samples:  Samples to keep in the maf file

    Returns:
        Configured maf dataframe

    """
    keep_mafdf = mafdf[
        mafdf['Tumor_Sample_Barcode'].isin(keep_samples.tolist())
    ]
    if not keep_mafdf.empty:
        fillnas = ['t_depth', 't_ref_count', 't_alt_count',
                   'n_depth', 'n_ref_count', 'n_alt_count']
        for col in fillnas:
            keep_mafdf[col].loc[keep_mafdf[col] == "."] = ""
        keep_mafdf["Validation_Status"] = ''
    return keep_mafdf


def change_days_to_months(days):
    """Changes date fields to months"""
    if math.isnan(days):
        return float('nan')
    else:
        return math.floor(days/30.4)


def change_days_to_years(days):
    """Changes date fields to years"""
    if math.isnan(days):
        return float('nan')
    else:
        return math.floor(days/365.25)


def create_regimens(syn, regimen_infodf, top_x_regimens=5, cohort="NSCLC"):
    """Create regimens to merge into the patient file

    Returns:
        dataframe: Expanded regimen data in patient file format
        dataframe: New regimen mappings for clinical headers
    """
    regimen_synid = regimen_infodf['id'].unique()[0]
    # regimen_synid = "syn22296818"
    regimens_to_exclude = ["Investigational Drug"]
    # HACK: exlude these regimens for test case
    regimens_to_exclude.extend(
        ["Fluorouracil, Irinotecan liposome, Leucovorin Calcium",
         "Fluorouracil, Irinotecan Hydrochloride, Leucovorin Calcium"]
    )
    regimen_ent = syn.get(regimen_synid)
    regimendf = pd.read_csv(regimen_ent.path)
    # Get only NSCLC cohort
    regimendf = regimendf[regimendf['cohort'] == cohort]
    # Use redcap_ca_index == Yes
    regimendf = regimendf[regimendf['redcap_ca_index'] == "Yes"]
    # Exclude regimens
    regimendf = regimendf[~regimendf['regimen_drugs'].isin(regimens_to_exclude)]
    regimendf = regimendf[
        ~regimendf['regimen_drugs'].str.contains("Investigational Drug")
    ]
    # Exclude all regimens with "Other"
    regimendf = regimendf[~regimendf['regimen_drugs'].str.contains("Other")]
    # sort file by regimen_number and drop rest of duplicates
    # (not all duplicates), if duplicated keep the first regimen
    regimendf.sort_values('regimen_number', inplace=True)
    regimendf.drop_duplicates(["record_id", "regimen_drugs"], inplace=True)

    count_of_regimens = regimendf['regimen_drugs'].value_counts()
    # Obtain top X number of regimens
    to_include_regimens = count_of_regimens[:top_x_regimens].index.tolist()

    subset_regimendf = regimendf[
        regimendf['regimen_drugs'].isin(to_include_regimens)
    ]
    regimen_groups = subset_regimendf.groupby("regimen_drugs")
    new_regimen_info = pd.DataFrame()
    # Create regimen clinical headers
    final_regimendf = pd.DataFrame()
    for regimen, df in regimen_groups:
        regimen_drug_info = regimen_infodf.copy()
        regimen_list = regimen.split(",")
        # Create regimen drug abbreviations
        regimen_words = [drug.split() for drug in regimen_list]
        drug_abbr = []
        for drug in regimen_words:
            drug_abbr.append("-".join(word.strip().upper()[:4] for word in drug))
        regimen_abbr = "_".join(drug_abbr)
        # Create correct column mappings for the clinical patient file
        regimen_drug_info['cbio'] = [
            value.format(regimen_abbr=regimen_abbr)
            for value in regimen_infodf['cbio']
        ]
        regimen_drug_info['labels'] = [
            value.format(regimen=regimen)
            for value in regimen_infodf['labels']
        ]
        regimen_drug_info['description'] = [
            value.format(regimen=regimen)
            for value in regimen_infodf['description']
        ]
        regimen_drug_info['priority'] = [
            int(value) for value in regimen_infodf['priority']
        ]
        new_regimen_info = new_regimen_info.append(regimen_drug_info)

        col_map = regimen_drug_info['cbio'].to_dict()
        col_map['record_id'] = "PATIENT_ID"
        regimen_patientdf = df[col_map.keys()].rename(columns=col_map)
        # Merge final regimen dataframe
        if final_regimendf.empty:
            final_regimendf = regimen_patientdf
        else:
            final_regimendf = final_regimendf.merge(
                regimen_patientdf, on="PATIENT_ID", how="outer"
            )
    return {'df': final_regimendf, 'info': new_regimen_info,
            'used': regimen_synid}


class BpcProjectRunner(metaclass=ABCMeta):
    """BPC redcap to cbioportal export"""
    # Sponsorted project name
    _SPONSORED_PROJECT = ''
    # BPC No PHI Data element catalog
    # Version 6 no longer has releaseScope, it has {SP}-sor
    # version 6 doesnt have MSI variables
    _DATA_ELEMENT_SYN_ID = "syn21431364"
    # Redcap codes to cbioportal mapping synid and form key is in
    _REDCAP_TO_CBIOMAPPING_SYNID = None
    # Mapping from Synapse Table to derived variables
    _DATA_TABLE_IDS = None
    # Storage of not found samples
    _SP_REDCAP_EXPORTS_SYNID = None
    # Run `git rev-parse HEAD` in Genie_processing directory to obtain shadigest
    _GITHUB_REPO = None

    def __init__(self, syn, cbiopath, release, staging=False):
        if not os.path.exists(cbiopath):
            raise ValueError("cbiopath doesn't exist")
        if self._SPONSORED_PROJECT == "":
            raise ValueError("Must configure _SPONSORED_PROJECT")
        self.syn = syn
        self.cbiopath = cbiopath
        self.staging = staging
        # Create case lists and release folder
        sp_data_folder = syn.store(Folder(self._SPONSORED_PROJECT,
                                          parentId="syn21241322"))
        release_folder = syn.store(Folder(release, parent=sp_data_folder))
        # Add cBioPortal files into cBioPortal files folder in release
        self._SP_SYN_ID = syn.store(Folder("cBioPortal_files", parent=release_folder))
        self._CASE_LIST_SYN_ID = syn.store(Folder("case_lists",
                                                  parent=self._SP_SYN_ID))
        self.genie_clinicaldf = self.get_main_genie_clinicaldf()

    def get_main_genie_clinicaldf(self) -> dict:
        """Get main GENIE clinical dataframe and perform retraction along
        the way
        """
        # Get all the samples/patients that should be uploaded to SP projects
        # Hard coded clinical database
        # need to use clinical database because SEQ_YEAR is not released in
        # public data
        genie_clinicaldb = self.syn.tableQuery(
            'select SAMPLE_ID, PATIENT_ID, ONCOTREE_CODE, SEQ_ASSAY_ID, '
            'SAMPLE_TYPE, SEQ_YEAR from syn7517674'
        )
        genie_clinicaldf = genie_clinicaldb.asDataFrame()
        # BPC retraction database
        bpc_retraction_db = self.syn.tableQuery(
            'select SAMPLE_ID from syn25779833 where '
            f'{self._SPONSORED_PROJECT} is true'
        )
        bpc_retractiondf = bpc_retraction_db.asDataFrame()
        keep_clinicaldf = genie_clinicaldf[
            ~genie_clinicaldf['SAMPLE_ID'].isin(bpc_retractiondf['SAMPLE_ID'])
        ]
        return keep_clinicaldf

    def download_metadata_files(self):
        """Downloads all the metadata files"""
        # TODO: need to edit the study file
        all_files = self.syn.getChildren(self._SP_SYN_ID)
        for genie_file in all_files:
            if 'meta' in genie_file['name']:
                self.syn.get(genie_file['id'],
                             downloadLocation=self._SPONSORED_PROJECT,
                             ifcollision="overwrite.local")

    def create_genematrixdf(self, clinicaldf, cna_samples, used_ent=None):
        """
        Create gene matrix dataframe
        """
        data_gene_panel = clinicaldf[["SAMPLE_ID", "SEQ_ASSAY_ID"]]
        data_gene_panel = data_gene_panel.rename(
            columns={"SEQ_ASSAY_ID": "mutations"})
        data_gene_panel = data_gene_panel[data_gene_panel['SAMPLE_ID'] != ""]
        data_gene_panel.drop_duplicates("SAMPLE_ID", inplace=True)
        cna_seqids = data_gene_panel['mutations'][
            data_gene_panel['SAMPLE_ID'].isin(cna_samples)].unique()
        data_gene_panel['cna'] = data_gene_panel['mutations']
        data_gene_panel['cna'][~data_gene_panel['cna'].isin(cna_seqids)] = "NA"
        gene_matrix_filepath = os.path.join(self._SPONSORED_PROJECT,
                                            "data_gene_matrix.txt")
        data_gene_panel.to_csv(gene_matrix_filepath, sep="\t", index=False)
        file_ent = File(gene_matrix_filepath,
                        parent=self._SP_SYN_ID)
        if not self.staging:
            self.syn.store(file_ent, used=used_ent,
                           executed=self._GITHUB_REPO)

    def configure_clinicaldf(self, clinicaldf, redcap_to_cbiomappingdf):
        """
        Create clinical file from sponsored project mapped dataframe

        Args:
            clinicaldf:  Can be patient or clinical dataframe
            redcap_to_cbiomappingdf: Synapse Table with the mapping
                                     between redcap and cbioportal
        """

        if not clinicaldf.columns.isin(redcap_to_cbiomappingdf['code']).all():
            raise ValueError("All column names must be in mapping dataframe")
        mapping = redcap_to_cbiomappingdf['cbio'].to_dict()
        clinicaldf = clinicaldf.rename(columns=mapping)
        clinicaldf = clinicaldf.drop_duplicates()
        if sum(clinicaldf['PATIENT_ID'].isnull()) > 0:
            raise ValueError("Must have no null patient ids")
        # Remove white spaces for PATIENT/SAMPLE ID
        clinicaldf['PATIENT_ID'] = [patient.strip()
                                    for patient in clinicaldf['PATIENT_ID']]

        if clinicaldf.get("SAMPLE_ID") is not None:
            # This line should not be here
            clinicaldf = clinicaldf[~clinicaldf['SAMPLE_ID'].isnull()]
            if sum(clinicaldf['SAMPLE_ID'].isnull()) > 0:
                raise ValueError("Must have no null sample ids")
            clinicaldf['SAMPLE_ID'] = [sample.strip()
                                       for sample in clinicaldf['SAMPLE_ID']]

        clinicaldf['SP'] = self._SPONSORED_PROJECT

        for col in clinicaldf:
            num_missing = sum(clinicaldf[col].isnull())
            if num_missing > 0:
                print(f"Number of missing {col}: {num_missing}")

        return clinicaldf

    def write_clinical_file(self, clinicaldf, redcap_to_cbiomappingdf,
                            filetype):
        """Writes out the clinical file

        params:
            clinicaldf: Can be patient or sample clinical dataframe
            redcap_to_cbiomappingdf: mapping dataframe between redcap
                                     and cbioportal
            filetype: "patient" or "sample"
        """
        if filetype not in ["patient", "sample", "supp_survival",
                            "supp_survival_treatment"]:
            raise ValueError(
                "sample type must be patient, sample, supp_survival or "
                "supp_survival_treatment"
            )
        # Must have this for the dict mappings after
        redcap_to_cbiomappingdf.index = redcap_to_cbiomappingdf['cbio']
        label_map = redcap_to_cbiomappingdf['labels'].to_dict()
        description_map = redcap_to_cbiomappingdf['description'].to_dict()
        coltype_map = redcap_to_cbiomappingdf['colType'].to_dict()
        # Priority column will determine which columns are shown on cBioPortal
        # Must columns should be shown on cBioPortal will have a 1
        # but survival columns should all be 0
        priority_map = redcap_to_cbiomappingdf['priority'].to_dict()

        labels = [str(label_map[col]) for col in clinicaldf]
        descriptions = [str(description_map[col]) for col in clinicaldf]
        coltype = [str(coltype_map[col]) for col in clinicaldf]
        priority = [str(int(priority_map[col])) for col in clinicaldf]

        clin_path = os.path.join(self._SPONSORED_PROJECT,
                                 f"data_clinical_{filetype}.txt")

        with open(clin_path, "w+") as clin_file:
            clin_file.write("#{}\n".format("\t".join(labels)))
            clin_file.write("#{}\n".format("\t".join(descriptions)))
            clin_file.write("#{}\n".format("\t".join(coltype)))
            # attributes in the supp file are PATIENT, so must
            # specify that.
            if filetype.startswith("supp_survival"):
                clin_file.write(
                    "#{}\n".format("\t".join(['PATIENT']*len(labels)))
                )
            clin_file.write("#{}\n".format("\t".join(priority)))

            clin_file.write(process_functions.removeStringFloat(
                clinicaldf.to_csv(index=False, sep="\t"))
            )
        return clin_path

    def make_timeline_treatmentdf(self, infodf, sample_type):
        """Make timeline treatment dataframe"""
        subset_infodf = infodf[infodf['sampleType'] == sample_type]
        # Exclude heme onc columns
        subset_infodf = subset_infodf[
            ~subset_infodf['data_type'].isin(['portal_value', 'heme'])
        ]
        synid = subset_infodf['id'].unique()[0]
        ent = self.syn.get(synid)
        used_entity = f'{synid}.{ent.versionNumber}'
        timelinedf = pd.read_csv(ent.path)
        # Only take lung cohort
        timelinedf = timelinedf[
            timelinedf['cohort'] == self._SPONSORED_PROJECT
        ]
        # Only take samples where redcap_ca_index is Yes
        timelinedf = timelinedf[
            timelinedf['redcap_ca_index'] == "Yes"
        ]
        # Flatten multiple columns values into multiple rows
        multiple_cols_idx = subset_infodf['code'].str.contains("[*]")
        final_timelinedf = pd.DataFrame()
        for _, row in subset_infodf[multiple_cols_idx].iterrows():
            code = row['code']
            # Make sure to use regex to get colname + integer
            new_code = code.replace("*", "[\d]")
            cols = timelinedf.columns.str.contains(new_code)
            wanted_cols = timelinedf.columns[cols].tolist()
            wanted_cols.extend(["record_id", "regimen_drugs",
                                'regimen_number'])
            # melt function creates multiple rows from multiple columns
            melted_df = pd.melt(
                timelinedf[wanted_cols],
                id_vars=["record_id", "regimen_drugs", 'regimen_number'],
                value_name=row['cbio']
            )
            del melted_df['variable']
            if final_timelinedf.empty:
                final_timelinedf = final_timelinedf.append(melted_df)
            else:
                final_timelinedf[row['cbio']] = melted_df[row['cbio']]
        final_timelinedf['EVENT_TYPE'] = "Treatment"
        final_timelinedf['TREATMENT_TYPE'] = "Medical Type"
        # Remove all START_DATE is NULL
        final_timelinedf = final_timelinedf[
            ~final_timelinedf['START_DATE'].isnull()
        ]
        non_multi_cols = subset_infodf[~multiple_cols_idx]['code'].tolist()
        non_multi_cols.append("record_id")

        # Merge final timeline
        final_timelinedf = final_timelinedf.merge(
            timelinedf[non_multi_cols],
            on=["record_id", "regimen_drugs", 'regimen_number']
        )
        # Make sure AGENT doesn't have parenthesis
        agents = []
        for index, agent in enumerate(final_timelinedf['AGENT']):
            if "(" in agent:
                agents.append(agent.split("(")[0].strip())
            else:
                # print(final_timelinedf['record_id'][index])
                agents.append(agent.split(",")[0].strip())

        final_timelinedf['AGENT'] = agents
        # Map timeline treatment columns
        mapping = subset_infodf['cbio'].to_dict()
        # Must add in PATIENT_ID
        mapping['record_id'] = 'PATIENT_ID'
        final_timelinedf = final_timelinedf.rename(columns=mapping)
        # timeline file must be in this order
        cols_to_order = ['PATIENT_ID', 'START_DATE', 'STOP_DATE',
                         'EVENT_TYPE', 'TREATMENT_TYPE', 'AGENT']
        cols_to_order.extend(
            final_timelinedf.columns.drop(cols_to_order).tolist()
        )
        final_timelinedf = self.filter_df(final_timelinedf)

        return {'df': final_timelinedf[cols_to_order].drop_duplicates(),
                'used': [used_entity]}

    def filter_df(self, df):
        """Make sure samples and patients exist in the main genie
        clinical samples with the exception of the removal of
        retraction database samples
        """
        if df.get("SAMPLE_ID") is not None:
            to_keep_samples_idx = df['SAMPLE_ID'].isin(
                self.genie_clinicaldf['SAMPLE_ID']
            )
            df = df[to_keep_samples_idx]
        elif df.get("PATIENT_ID") is not None:
            to_keep_patient_idx = df['PATIENT_ID'].isin(
                self.genie_clinicaldf['PATIENT_ID']
            )
            df = df[to_keep_patient_idx]
        return df

    def write_and_storedf(self, df, filepath, used_entities=[]):
        """Write and store dataframe

        Args:
            df: Dataframe to store
            filepath: Path to store dataframe
            used_entities: Used entities
        """
        df_text = process_functions.removePandasDfFloat(df)
        with open(filepath, 'w') as file_f:
            file_f.write(df_text)

        if not self.staging:
            ent = File(filepath, parent=self._SP_SYN_ID)
            self.syn.store(ent, executed=self._GITHUB_REPO, used=used_entities)

    def create_fixed_timeline_files(self, timeline_infodf, timeline_type):
        """Create timeline files straight from derived variables

        Args:
            timeline_infodf: Timeline column mapping dataframe
            timeline_type: Type of timeline

        Returns:
            dict: mapped dataframe,
                  list of used entities
        """
        # Remove portal values
        subset_infodf = timeline_infodf[
            timeline_infodf['sampleType'] == timeline_type
        ]
        # Obtain portal value (EVENT_TYPE)
        portal_value_idx = subset_infodf['data_type'] == 'portal_value'
        portal_value = subset_infodf['code'][portal_value_idx].values[0]
        subset_infodf = subset_infodf[
            subset_infodf['data_type'] != 'portal_value'
        ]
        timeline_data = get_file_data(self.syn, subset_infodf,
                                      timeline_type,
                                      cohort=self._SPONSORED_PROJECT)
        timelinedf = timeline_data['df']

        used_entities = timeline_data['used']

        timelinedf['EVENT_TYPE'] = portal_value

        mapping = subset_infodf['cbio'].to_dict()
        # Must add in PATIENT_ID
        mapping['record_id'] = 'PATIENT_ID'
        timelinedf = timelinedf.rename(columns=mapping)
        timelinedf['STOP_DATE'] = ''
        # timeline file must be in this order
        cols_to_order = ['PATIENT_ID', 'START_DATE', 'STOP_DATE',
                         'EVENT_TYPE']
        cols_to_order.extend(
            timelinedf.columns.drop(cols_to_order).tolist()
        )
        timelinedf = self.filter_df(timelinedf)
        return {'df': timelinedf[cols_to_order].drop_duplicates(),
                'used': used_entities}

    def create_maf(self, keep_samples):
        """Create maf file from release maf

        Args:
            keep_samples: List of samples to keep
        """
        mafpath = os.path.join(self._SPONSORED_PROJECT,
                               "data_mutations_extended.txt")
        # 8.0 maf ent
        maf_synid = "syn22228700"
        maf_ent = self.syn.get(maf_synid)
        maf_chunks = pd.read_table(maf_ent.path, chunksize=50000)
        index = 0
        for maf_chunk in maf_chunks:
            mafdf = configure_mafdf(maf_chunk, keep_samples)
            # Skip to next chunk if empty
            if mafdf.empty:
                continue
            # If maf file has not been created
            if index == 0:
                maf_text = process_functions.removePandasDfFloat(
                    mafdf
                )
                with open(mafpath, 'w') as maf_f:
                    maf_f.write(maf_text)
            else:
                maf_text = mafdf.to_csv(sep="\t", header=None,
                                        index=False)
                maf_text = process_functions.removeStringFloat(
                    maf_text
                )
                with open(mafpath, 'a') as maf_f:
                    maf_f.write(maf_text)
            index += 1
        if not self.staging:
            file_ent = File(mafpath, parent=self._SP_SYN_ID)
            self.syn.store(file_ent, used=[maf_synid],
                           executed=self._GITHUB_REPO)

    def create_cna(self, keep_samples):
        """Create CNA file

        Args:
            keep_samples: List of samples to keep
        """
        cna_synid = "syn22228693"
        cna_path = os.path.join(self._SPONSORED_PROJECT, "data_CNA.txt")
        cna_ent = self.syn.get(cna_synid)
        cnadf = pd.read_table(cna_ent.path)
        keep_cols = ['Hugo_Symbol']
        keep_cols.extend(
            cnadf.columns[cnadf.columns.isin(keep_samples)].tolist()
        )
        cnadf = cnadf[keep_cols]
        cna_text = process_functions.removePandasDfFloat(cnadf)
        # Must do this replace twice because \t\t\t ->
        # \tNA\t\t -> \tNA\tNA\t
        cna_text = cna_text.replace(
            "\t\t", "\tNA\t").replace(
            "\t\t", "\tNA\t").replace(
            '\t\n', "\tNA\n")

        with open(cna_path, "w") as cna_file:
            cna_file.write(cna_text)

        if not self.staging:
            file_ent = File(cna_path, parent=self._SP_SYN_ID)
            self.syn.store(file_ent, used=[cna_synid],
                           executed=self._GITHUB_REPO)
        return cnadf.columns.tolist()

    def create_fusion(self, keep_samples):
        """Create fusion file

        Args:
            keep_samples: List of samples to keep
        """
        fusion_synid = "syn22228696"
        fusion_ent = self.syn.get(fusion_synid)
        fusiondf = pd.read_table(fusion_ent.path)
        fusiondf = fusiondf[
            fusiondf['Tumor_Sample_Barcode'].isin(keep_samples)
        ]
        # cBioPortal validation fails when Hugo Symbol is null
        fusiondf = fusiondf[~fusiondf['Hugo_Symbol'].isnull()]
        fusion_path = os.path.join(self._SPONSORED_PROJECT,
                                   "data_fusions.txt")
        self.write_and_storedf(fusiondf, fusion_path,
                               used_entities=[fusion_synid])

    def create_seg(self, keep_samples):
        """Create seg file

        Args:
            keep_samples: List of samples to keep
        """
        seg_synid = "syn22228729"
        seg_ent = self.syn.get(seg_synid)
        segdf = pd.read_table(seg_ent.path)
        segdf = segdf[segdf['ID'].isin(keep_samples)]
        seg_path = "{}/genie_{}_data_cna_hg19.seg".format(
            self._SPONSORED_PROJECT, self._SPONSORED_PROJECT.lower()
        )
        self.write_and_storedf(segdf, seg_path,
                               used_entities=[seg_synid])

    def create_gene_panels(self, keep_seq_assay_ids):
        """Create gene panels"""
        genomic_info_synid = "syn22228730"
        genomic_info_ent = self.syn.get(genomic_info_synid)
        genomic_infodf = pd.read_table(genomic_info_ent.path)
        # Filter by SEQ_ASSAY_ID and only exonic regions
        genomic_infodf = genomic_infodf[
            (genomic_infodf['SEQ_ASSAY_ID'].isin(keep_seq_assay_ids)) &
            (genomic_infodf['Feature_Type'] == "exon") &
            (~genomic_infodf['Hugo_Symbol'].isnull()) &
            (genomic_infodf['includeInPanel'])
        ]

        seq_assay_groups = genomic_infodf.groupby('SEQ_ASSAY_ID')
        for seq_assay_id, seqdf in seq_assay_groups:
            unique_genes = seqdf.Hugo_Symbol.unique()
            gene_panel_text = (
                "stable_id: {seq_assay_id}\n"
                "description: {seq_assay_id}, "
                "Number of Genes - {num_genes}\n"
                "gene_list:\t{genelist}".format(
                    seq_assay_id=seq_assay_id,
                    num_genes=len(unique_genes),
                    genelist="\t".join(unique_genes)
                )
            )
            gene_panel_name = f"data_gene_panel_{seq_assay_id}.txt"
            gene_panel_path = os.path.join(self._SPONSORED_PROJECT,
                                        gene_panel_name)
            with open(gene_panel_path, "w+") as f:
                f.write(gene_panel_text)
            if not self.staging:
                fileEnt = File(gene_panel_path, parent=self._SP_SYN_ID)
                self.syn.store(fileEnt, used=[genomic_info_synid],
                               executed=self._GITHUB_REPO)

    def run(self):
        """Runs the redcap export to export all files"""
        # Create folder to house release files
        if not os.path.exists(self._SPONSORED_PROJECT):
            os.mkdir(self._SPONSORED_PROJECT)
        else:
            filelists = os.listdir(self._SPONSORED_PROJECT)
            for each_file in filelists:
                if each_file != "case_lists":
                    os.remove(os.path.join(self._SPONSORED_PROJECT,
                                           each_file))
        # Obtain mappings
        # Create full mapping table to get the values of the data model
        # Must use this to determine release scope!
        data_elements = self.syn.tableQuery(
            f"select distinct variable from {self._DATA_ELEMENT_SYN_ID} "
            f'where "{self._SPONSORED_PROJECT}-sor" in '
            "('project', 'consortium')"
        )
        data_elementsdf = data_elements.asDataFrame()
        data_elements_str = "','".join(data_elementsdf['variable'])
        redcap_to_cbiomapping = self.syn.tableQuery(
            f"SELECT * FROM {self._REDCAP_TO_CBIOMAPPING_SYNID} where "
            f"(code in ('{data_elements_str}') or "
            "cbio = 'EVENT_TYPE' or "
            "sampleType in ('TIMELINE-TREATMENT', 'TIMELINE-TREATMENT-RT') and "
            f"data_type <> 'heme') and {self._SPONSORED_PROJECT} is true"
        )
        redcap_to_cbiomappingdf = redcap_to_cbiomapping.asDataFrame()
        data_tables = self.syn.tableQuery(
            f"SELECT id, dataset FROM {self._DATA_TABLE_IDS} "
        )
        data_tablesdf = data_tables.asDataFrame()

        patient_sample_idx = redcap_to_cbiomappingdf['sampleType'].isin(
            ['PATIENT', 'SAMPLE', 'SURVIVAL']
        )
        # Patient and Sample mapping values
        infodf = redcap_to_cbiomappingdf[patient_sample_idx].merge(
            data_tablesdf, on="dataset", how='left'
        )
        infodf.index = infodf['code']

        # Regimen mapping values
        regimen_idx = redcap_to_cbiomappingdf['sampleType'].isin(['REGIMEN'])
        regimen_infodf = redcap_to_cbiomappingdf[regimen_idx].merge(
            data_tablesdf, on="dataset", how='left'
        )
        regimen_infodf.index = regimen_infodf['code']
        # Create timeline column mapping, merges _REDCAP_TO_CBIOMAPPING_SYNID
        # with _DATA_TABLE_IDS
        timeline_infodf = redcap_to_cbiomappingdf[
            ~patient_sample_idx & ~regimen_idx
        ].merge(data_tablesdf, on="dataset", how='left')
        # Add in rt_rt_int for TIMELINE-TREATMENT-RT STOP_DATE
        timeline_infodf = timeline_infodf.append(
            pd.DataFrame({"code": 'rt_rt_int',
                          'sampleType': 'TIMELINE-TREATMENT-RT',
                          'dataset': 'Cancer-Directed Radiation Therapy dataset',
                          'cbio': 'TEMP'},
                         index=['rt_rt_int'])
        )
        # Must do this, because index gets reset after appending
        timeline_infodf.index = timeline_infodf['code']
        # TODO: Must add sample retraction here, also check against main
        # GENIE samples for timeline files...
        print("TREATMENT")
        # Create timeline treatment
        treatment_data = self.make_timeline_treatmentdf(
            timeline_infodf, "TIMELINE-TREATMENT"
        )
        print("TREATMENT-RAD")
        # TODO: Add rt_rt_int
        treatment_rad_data = self.create_fixed_timeline_files(
            timeline_infodf, "TIMELINE-TREATMENT-RT"
        )
        rad_df = treatment_rad_data['df']
        rad_df['STOP_DATE'] = rad_df['START_DATE'] + rad_df['TEMP']
        rad_df = rad_df[rad_df['REDCAP_CA_INDEX'] == "Yes"]
        del rad_df['REDCAP_CA_INDEX']
        del rad_df['TEMP']
        treatment_data['df'] = treatment_data['df'].append(
            rad_df
        )
        treatment_path = os.path.join(self._SPONSORED_PROJECT,
                                      "data_timeline_treatment.txt")
        self.write_and_storedf(treatment_data['df'], treatment_path,
                               used_entities=treatment_data['used'])

        # Create static timeline files
        # Cancer dx
        print("DX")
        cancerdx_data = self.create_fixed_timeline_files(
            timeline_infodf, "TIMELINE-DX"
        )
        cancerdx_data['df'] = fill_cancer_dx_start_date(cancerdx_data['df'])
        cancerdx_path = os.path.join(self._SPONSORED_PROJECT,
                                     "data_timeline_cancer_diagnosis.txt")
        # There are patients >89 that don't have START_DATE.
        # These must be removed
        cancerdx_data['df'] = cancerdx_data['df'][
            ~cancerdx_data['df']['START_DATE'].isnull()
        ]
        self.write_and_storedf(cancerdx_data['df'], cancerdx_path,
                               used_entities=cancerdx_data['used'])

        # Pathology Data
        print("PATHOLOGY")
        pathology_data = self.create_fixed_timeline_files(
            timeline_infodf, "TIMELINE-PATHOLOGY"
        )
        pathology_path = os.path.join(self._SPONSORED_PROJECT,
                                      "data_timeline_pathology.txt")
        self.write_and_storedf(pathology_data['df'], pathology_path,
                               used_entities=pathology_data['used'])

        # # Sample acquisition
        # This is important because path_num_spec is needed
        # timeline_infodf = timeline_infodf.append(
        #     pd.DataFrame({"code": 'path_num_spec',
        #                   'sampleType': 'TIMELINE-SAMPLE',
        #                   'dataset': 'Pathology-report level dataset',
        #                   'cbio': 'TEMP'},
        #                  index=['path_num_spec'])
        # )
        print("SAMPLE-ACQUISITION")
        acquisition_data = self.create_fixed_timeline_files(
            timeline_infodf, "TIMELINE-SAMPLE"
        )
        acquisition_path = os.path.join(
            self._SPONSORED_PROJECT, "data_timeline_sample_acquisition.txt"
        )
        # acquisitiondf = acquisition_data['df']
        # keep_idx = acquisitiondf['TEMP'] == acquisitiondf['PATH_PROC_NUMBER']
        # acquisitiondf.drop(columns="TEMP", inplace=True)
        # acquisition_data['df'] = acquisitiondf[keep_idx]

        # TODO: Can add getting of samples with NULL start dates in
        # self.create_fixed_timeline_files
        null_dates_idx = acquisition_data['df']['START_DATE'].isnull()
        if null_dates_idx.any():
            print("timeline sample with null START_DATE: {}".format(
                ", ".join(acquisition_data['df']['SAMPLE_ID'][null_dates_idx])
            ))
        self.write_and_storedf(acquisition_data['df'][~null_dates_idx],
                               acquisition_path,
                               used_entities=acquisition_data['used'])

        # Medonc
        print("MEDONC")
        medonc_data = self.create_fixed_timeline_files(
            timeline_infodf, "TIMELINE-MEDONC"
        )
        medonc_path = os.path.join(
            self._SPONSORED_PROJECT, "data_timeline_medonc.txt"
        )
        self.write_and_storedf(medonc_data['df'], medonc_path,
                               used_entities=medonc_data['used'])

        # Imaging
        print("IMAGING")

        imaging_data = self.create_fixed_timeline_files(
            timeline_infodf, "TIMELINE-IMAGING"
        )
        imaging_path = os.path.join(
            self._SPONSORED_PROJECT, "data_timeline_imaging.txt"
        )
        self.write_and_storedf(imaging_data['df'], imaging_path,
                               used_entities=imaging_data['used'])

        # Sequencing
        print("SEQUENCE")
        sequence_data = self.create_fixed_timeline_files(
            timeline_infodf, "TIMELINE-SEQUENCE"
        )
        sequence_path = os.path.join(
            self._SPONSORED_PROJECT, "data_timeline_sequencing.txt"
        )
        self.write_and_storedf(sequence_data['df'], sequence_path,
                               used_entities=sequence_data['used'])

        if self._SPONSORED_PROJECT in ["CRC", "BrCa"]:
            # Lab test
            print("LABTEST")
            lab_data = self.create_fixed_timeline_files(
                timeline_infodf, "TIMELINE-LAB"
            )
            lab_path = os.path.join(
                self._SPONSORED_PROJECT, "data_timeline_labtest.txt"
            )
            self.write_and_storedf(lab_data['df'], lab_path,
                                   used_entities=sequence_data['used'])

        # supplemental clinical file
        print("SURVIVAL")
        # This is important because tt_first_index_ca is needed
        # For filtering
        infodf = infodf.append(
            pd.DataFrame({"code": 'tt_first_index_ca',
                          'sampleType': 'SURVIVAL',
                          'dataset': 'Cancer-level index dataset',
                          'cbio': 'CANCER_INDEX'},
                         index=['tt_first_index_ca'])
        )
        # Must do this because index gets reset
        infodf.index = infodf['code']
        survival_infodf = infodf[infodf['sampleType'] == "SURVIVAL"]
        survival_data = get_file_data(self.syn, survival_infodf, "SURVIVAL",
                                      cohort=self._SPONSORED_PROJECT)
        survivaldf = survival_data['df']
        final_survivaldf = self.configure_clinicaldf(survivaldf,
                                                     survival_infodf)
        # Only take rows where cancer index is null
        final_survivaldf = final_survivaldf[
            final_survivaldf['CANCER_INDEX'].isnull()
        ]
        # Remove cancer index column
        del final_survivaldf['CANCER_INDEX']
        print("PATIENT")
        # Patient and sample files
        patient_infodf = infodf[infodf['sampleType'] == "PATIENT"]
        # Must get redcap_ca_index to grab only the index cancers
        patient_infodf = patient_infodf.append(
            {"code": 'redcap_ca_index', 'sampleType': 'PATIENT',
             'dataset': 'Cancer-level dataset'},
            ignore_index=True
        )
        patient_infodf.index = patient_infodf['code']
        patient_data = get_file_data(self.syn, patient_infodf, "PATIENT",
                                     cohort=self._SPONSORED_PROJECT)

        patientdf = patient_data['df']
        # Subset by redcap_ca_index == Yes
        patientdf = patientdf[patientdf['redcap_ca_index'] == "Yes"]
        # remove columns after done with it
        patientdf.drop(columns="redcap_ca_index", inplace=True)
        final_patientdf = self.configure_clinicaldf(patientdf, infodf)

        print('SAMPLE')
        sample_infodf = infodf[infodf['sampleType'] == "SAMPLE"]
        sample_data = get_file_data(self.syn, sample_infodf, "SAMPLE",
                                    cohort=self._SPONSORED_PROJECT)
        sampledf = sample_data['df']
        del sampledf['path_proc_number']
        # SAMPLE FILE
        final_sampledf = self.configure_clinicaldf(sampledf, infodf)

        # Only patients and samples that exist in the
        # sponsored project uploads are going to be pulled into the SP project

        to_keep_patient_idx = final_patientdf['PATIENT_ID'].isin(
            self.genie_clinicaldf['PATIENT_ID']
        )
        subset_patientdf = final_patientdf[final_patientdf['PATIENT_ID'].isin(
            self.genie_clinicaldf['PATIENT_ID']
        )]

        # Fix patient duplicated values due to cancer index DOB
        # Take the larger DX_LASTALIVE_INT_MOS value for all records
        # TODO: check if its fine to drop this column
        # subset_patientdf.sort_values("DX_LASTALIVE_INT_MOS", inplace=True,
        #                              ascending=False)
        subset_patientdf.drop_duplicates("PATIENT_ID", inplace=True)

        duplicated = subset_patientdf.PATIENT_ID.duplicated()
        if duplicated.any():
            print("DUPLICATED PATIENT_IDs: {}".format(
                ",".join(subset_patientdf['PATIENT_ID'][duplicated])
            ))

        del subset_patientdf['SP']
        cols_to_order = ['PATIENT_ID']
        cols_to_order.extend(
            subset_patientdf.columns.drop(cols_to_order).tolist()
        )
        # HACK: Temporary remapping of specific values in a column
        laterality_mapping = {
            0: "Not a paired site",
            1: "Right: origin of primary",
            2: "Left: origin of primary",
            3: "Only one side involved, right or left origin unspecified",
            4: "Bilateral involvement at time of diagnosis, lateral origin "
               "unknown for a single primary; or both ovaries involved "
               "simultaneously, single histology; bilateral retinoblastomas; "
               "bilateral Wilms' tumors",
            5: "Paired site: midline tumor",
            9: "Paired site, but no information concerning laterality"
        }
        if subset_patientdf.get("NAACCR_LATERALITY_CD") is not None:
            remapped_values = subset_patientdf['NAACCR_LATERALITY_CD'].map(
                laterality_mapping
            )
            subset_patientdf['NAACCR_LATERALITY_CD'] = remapped_values

        # Write patient file out
        patient_path = self.write_clinical_file(
            subset_patientdf[cols_to_order], infodf, "patient"
        )
        # Create regimens data for patient file
        regimens_data = create_regimens(self.syn, regimen_infodf,
                                        top_x_regimens=20,
                                        cohort=self._SPONSORED_PROJECT)

        survival_info = infodf.append(regimens_data['info'])

        # Create survival data
        subset_survivaldf = final_survivaldf[
            final_survivaldf['PATIENT_ID'].isin(
                self.genie_clinicaldf['PATIENT_ID']
            )
        ]
        del subset_survivaldf['SP']
        # subset_survivaldf = subset_survivaldf.merge(
        #     regimens_data['df'], on="PATIENT_ID", how="left"
        # )
        # Must change the values in OS and PFS columns from
        # integer to (status:label)
        # TODO: put remapping of os values in helper function
        os_pfs_cols = [col for col in subset_survivaldf.columns
                       if col.startswith(('OS', 'PFS')) and
                       col.endswith("STATUS")]
        remap_os_values = {col: {0: "0:LIVING", 1: "1:DECEASED"}
                           for col in os_pfs_cols}
        subset_survivaldf.replace(remap_os_values, inplace=True)
        cols_to_order = ['PATIENT_ID']
        cols_to_order.extend(
            subset_survivaldf.columns.drop(cols_to_order).tolist()
        )
        # Order is maintained in the derived variables file so just drop
        # Duplicates
        # subset_survivaldf.drop_duplicates("PATIENT_ID", inplace=True)
        survival_path = self.write_clinical_file(
            subset_survivaldf[cols_to_order], survival_info, "supp_survival"
        )
        # survival treatment
        survival_treatmentdf = regimens_data['df']
        os_pfs_cols = [col for col in survival_treatmentdf.columns
                       if col.startswith(('OS', 'PFS')) and
                       col.endswith("STATUS")]
        remap_os_values = {col: {0: "0:LIVING", 1: "1:DECEASED"}
                           for col in os_pfs_cols}
        survival_treatmentdf.replace(remap_os_values, inplace=True)
        cols_to_order = ['PATIENT_ID']
        cols_to_order.extend(
            survival_treatmentdf.columns.drop(cols_to_order).tolist()
        )

        # Order is maintained in the derived variables file so just drop
        # Duplicates
        # survival_treatmentdf.drop_duplicates("PATIENT_ID", inplace=True)
        surv_treatment_path = self.write_clinical_file(
            survival_treatmentdf[cols_to_order], survival_info,
            "supp_survival_treatment"
        )
        # Fill in ONCOTREE_CODE
        final_sampledf['ONCOTREE_CODE'] = [
            self.genie_clinicaldf['ONCOTREE_CODE'][
                self.genie_clinicaldf['SAMPLE_ID'] == sample].values[0]
            if sum(self.genie_clinicaldf['SAMPLE_ID'] == sample) > 0
            else float('nan') for sample in final_sampledf['SAMPLE_ID']]
        # Fill in SEQ_ASSAY_ID
        final_sampledf['SEQ_ASSAY_ID'] = [
            self.genie_clinicaldf['SEQ_ASSAY_ID'][
                self.genie_clinicaldf['SAMPLE_ID'] == sample].values[0]
            if sum(self.genie_clinicaldf['SAMPLE_ID'] == sample) > 0
            else float('nan') for sample in final_sampledf['SAMPLE_ID']]

        subset_sampledf = final_sampledf[final_sampledf['SAMPLE_ID'].isin(
            self.genie_clinicaldf['SAMPLE_ID']
        )]
        del subset_sampledf['SP']
        # TODO: add YEARS
        # days_to_years_col = ['AGE_AT_SEQ_REPORT_YEARS',
        #                      'CPT_ORDER_INT', 'CPT_REPORT_INT']
        days_to_years_col = ['AGE_AT_SEQ_REPORT_YEARS', 'CPT_ORDER_INT',
                             'CPT_REPORT_INT', 'AGE_AT_SEQUENCING']
        for col in days_to_years_col:
            # not all columns could exist, so check if column exists
            if col in subset_sampledf:
                years = subset_sampledf[col].apply(change_days_to_years)
                subset_sampledf[col] = years

        # Remove SAMPLE_TYPE and CPT_SEQ_DATE because the values are incorrect
        del subset_sampledf['SAMPLE_TYPE']
        del subset_sampledf['CPT_SEQ_DATE']
        # Obtain this information from the main GENIE cohort
        subset_sampledf = subset_sampledf.merge(
            self.genie_clinicaldf[['SAMPLE_ID', 'SAMPLE_TYPE', 'SEQ_YEAR']],
            on="SAMPLE_ID", how="left"
        )
        subset_sampledf.rename(columns={"SEQ_YEAR": "CPT_SEQ_DATE"},
                               inplace=True)
        # Remove duplicated samples due to PDL1
        # Keep only one sample in this priority
        # PDL1_POSITIVE_ANY: Yes
        # PDL1_POSITIVE_ANY: No
        # PDL1_POSITIVE_ANY: <blank>
        subset_sampledf.sort_values(
            "PDL1_POSITIVE_ANY", ascending=False, inplace=True
        )
        subset_sampledf.drop_duplicates("SAMPLE_ID", inplace=True)
        # duplicated = subset_sampledf.SAMPLE_ID.duplicated()
        # if duplicated.any():
        #     # TODO: Add in duplicated ids
        #     print("DUPLICATED SAMPLE_IDs")
            # There are duplicated samples
            # subset_sampledf = subset_sampledf[~duplicated]
        sample_path = self.write_clinical_file(subset_sampledf, infodf,
                                               "sample")

        # Remove oncotree code here, because no longer need it
        merged_clinicaldf = subset_sampledf.merge(subset_patientdf,
                                                  on="PATIENT_ID",
                                                  how="outer")
        missing_sample_idx = merged_clinicaldf['SAMPLE_ID'].isnull()
        # Make sure there are no missing sample ids
        if sum(missing_sample_idx) > 0:
            print("MISSING SAMPLE_ID for: {}".format(
                ",".join(merged_clinicaldf['PATIENT_ID'][missing_sample_idx])
            ))
            merged_clinicaldf = merged_clinicaldf[~missing_sample_idx]

        # upload samples that are not part of the main GENIE cohort
        if merged_clinicaldf.get("SAMPLE_ID") is not None:
            print("Samples not in GENIE clinical databases (SP and normal)")
            not_found_samples = merged_clinicaldf[
                'SAMPLE_ID'][~merged_clinicaldf['SAMPLE_ID'].isin(
                    self.genie_clinicaldf['SAMPLE_ID']
                )
            ]
            if not not_found_samples.empty:
                print(not_found_samples[~not_found_samples.isnull()])
                not_found_samples.to_csv("notfoundsamples.csv")
                if not self.staging:
                    self.syn.store(synapseclient.File(
                        "notfoundsamples.csv",
                        parent=self._SP_REDCAP_EXPORTS_SYNID
                    ))

        # Hard coded most up to date oncotree version
        oncotreelink = self.syn.get("syn13890902").externalURL
        # Use the old oncotree link for now
        oncotreelink = 'http://oncotree.mskcc.org/api/tumorTypes/tree?version=oncotree_2018_06_01'
        oncotree_dict = genie.process_functions.get_oncotree_code_mappings(oncotreelink)
        # Map cancer type and cancer type detailed
        # This is to create case list files
        merged_clinicaldf['CANCER_TYPE'] = [
            oncotree_dict[code.upper()].get("CANCER_TYPE", float('nan'))
            for code in merged_clinicaldf['ONCOTREE_CODE']
        ]
        merged_clinicaldf['CANCER_TYPE_DETAILED'] = [
            oncotree_dict[code.upper()].get("CANCER_TYPE_DETAILED",
                                            float('nan'))
            for code in merged_clinicaldf['ONCOTREE_CODE']
        ]
        merged_clinicaldf['ONCOTREE_PRIMARY_NODE'] = [
            oncotree_dict[code.upper()].get("ONCOTREE_PRIMARY_NODE",
                                            float('nan'))
            for code in merged_clinicaldf['ONCOTREE_CODE']
        ]
        merged_clinicaldf['ONCOTREE_SECONDARY_NODE'] = [
            oncotree_dict[code.upper()].get("ONCOTREE_SECONDARY_NODE",
                                            float('nan'))
            for code in merged_clinicaldf['ONCOTREE_CODE']
        ]
        # Remove duplicated sample ids (there shouldn't be any)
        merged_clinicaldf = merged_clinicaldf.drop_duplicates("SAMPLE_ID")
        merged_clinicaldf.to_csv(os.path.join(self._SPONSORED_PROJECT,
                                              "data_clinical.txt"),
                                 index=False, sep="\t")
        if not self.staging:
            patient_fileent = File(patient_path, parent=self._SP_SYN_ID)
            patient_ent = self.syn.store(patient_fileent,
                                         used=patient_data['used'],
                                         executed=self._GITHUB_REPO)

            sample_fileent = File(sample_path, parent=self._SP_SYN_ID)
            sample_ent = self.syn.store(sample_fileent,
                                        used=sample_data['used'],
                                        executed=self._GITHUB_REPO)

            survival_fileent = File(survival_path, parent=self._SP_SYN_ID)
            used = survival_data['used']
            used.append(regimens_data['used'])
            survival_ent = self.syn.store(survival_fileent, used=used,
                                          executed=self._GITHUB_REPO)

            survival_treatment_fileent = File(
                surv_treatment_path, parent=self._SP_SYN_ID
            )
            used = survival_data['used']
            used.append(regimens_data['used'])
            survival_treatment_fileent = self.syn.store(
                survival_treatment_fileent, used=used,
                executed=self._GITHUB_REPO
            )
        self.create_maf(subset_sampledf['SAMPLE_ID'])

        cna_samples = self.create_cna(subset_sampledf['SAMPLE_ID'])

        self.create_genematrixdf(subset_sampledf, cna_samples)

        self.create_fusion(subset_sampledf['SAMPLE_ID'])

        self.create_seg(subset_sampledf['SAMPLE_ID'])

        # Create case lists
        case_list_path = os.path.join(self._SPONSORED_PROJECT,
                                      'case_lists')

        if not os.path.exists(case_list_path):
            os.mkdir(case_list_path)
        else:
            caselists = os.listdir(case_list_path)
            for caselist in caselists:
                os.remove(os.path.join(case_list_path, caselist))

        # Write out cases sequenced so people can tell
        # which samples were sequenced
        assay_info = self.syn.tableQuery("select * from syn17009222",
                                         includeRowIdAndRowVersion=False,
                                         separator="\t")
        create_case_lists.main(
            os.path.join(self._SPONSORED_PROJECT, "data_clinical.txt"),
            assay_info.filepath,
            case_list_path,
            f"{self._SPONSORED_PROJECT.lower()}_genie_bpc"
        )

        case_list_files = os.listdir(case_list_path)
        for casepath in case_list_files:
            casepath = os.path.join(case_list_path, casepath)
            if not self.staging:
                file_ent = File(casepath, parent=self._CASE_LIST_SYN_ID)
                self.syn.store(file_ent,
                               used=[patient_ent.id, sample_ent.id],
                               executed=self._GITHUB_REPO)

        self.create_gene_panels(subset_sampledf['SEQ_ASSAY_ID'].unique())
        # Make sure to re download all the metadata files again
        self.download_metadata_files()

        cmd = ['python',
               os.path.join(self.cbiopath,
                            "core/src/main/scripts/importer/validateData.py"),
               "-s", self._SPONSORED_PROJECT, "-n"]
        subprocess.call(cmd)
