"""
Sponsored project configuration classes

USAGE:
git clone https://github.com/cBioPortal/cbioportal.git
python runSP.py AKT1 ../cbioportal/ --staging
"""
import os
import random
import string

import pandas as pd

import new_redcap_export_mapping
import sp_redcap_export_mapping


class Akt1(SP_redcap_export_mapping.SponsoredProjectRunner):
    """
    ########################################################################
    # AKT1 PROCESSES
    # - ONE TIMELINE FILE
    # - CLINICAL FILE
    #   OS_MONTHS = death_date_int - mets_disease_date_int
    #   OS_MONTHS_PRIMARY = death_date_int - primary_dx_date_int 
    #   All dates are converted from days to months (days/30.4)
    #   Add headers
    #   REMOVE PATIENTS/SAMPLES THAT DON'T HAVE GENIE SAMPLE IDS
    ########################################################################
    """
    _SPONSORED_PROJECT = "AKT1"
    _DATES = ["death_date_int","follow_up_date_int","primary_dx_date_int","lrr_date_int","mets_disease_date_int","sample_date_int_1",
              "sequence_report_date_int_1","sequence_report_date_int_1_static","sample_date_int_2","sample_date_int_2_static",
              "sequence_report_date_int_2","sequence_report_date_int_2_static","sequence_report_date_int_3_static",
              "OS_MONTHS","OS_MONTHS_PRIMARY"]
    _CASE_LIST_MAF_SAMPLES_TEMPLATE = "cancer_study_identifier: genie_akt1\nstable_id: genie_akt1_sequenced\ncase_list_category: all_cases_with_mutation_data\ncase_list_name: Sequenced Tumors\ncase_list_description: All sequenced samples (%s samples)\ncase_list_ids: %s"
    _CASE_LIST_PATH = os.path.join(_SPONSORED_PROJECT,'case_lists')
    _UNMAPPED_SYN_ID = "syn11066652"
    _MAPPED_SYN_ID = "syn8404878"
    _CASE_LIST_SYN_ID = "syn10145838"
    _SP_SYN_ID = "syn8363325"
    _REDCAP_TO_CBIOMAPPING_SYNID = "syn8220815"
    _SP_REDCAP_EXPORTS_SYNID = "syn8404875" #Storage of not found samples
    _NUM_SAMPLE_COLS = 3

    def addOSMonths(self, sponsoredProject_mapped_df):
        #Must add new date fields to the DATE variable along with add to the mapping table: syn8220815
        sponsoredProject_mapped_df['OS_MONTHS'] = sponsoredProject_mapped_df['death_date_int'] - sponsoredProject_mapped_df['mets_disease_date_int']  
        sponsoredProject_mapped_df['OS_MONTHS_PRIMARY'] = sponsoredProject_mapped_df['death_date_int'] - sponsoredProject_mapped_df['primary_dx_date_int']  
        return(sponsoredProject_mapped_df)

    def createTemporaryGenieId(self, x, tempIdMapping):
        uniqId = x['record_id'] + x['redcap_data_access_group']
        tempIdMap = tempIdMapping['patientId'][tempIdMapping['uniqueId'] == uniqId]
        tempId = 'GENIE-%s-%s' % (x['redcap_data_access_group'],''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(10)))
        if len(tempIdMap) == 0:
            return(tempId)
        else:
            return(tempIdMap.values[0])
        # if sum(tempIdMapping['uniqueId'] == uniqId) == 0:
        #   #syn.store(synapseclient.Table(syn.get("syn10164044"),[[uniqId, tempId, SPONSORED_PROJECT]]))
        #   return(tempId)
        # elif pd.np.isnan(temp['tempPatientId'][tempIdMapping['uniqueId'] == uniqId].values[0]):

        # else:
        #   return(tempIdMapping[tempIdMapping['uniqueId'] == uniqId]['tempPatientId'].values[0])

    def createNullPatients(self, sponsoredProject_mapped_df, tempIdMappingDf):
        print("RENAMING %s NULL PATIENTS" % sum(sponsoredProject_mapped_df['genie_patient_id'].isnull()))
        #Create temp patient Id
        allNullPatients = sponsoredProject_mapped_df[['record_id','redcap_data_access_group','genie_patient_id']][sponsoredProject_mapped_df['genie_patient_id'].isnull()]
        temporaryIds = allNullPatients.apply(lambda x: self.createTemporaryGenieId(x, tempIdMappingDf), axis =1)
        if sponsoredProject_mapped_df['genie_patient_id'].isnull().any():
            sponsoredProject_mapped_df['genie_patient_id'][sponsoredProject_mapped_df['genie_patient_id'].isnull()] = temporaryIds
        assert sum(sponsoredProject_mapped_df['genie_patient_id'].isnull()) ==0, "Make sure there are no null genie patient Ids"

        sponsoredProject_mapped_df['genie_patient_id'] = sponsoredProject_mapped_df.apply(lambda x: self.checkGenieId(x, 'redcap_data_access_group','genie_patient_id'), axis=1)
        sponsoredProject_mapped_df.reset_index(inplace=True,drop=True)
        return(sponsoredProject_mapped_df, temporaryIds)

    def makeTimeLineDf(self, redCapExportDf, therapyRange = 18):
        START_DATE = []
        STOP_DATE = []
        TREATMENT_TYPE = []
        SUBTYPE = []
        AGENT = []
        THERAPY_DRUG_CLINTRIAL = []
        THERAPY_DRUG_AZD5363 = []
        THERAPY_DRUG_OTHER = []
        THERAPY_DRUG_DISCONTINUE = []
        THERAPY_DRUG_REASON = []
        THERAPY_COMBO_YN = []
        THERAPY_COMBO_NUM = []
        #THERAPY NUMBER
        for therapyNumber in range(1,therapyRange):
            therapyCols = [i for i in redCapExportDf if "therapy%d_" % therapyNumber in i]
            START_DATE.extend([i for i in therapyCols if "start_int" in i])
            STOP_DATE.extend([i for i in therapyCols if "end_int" in i])
            AGENT.extend([i for i in therapyCols if len(i.split("_")) == 2])
            THERAPY_DRUG_CLINTRIAL.extend([i for i in therapyCols if "clintrial" in i])
            THERAPY_DRUG_AZD5363.extend([i for i in therapyCols if "azd" in i])
            THERAPY_DRUG_OTHER.extend([i for i in therapyCols if "other" in i])
            THERAPY_DRUG_DISCONTINUE.extend([i for i in therapyCols if "discontinue" in i])
            THERAPY_DRUG_REASON.extend([i for i in therapyCols if "reason" in i])
            THERAPY_COMBO_YN.extend([i for i in therapyCols if "combo_yn" in i] * len([i for i in therapyCols if "start_int" in i]))
            THERAPY_COMBO_NUM.extend([i for i in therapyCols if "combo_num" in i]* len([i for i in therapyCols if "start_int" in i]))
            TREATMENT_TYPE.extend(["Medical Therapy %d" % therapyNumber]* len([i for i in therapyCols if "start_int" in i]))
            SUBTYPE.extend(["Chemo/Target/Immuno etc."] * len([i for i in therapyCols if "start_int" in i]))
        #OVARIAN
        ovarian = [i for i in redCapExportDf if "ovariansup" in i]
        ovarian_len = len([i for i in ovarian if "start_int" in i])
        START_DATE.extend([i for i in ovarian if "start_int" in i])
        STOP_DATE.extend([i for i in ovarian if "end_int" in i])
        TREATMENT_TYPE.extend(["Ovarian Suppression At Primary"] * ovarian_len)
        SUBTYPE.extend(["Ovarian Suppression"] * ovarian_len)
        AGENT.extend(['']*ovarian_len)
        THERAPY_DRUG_CLINTRIAL.extend(['']*ovarian_len)
        THERAPY_DRUG_AZD5363.extend(['']*ovarian_len)
        THERAPY_DRUG_OTHER.extend(['']*ovarian_len)
        THERAPY_DRUG_DISCONTINUE.extend(['']*ovarian_len)
        THERAPY_DRUG_REASON.extend(['']*ovarian_len)
        THERAPY_COMBO_YN.extend(['']*ovarian_len)
        THERAPY_COMBO_NUM.extend(['']*ovarian_len)
        #HORMONE
        hormo = [i for i in redCapExportDf if "hormo" in i]
        hormo_len = len([i for i in hormo if "start_int" in i])
        START_DATE.extend([i for i in hormo if "start_int" in i])
        STOP_DATE.extend([i for i in hormo if "end_int" in i])
        THERAPY_DRUG_CLINTRIAL.extend([i for i in hormo if "clintrial" in i])
        THERAPY_DRUG_AZD5363.extend(['']*hormo_len)
        THERAPY_DRUG_OTHER.extend([i for i in hormo if "other" in i])
        THERAPY_DRUG_DISCONTINUE.extend([i for i in hormo if "discon" in i])
        THERAPY_DRUG_REASON.extend([i for i in hormo if "reason" in i])
        AGENT.extend([i for i in hormo if "reason" not in i and "discon" not in i and "other" not in i and "clintrial" not in i and "start_int" not in i and "end_int" not in i and "therapy" not in i])
        THERAPY_COMBO_YN.extend(['']*hormo_len)
        THERAPY_COMBO_NUM.extend(['']*hormo_len)
        SUBTYPE.extend(["Hormone Therapy"] * hormo_len)
        TREATMENT_TYPE.extend(["Medical Therapy 1"] * hormo_len)
        EVENT_TYPE = ["TREATMENT"]*len(AGENT)

        #METASTATIC DIAGNOSIS
        metaDiagnosis = pd.DataFrame()
        metaDiagnosis['PATIENT_ID'] = redCapExportDf['genie_patient_id']
        #MET DISEASE IS TIMEPOINT 0
        metaDiagnosis['START_DATE'] = 0
        #metaDiagnosis['START_DATE'] = redCapExportDf['mets_disease_date_int']
        metaDiagnosis['EVENT_TYPE'] = 'STATUS'
        metaDiagnosis['STATUS'] = 'Metastatic Diagnosis'
        metaDiagnosis = metaDiagnosis[~metaDiagnosis['START_DATE'].isnull()]

        removeCols = START_DATE+STOP_DATE+AGENT+THERAPY_DRUG_CLINTRIAL+THERAPY_DRUG_AZD5363+THERAPY_DRUG_OTHER+THERAPY_DRUG_DISCONTINUE+THERAPY_DRUG_REASON+THERAPY_COMBO_YN+THERAPY_COMBO_NUM
        lengths = set([
        len(START_DATE),
        len(STOP_DATE),
        len(TREATMENT_TYPE),
        len(SUBTYPE),
        len(AGENT),
        len(THERAPY_DRUG_CLINTRIAL),
        len(THERAPY_DRUG_AZD5363),
        len(THERAPY_DRUG_OTHER),
        len(THERAPY_DRUG_DISCONTINUE),
        len(THERAPY_DRUG_REASON),
        len(THERAPY_COMBO_YN),
        len(THERAPY_COMBO_NUM),
        len(EVENT_TYPE)])
        assert len(lengths) == 1,"Lengths must all be the same"

        total = pd.DataFrame()
        for i in range(len(redCapExportDf)):
            timelineDF = pd.DataFrame()
            timelineDF['PATIENT_ID'] = [redCapExportDf['genie_patient_id'][i]]*len(START_DATE)
            #timelineDF['START_DATE'] = redCapExportDf.ix[i][START_DATE].reset_index(drop=True) - redCapExportDf.ix[i]['primary_dx_date_int']
            #timelineDF['STOP_DATE'] = redCapExportDf.ix[i][STOP_DATE].reset_index(drop=True) - redCapExportDf.ix[i]['primary_dx_date_int']
            #MET DISEASE IS TIMEPOINT 0
            timelineDF['START_DATE'] = redCapExportDf.iloc[i][START_DATE].reset_index(drop=True) - redCapExportDf.iloc[i]['mets_disease_date_int']
            timelineDF['STOP_DATE'] = redCapExportDf.iloc[i][STOP_DATE].reset_index(drop=True) - redCapExportDf.iloc[i]['mets_disease_date_int']
            timelineDF['EVENT_TYPE'] = EVENT_TYPE
            timelineDF['TREATMENT_TYPE'] = TREATMENT_TYPE
            timelineDF['SUBTYPE'] = SUBTYPE
            timelineDF['AGENT'] = redCapExportDf.iloc[i][AGENT].reset_index(drop=True)
            timelineDF['THERAPY_DRUG_CLINTRIAL'] = redCapExportDf.iloc[i][THERAPY_DRUG_CLINTRIAL].reset_index(drop=True)
            timelineDF['THERAPY_DRUG_AZD5363'] = redCapExportDf.iloc[i][THERAPY_DRUG_AZD5363].reset_index(drop=True)
            timelineDF['THERAPY_DRUG_OTHER'] = redCapExportDf.iloc[i][THERAPY_DRUG_OTHER].reset_index(drop=True)
            timelineDF['THERAPY_DRUG_DISCONTINUE'] = redCapExportDf.iloc[i][THERAPY_DRUG_DISCONTINUE].reset_index(drop=True)
            timelineDF['THERAPY_DRUG_REASON'] = redCapExportDf.iloc[i][THERAPY_DRUG_REASON].reset_index(drop=True)
            timelineDF['THERAPY_COMBO_YN'] = redCapExportDf.iloc[i][THERAPY_COMBO_YN].reset_index(drop=True)
            timelineDF['THERAPY_COMBO_NUM'] = redCapExportDf.iloc[i][THERAPY_COMBO_NUM].reset_index(drop=True)
            total = total.append(timelineDF)
        total['STATUS'] = ''
        ordering = total.columns
        total = total.append(metaDiagnosis)
        total = total[ordering]
        return(total,removeCols)

    def getSpecimen(self, getTimelineSpecimen):
        specimen = pd.DataFrame()
        specimen['PATIENT_ID'] = getTimelineSpecimen['PATIENT_ID']
        specimen['START_DATE'] = getTimelineSpecimen.SEQUENCE_REPORT_DATE_INT_STATIC - getTimelineSpecimen.METS_DISEASE_DATE_INT
        specimen['EVENT_TYPE'] = 'SPECIMEN'
        specimen['SAMPLE_ID'] = getTimelineSpecimen['SAMPLE_ID']
        specimen['SAMPLE_NOTES'] = getTimelineSpecimen.SEQUENCE_REPORT_DATE_INT_STATIC
        specimen = specimen[~specimen['START_DATE'].isnull()]
        return(specimen)


class Erbb2(SP_redcap_export_mapping.SponsoredProjectRunner):

    _SPONSORED_PROJECT = "ERBB2"
    _DATES = ['follow_up_date_int','date_death_int','primary_dx_date_int','lrr_date_int','date_first_met_int',
             'sample_date_int_1','seq_report_date_int_1','sample_date_int_2','seq_report_date_int_2','sample_date_int_3',
             'sequence_report_date_int_3','sample_date_int_4','sequence_report_date_int_4','sample_date_int_5','sequence_report_date_int_5',
             'sample_date_int_6','seq_report_date_int_6','sample_date_int_7','seq_report_date_int_7','sample_date_int_8',
             'sequence_report_date_int_8','sample_date_int_9','sequence_report_date_int_9','sample_date_int_10',
             'sequence_report_date_int_10','date_bso_int','OS_MONTHS','OS_MONTHS_PRIMARY']

    _CASE_LIST_MAF_SAMPLES_TEMPLATE = "cancer_study_identifier: genie_erbb2\nstable_id: genie_erbb2_sequenced\ncase_list_category: all_cases_with_mutation_data\ncase_list_name: Sequenced Tumors\ncase_list_description: All sequenced samples (%s samples)\ncase_list_ids: %s"
    _CASE_LIST_PATH = os.path.join(_SPONSORED_PROJECT,'case_lists')
    _UNMAPPED_SYN_ID = "syn8356977"
    _MAPPED_SYN_ID = "syn8367692"
    _CASE_LIST_SYN_ID = "syn10145925"
    _SP_SYN_ID = "syn8363326"
    _REDCAP_TO_CBIOMAPPING_SYNID = "syn8363731"
    _SP_REDCAP_EXPORTS_SYNID = "syn8322425" #Storage of not found samples
    _NUM_SAMPLE_COLS = 10

    def addOSMonths(self, sponsoredProject_mapped_df):
        #Must add new date fields to the DATE variable along with add to the mapping table: syn8220815
        sponsoredProject_mapped_df['OS_MONTHS'] = sponsoredProject_mapped_df['date_death_int'] - sponsoredProject_mapped_df['date_first_met_int']  
        sponsoredProject_mapped_df['OS_MONTHS_PRIMARY'] = sponsoredProject_mapped_df['date_death_int'] - sponsoredProject_mapped_df['primary_dx_date_int']  
        return(sponsoredProject_mapped_df)

    def createTemporaryGenieId(self, x, tempIdMapping, patientIdCol):
        """
        Create temporary genie id for those that don't have 
        """
        uniqId = x['record_id_patient_id'] + x['redcap_data_access_group']
        if sum(tempIdMapping['uniqueId'] == uniqId) == 0:
            tempId = 'GENIE-%s-%s' % (x['redcap_data_access_group'],''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(10)))
            self.syn.store(synapseclient.Table(self.syn.get("syn10164044"),[[uniqId, tempId]]))
            return(tempId)
        else:
            return(tempIdMapping[tempIdMapping['uniqueId'] == uniqId]['temporaryId'].values[0])
    
    def createNullPatients(self, sponsoredProject_mapped_df, tempIdMappingDf):
        #### TIMELINE FILE
        sponsoredProject_mapped_df['redcap_data_access_group'] = [i.upper() for i in sponsoredProject_mapped_df['redcap_data_access_group']]
        allNullPatients = sponsoredProject_mapped_df[['record_id_patient_id','redcap_data_access_group']][sponsoredProject_mapped_df['record_id_patient_id'].isnull()]
        temporaryIds = allNullPatients.apply(lambda x: self.createTemporaryGenieId(x, tempIdMappingDf, 'record_id_patient_id'), axis =1)
        if not temporaryIds.empty:
            sponsoredProject_mapped_df['record_id_patient_id'][sponsoredProject_mapped_df['record_id_patient_id'].isnull()] = temporaryIds
        assert sum(sponsoredProject_mapped_df['record_id_patient_id'].isnull()) == 0, "Make sure there are no null genie patient Ids"
        sponsoredProject_mapped_df['record_id_patient_id'] = sponsoredProject_mapped_df.apply(lambda x: self.checkGenieId(x, 'redcap_data_access_group','record_id_patient_id'), axis=1)
        return(sponsoredProject_mapped_df, temporaryIds)

    def makeTimeLineDf(self, redCapExportDf, therapyRange = 16):
        START_DATE = []
        STOP_DATE = []
        TREATMENT_TYPE = []
        SUBTYPE = []
        AGENT = []
        THERAPY_RESPONSE = []
        THERAPY_DRUG_OTHER = []
        THERAPY_DRUG_DISCONTINUE = []
        THERAPY_DRUG_REASON = []
        THERAPY_COMBO_YN = []
        THERAPY_COMBO_NUM = []
        ADD_TREATMENT = []
        TREATMENT_SETTING = []
        for therapyNumber in range(1,therapyRange):
            therapyCols = [i for i in redCapExportDf if ("therapy%d_" % therapyNumber in i or "combo_therapy_yn_%d" %therapyNumber == i or "add_treatment_%d" % therapyNumber == i or "treatment_setting_%d" % therapyNumber == i)]
            START_DATE.extend([i for i in therapyCols if "start_int" in i])
            STOP_DATE.extend([i for i in therapyCols if "end_int" in i])
            AGENT.extend([i for i in therapyCols if len(i.split("_")) == 2 and "response" not in i and "ctdrug" not in i])
            THERAPY_DRUG_OTHER.extend([i for i in therapyCols if "other" in i])
            THERAPY_DRUG_DISCONTINUE.extend([i for i in therapyCols if "discon" in i])
            THERAPY_DRUG_REASON.extend([i for i in therapyCols if "reason" in i])
            THERAPY_COMBO_YN.extend([i for i in therapyCols if "combo_therapy_yn" in i] * len([i for i in therapyCols if "start_int" in i]))
            THERAPY_COMBO_NUM.extend([i for i in therapyCols if "combo_num" in i]* len([i for i in therapyCols if "start_int" in i]))
            TREATMENT_TYPE.extend(["Medical Therapy %d" % therapyNumber]* len([i for i in therapyCols if "start_int" in i]))
            SUBTYPE.extend(["Chemo/Target/Immuno etc."] * len([i for i in therapyCols if "start_int" in i]))
            THERAPY_RESPONSE.extend([i for i in therapyCols if "response" in i] *len([i for i in therapyCols if "start_int" in i]))
            ADD_TREATMENT.extend([i for i in therapyCols if "add_treatment" in i] * len([i for i in therapyCols if "start_int" in i]))
            TREATMENT_SETTING.extend([i for i in therapyCols if "treatment_setting" in i] * len([i for i in therapyCols if "start_int" in i]))
        EVENT_TYPE = ["TREATMENT"]*len(AGENT)
        ADD_TREATMENT.extend(['']*4)

        #METASTATIC DIAGNOSIS
        metaDiagnosis = pd.DataFrame()
        #MET DISEASE IS TIMEPOINT 0
        metaDiagnosis['PATIENT_ID'] = redCapExportDf['record_id_patient_id']
        metaDiagnosis['START_DATE'] = 0
        #metaDiagnosis['START_DATE'] = redCapExportDf['date_first_met_int']
        metaDiagnosis['EVENT_TYPE'] = 'STATUS'
        metaDiagnosis['STATUS'] = 'Metastatic Diagnosis'
        metaDiagnosis = metaDiagnosis[~metaDiagnosis['START_DATE'].isnull()]

        removeCols = START_DATE+STOP_DATE+AGENT+THERAPY_DRUG_OTHER+THERAPY_RESPONSE+THERAPY_DRUG_DISCONTINUE+THERAPY_DRUG_REASON+THERAPY_COMBO_YN+THERAPY_COMBO_NUM+ADD_TREATMENT + TREATMENT_SETTING

        lengths = set([
        len(START_DATE),
        len(STOP_DATE),
        len(TREATMENT_TYPE),
        len(SUBTYPE),
        len(AGENT),
        len(THERAPY_RESPONSE),
        len(THERAPY_DRUG_OTHER),
        len(TREATMENT_SETTING),
        len(ADD_TREATMENT),
        len(THERAPY_DRUG_DISCONTINUE),
        len(THERAPY_DRUG_REASON),
        len(THERAPY_COMBO_YN),
        len(THERAPY_COMBO_NUM),
        len(EVENT_TYPE)])
        assert len(lengths) == 1,"Lengths must all be the same"

        total = pd.DataFrame()
        for i in range(len(redCapExportDf)):
            timelineDF = pd.DataFrame()
            timelineDF['PATIENT_ID'] = [redCapExportDf['record_id_patient_id'][i]]*len(START_DATE)
            if not pd.isnull(redCapExportDf.iloc[i]['date_first_met_int']):
                timelineDF['START_DATE'] = [start if pd.isnull(start) else int(start) - int(redCapExportDf.iloc[i]['date_first_met_int']) for start in redCapExportDf.iloc[i][START_DATE].reset_index(drop=True)]
                timelineDF['STOP_DATE'] =  [end if pd.isnull(end) else int(end) - int(redCapExportDf.iloc[i]['date_first_met_int']) for end in redCapExportDf.iloc[i][STOP_DATE].reset_index(drop=True)]
            else:
                timelineDF['START_DATE'] = pd.np.nan
                timelineDF['STOP_DATE'] = pd.np.nan
            timelineDF['EVENT_TYPE'] = EVENT_TYPE
            timelineDF['TREATMENT_TYPE'] = TREATMENT_TYPE
            timelineDF['SUBTYPE'] = SUBTYPE
            timelineDF['AGENT'] = redCapExportDf.iloc[i][AGENT].reset_index(drop=True)
            timelineDF['THERAPY_DRUG_OTHER'] = redCapExportDf.iloc[i][THERAPY_DRUG_OTHER].reset_index(drop=True)
            timelineDF['THERAPY_DRUG_DISCONTINUE'] = redCapExportDf.iloc[i][THERAPY_DRUG_DISCONTINUE].reset_index(drop=True)
            timelineDF['THERAPY_DRUG_REASON'] = redCapExportDf.iloc[i][THERAPY_DRUG_REASON].reset_index(drop=True)
            timelineDF['THERAPY_COMBO_YN'] = redCapExportDf.iloc[i][THERAPY_COMBO_YN].reset_index(drop=True)
            timelineDF['THERAPY_COMBO_NUM'] = redCapExportDf.iloc[i][THERAPY_COMBO_NUM].reset_index(drop=True)
            total = total.append(timelineDF)
        total['STATUS'] = ''
        ordering = total.columns
        total = total.append(metaDiagnosis)
        total = total[ordering]
        return(total, removeCols)

    def getSpecimen(self, getTimelineSpecimen):
        specimen = pd.DataFrame()
        specimen['PATIENT_ID'] = getTimelineSpecimen['PATIENT_ID']
        getTimelineSpecimen = getTimelineSpecimen[~getTimelineSpecimen.SEQUENCE_REPORT_DATE_INT_STATIC.isnull()]
        getTimelineSpecimen = getTimelineSpecimen[~getTimelineSpecimen.METS_DISEASE_DATE_INT.isnull()]
        specimen['START_DATE'] = getTimelineSpecimen.SEQUENCE_REPORT_DATE_INT_STATIC.astype(int) - getTimelineSpecimen.METS_DISEASE_DATE_INT.astype(int)
        specimen['EVENT_TYPE'] = 'SPECIMEN'
        specimen['SAMPLE_ID'] = getTimelineSpecimen['SAMPLE_ID']
        specimen['SAMPLE_NOTES'] = getTimelineSpecimen.SEQUENCE_REPORT_DATE_INT_STATIC
        specimen = specimen[~specimen['START_DATE'].isnull()]
        return(specimen)


class Fgfr4(new_redcap_export_mapping.SponsoredProjectRunner):

    _DATA_ELEMENT_SYN_ID = "syn12032922"
    _SPONSORED_PROJECT = 'FGFR4'
    # No need to define in class
    _CASE_LIST_PATH = os.path.join(_SPONSORED_PROJECT, 'case_lists')
    _NUM_COUNTS = 4
    _REDCAP_TO_CBIOMAPPING_SYNID = "syn15572052"
    _UNLABELLED_SYN_ID = "syn15341849"
    _LABELLED_SYN_ID = "syn15341838"
    # Storage of not found samples
    _SP_REDCAP_EXPORTS_SYNID = "syn11812526"
    _SP_SYN_ID = "syn14721789"
    _CASE_LIST_MAF_SAMPLES_TEMPLATE = (
        "cancer_study_identifier: genie_fgfr4\n"
        "stable_id: genie_fgfr4_sequenced\n"
        "case_list_category: all_cases_with_mutation_data\n"
        "case_list_name: Sequenced Tumors\n"
        "case_list_description: All sequenced samples "
        "(%s samples)\ncase_list_ids: %s")
    _CASE_LIST_SYN_ID = "syn14721794"

    # def addOSMonths(self, sponsoredProject_mapped_df):
    #     '''
    #     Must add new date fields to the DATE variable along with add
    #     to the mapping table: syn8220815
    #     '''
    #     sponsoredProject_mapped_df['OS_MONTHS'] = \
    #         sponsoredProject_mapped_df['death_date_int'] - \
    #         sponsoredProject_mapped_df['date_first_met_int']
    #     sponsoredProject_mapped_df['OS_MONTHS_PRIMARY'] = \
    #         sponsoredProject_mapped_df['death_date_int'] - \
    #         sponsoredProject_mapped_df['primary_dx_date_int']
    #     return(sponsoredProject_mapped_df)

    def makeTimeLineDf(
            self, treatmentDf, finalPatientDf, therapyRange=5):
        # These variables are capitalized to match with the column headers
        START_DATE = []
        STOP_DATE = []
        TREATMENT_TYPE = []
        SUBTYPE = []
        AGENT = []
        THERAPY_RESPONSE = []
        # Name of Chemotherapeutic Agent or Hormone Therapy - Experimental or
        # OTHER (NCIT ID)
        THERAPY_DRUG_OTHER = []
        THERAPY_DRUG_DISCONTINUE = []
        THERAPY_DRUG_REASON = []
        TREATMENT_SETTING = []
        RXNORM_ID = []
        # Name of Chemotherapeutic Agent or Hormone Therapy - Experimental or
        # OTHER
        THERAPY_DRUG_START_ESTIMATED = []
        THERAPY_DRUG_OTHER_NAME = []
        THERAPY_DRUG_END_ESTIMATED = []
        for therapyNumber in range(1, therapyRange):
            therapyCols = [
                i for i in treatmentDf
                if "therapy_drug%d" % therapyNumber in i]
            startCols = [i for i in therapyCols if "start_int" in i]
            START_DATE.extend(startCols)
            STOP_DATE.extend([i for i in therapyCols if "end_int" in i])
            AGENT.extend([
                i for i in therapyCols if "name" in i and "other" not in i])
            RXNORM_ID.extend([
                i for i in therapyCols
                if i == "therapy_drug%d" % therapyNumber])
            THERAPY_DRUG_OTHER.extend([
                i for i in therapyCols if "other" in i and 'name' not in i])
            THERAPY_DRUG_DISCONTINUE.extend([
                i for i in therapyCols if "discon" in i])
            THERAPY_DRUG_REASON.extend([
                i for i in therapyCols if "reason" in i])
            THERAPY_DRUG_OTHER_NAME.extend([
                i for i in therapyCols if "other_name" in i])
            THERAPY_DRUG_START_ESTIMATED.extend([
                i for i in therapyCols if "start_estimated" in i])
            THERAPY_DRUG_END_ESTIMATED.extend([
                i for i in therapyCols if "end_estimated" in i])
            # Value
            TREATMENT_TYPE.extend([
                "Medical Therapy %d" % therapyNumber] * len(startCols))
        # Value
        SUBTYPE = ["Chemo/Target/Immuno etc."] * len(AGENT)
        TREATMENT_SETTING = ['treatment_setting'] * len(AGENT)
        THERAPY_RESPONSE = ['therapy_response'] * len(AGENT)
        # Value
        EVENT_TYPE = ["TREATMENT"]*len(AGENT)
        LINE_START = ['line_start_int'] * len(AGENT)
        REGIMEN_NAME = ['regimen_name'] * len(AGENT)
        CLINICAL_TRIAL = ['clinical_trial'] * len(AGENT)
        CENTER = ['redcap_data_access_group'] * len(AGENT)

        lengths = [
            len(START_DATE),
            len(STOP_DATE),
            len(TREATMENT_TYPE),
            len(AGENT),
            len(THERAPY_DRUG_OTHER),
            len(THERAPY_DRUG_DISCONTINUE),
            len(THERAPY_DRUG_REASON),
            len(RXNORM_ID),
            len(THERAPY_DRUG_OTHER_NAME),
            len(THERAPY_DRUG_START_ESTIMATED),
            len(THERAPY_DRUG_END_ESTIMATED),
            len(TREATMENT_TYPE)]
        assert len(set(lengths)) == 1, "Lengths must all be the same"

        total = pd.DataFrame()
        for i in range(len(treatmentDf)):
            timelineDF = pd.DataFrame()
            timelineDF['PATIENT_ID'] = \
                [treatmentDf['patient_id'].iloc[i]]*len(START_DATE)
            timelineDF['START_DATE'] = \
                treatmentDf.iloc[i][START_DATE].reset_index(drop=True)
            timelineDF['STOP_DATE'] = \
                treatmentDf.iloc[i][STOP_DATE].reset_index(drop=True)

            timelineDF['EVENT_TYPE'] = EVENT_TYPE
            # has to be in this order of PATIENT_ID, START, STOP and EVENT_TYPE
            timelineDF['TREATMENT_TYPE'] = TREATMENT_TYPE
            timelineDF['SUBTYPE'] = SUBTYPE
            timelineDF['AGENT'] = \
                treatmentDf.iloc[i][AGENT].reset_index(drop=True)
            timelineDF['RXNORM_ID'] = \
                treatmentDf.iloc[i][RXNORM_ID].reset_index(drop=True)
            timelineDF['THERAPY_DRUG_OTHER'] = \
                treatmentDf.iloc[i][THERAPY_DRUG_OTHER].reset_index(drop=True)
            timelineDF['THERAPY_DRUG_DISCONTINUE'] = treatmentDf.iloc[i][
                THERAPY_DRUG_DISCONTINUE].reset_index(drop=True)
            timelineDF['THERAPY_DRUG_REASON'] = \
                treatmentDf.iloc[i][THERAPY_DRUG_REASON].reset_index(drop=True)
            timelineDF['THERAPY_DRUG_OTHER_NAME'] = treatmentDf.iloc[i][
                THERAPY_DRUG_OTHER_NAME].reset_index(drop=True)
            timelineDF['THERAPY_DRUG_START_ESTIMATED'] = treatmentDf.iloc[i][
                THERAPY_DRUG_START_ESTIMATED].reset_index(drop=True)
            timelineDF['THERAPY_DRUG_END_ESTIMATED'] = treatmentDf.iloc[i][
                THERAPY_DRUG_END_ESTIMATED].reset_index(drop=True)
            timelineDF['TREATMENT_SETTING'] = \
                treatmentDf.iloc[i][TREATMENT_SETTING].reset_index(drop=True)
            timelineDF['THERAPY_RESPONSE'] = \
                treatmentDf.iloc[i][THERAPY_RESPONSE].reset_index(drop=True)
            timelineDF['LINE_START'] = \
                treatmentDf.iloc[i][LINE_START].reset_index(drop=True)
            timelineDF['REGIMEN_NAME'] = \
                treatmentDf.iloc[i][REGIMEN_NAME].reset_index(drop=True)
            timelineDF['CLINICAL_TRIAL'] = \
                treatmentDf.iloc[i][CLINICAL_TRIAL].reset_index(drop=True)
            timelineDF['CENTER'] = \
                treatmentDf.iloc[i][CENTER].reset_index(drop=True)
            total = total.append(timelineDF, sort=False)
        # remove all without START dates
        total = total[~total['START_DATE'].isnull()]
        total['SP'] = self._SPONSORED_PROJECT
        total['STATUS'] = ''
        total['START_DATE'] = total['START_DATE'].astype('float')
        total['STOP_DATE'] = total['STOP_DATE'].astype('float')
        total['RXNORM_ID'] = total['RXNORM_ID'].astype('float')
        total['LINE_START'] = total['LINE_START'].astype('float')
        total.drop_duplicates(inplace=True)
        # Anchor point is MET_DX_DATE_INT
        date_met_int = [
            float(finalPatientDf['MET_DX_DATE_INT'][
                finalPatientDf['PATIENT_ID'] == patient].values[0])
            for patient in total['PATIENT_ID']]
        total['START_DATE'] = total['START_DATE'] - date_met_int
        total['STOP_DATE'] = total['STOP_DATE'] - date_met_int
        total['LINE_START'] = total['LINE_START'] - date_met_int

        return(total)

    def createSpecimenDf(self, sampleDf, patientDf):
        clinicalDf = sampleDf.merge(patientDf, on="PATIENT_ID", how="outer")
        clinicalDf = clinicalDf[~clinicalDf.AGE_AT_SEQ_REPORT.isnull()]
        clinicalDf = \
            clinicalDf[~clinicalDf.DATE_FIRST_DISTANT_MET_INT.isnull()]
        specimen = pd.DataFrame()
        specimen['PATIENT_ID'] = clinicalDf['PATIENT_ID']
        specimen['SAMPLE_ID'] = clinicalDf['SAMPLE_ID']
        specimen['START_DATE'] = \
            clinicalDf.AGE_AT_SEQ_REPORT.astype(int) - \
            clinicalDf.DATE_FIRST_DISTANT_MET_INT.astype(int)
        specimen['EVENT_TYPE'] = 'SPECIMEN'
        specimen['SAMPLE_NOTES'] = clinicalDf.AGE_AT_SEQ_REPORT
        specimen = specimen[~specimen['START_DATE'].isnull()]
        return(specimen)

