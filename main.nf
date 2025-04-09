#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
Run cBioPortal Export
*/
process cBioPortalExport {
   container "$params.geniesp_docker"
   secret 'SYNAPSE_AUTH_TOKEN'

   input:
   val cohort
   val release
   val production
   val use_grs
   val cpt_seq_date_replacement_type

   output:
   stdout

   script:
   if (production && use_grs) {
      """
      geniesp $cohort $release \
         --upload \
         --cbioportal /usr/src/cbioportal \
         --production \
         --use-grs \
         --cpt_seq_date_replacement_type $cpt_seq_date_replacement_type
      """
   } else if (production && !use_grs){
      """
      geniesp $cohort $release \
         --upload \
         --cbioportal /usr/src/cbioportal \
         --production \
         --cpt_seq_date_replacement_type $cpt_seq_date_replacement_type
      """
   } else if (!production && use_grs){
      """
      geniesp $cohort $release \
         --upload \
         --cbioportal /usr/src/cbioportal \
         --use-grs \
         --cpt_seq_date_replacement_type $cpt_seq_date_replacement_type
      """
   } else {
      """
      geniesp $cohort $release \
         --upload \
         --cbioportal /usr/src/cbioportal \
         --cpt_seq_date_replacement_type $cpt_seq_date_replacement_type
      """
   }
}

workflow {
   params.cohort = 'NSCLC' // Default
   params.release = '1.1-consortium'  // Default
   params.production = false
   params.use_grs = false
   params.cpt_seq_date_replacement_type = 'derived_variable'

   // Check if cohort is part of allowed cohort list
   def allowed_cohorts = ["BLADDER", "BrCa", "CRC", "ESOPHAGO", "MELANOMA", "NSCLC", "OVARIAN", "PANC", "Prostate", "RENAL"]
   if (!allowed_cohorts.contains(params.cohort)) {exit 1, 'Invalid cohort name'}

   ch_cohort = Channel.value(params.cohort)
   ch_release = Channel.value(params.release)
   ch_production = Channel.value(params.production)
   ch_use_grs = Channel.value(params.use_grs)
   ch_cpt_seq_date_replacement_type = Channel.value(params.cpt_seq_date_replacement_type)

   cBioPortalExport(ch_cohort, ch_release, ch_production, ch_use_grs, ch_cpt_seq_date_replacement_type)
}
