#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.cohort = 'NSCLC' // Default
params.release = '1.1-consortium'  // Default

// Check if cohort is part of allowed cohort list
def allowed_cohorts = ["BLADDER", "BrCa", "CRC", "NSCLC", "PANC", "Prostate"]
if (!allowed_cohorts.contains(params.cohort)) {exit 1, 'Invalid cohort name'}

ch_cohort = Channel.value(params.cohort)
ch_release = Channel.value(params.release)

/*
Run cBioPortal Export
*/
process cBioPortalExport {
   container 'sagebionetworks/geniesp'
   secret 'SYNAPSE_AUTH_TOKEN'

   input:
   val cohort
   val release

   output:
   stdout

   script:
   """
   geniesp $cohort $release --staging
   """
}


workflow {
   cBioPortalExport(ch_cohort, ch_release)
}