# GENIE-Sponsored-Projects
This repository will contain processing code for GENIE sponsored projects used to for creating cBioPortal files


## Installation

1. Install Docker: https://docs.docker.com/get-docker/

2. Clone this repository and navigate to the directory
```
git clone git@github.com:Sage-Bionetworks/GENIE-Sponsored-Projects.git
cd GENIE-Sponsored-Projects
```

3. Build the container 
```
docker build -t geniesp .
```

## Synapse credentials

Cache your Synapse personal access token (PAT) as an environmental variable:
```
export SYNAPSE_AUTH_TOKEN={your_personal_access_token_here}
```

## Usage

To view usage details, run
```
docker run --rm geniesp -h
```

Output will be as follows

```
usage: geniesp [-h] [--staging]
               {NSCLC,CRC,BrCa,PANC,Prostate,AKT1,ERRB2,FGFR4} release

Run GENIE sponsored projects

positional arguments:
  {NSCLC,CRC,BrCa,PANC,Prostate,AKT1,ERRB2,FGFR4}
                        Specify sponsored project to run
  release               Specify bpc release

optional arguments:
  -h, --help            show this help message and exit
  --staging             If true, files aren't uploaded onto synapse
```

Example command line for 'PANC' (pancreas) cohort and the '1.1-consortium' release without loading files to Synapse.  Only the output to the console will be accessible.
```
docker run --rm -e SYNAPSE_AUTH_TOKEN=$SYNAPSE_AUTH_TOKEN geniesp PANC 1.1-consortium --staging
```

To load the files directly to Synapse, remove the --staging parameter
```
docker run --rm -e SYNAPSE_AUTH_TOKEN=$SYNAPSE_AUTH_TOKEN geniesp PANC 1.1-consortium
```

## Scripts

To validate a cBioPortal mapping file stored on synapse:
```
python validate_map.py -s {synapse_id} -c {cohort} -r {release} -l info
```

or stored in a local file:
```
python validate_map.py -f {/path/to/file.csv} -c {cohort} -r {release} -l info
```

To view full usage details:
```
python validate_map.py -h
```

which outputs:
```
usage: validate_map.py [-h] [--synapse_id SYNAPSE_ID | --file FILE] [--version VERSION] [--cohort {BLADDER,BRCA,CRC,NSCLC,PANC,PROSTATE}]
                       [--release {1.1-consortium,1.2-consortium,2.0-public,2.1-consortium}] [--outfile OUTFILE] [--log {debug,info,warning,error}]

Checks validity of BPC to cBioPortal mapping file

optional arguments:
  -h, --help            show this help message and exit
  --synapse_id SYNAPSE_ID, -s SYNAPSE_ID
                        Synapse ID of mapping file
  --file FILE, -f FILE  Local path to mapping file
  --version VERSION, -v VERSION
                        Synapse entity version number (default: current)
  --cohort {BLADDER,BRCA,CRC,NSCLC,PANC,PROSTATE}, -c {BLADDER,BRCA,CRC,NSCLC,PANC,PROSTATE}
                        BPC cohort label (default: BLADDER)
  --release {1.1-consortium,1.2-consortium,2.0-public,2.1-consortium}, -r {1.1-consortium,1.2-consortium,2.0-public,2.1-consortium}
                        Release label (default: 1.1-consortium)
  --outfile OUTFILE, -o OUTFILE
                        Name of output file (default: output.csv)
  --log {debug,info,warning,error}, -l {debug,info,warning,error}
                        Set logging output level (default: error)
```

## Troubleshooting
The most common issues when running GENIE-Sponsored-Projects code for BPC involve changes to variable names of the underlying source data and outdated or incorrect references.  

Variable references
1. syn22294851: the Scope of Release maintains a running log of changes to derived variable names  
2. syn17011602: curated variable names are listed in data dictionaries for each cohort

If obtaining an error that a variable cannot be found, check and update the following references on Synapse:
1. syn21431364: GENIE BPC No PHI Data Elements Catalog
1. syn25712693: BPC REDCap to cbio mapping
1. syn22296821: Dataset labels in Data files for derived variables

On rare occassions, variable name changes may also require changes in the codebase:
1. Check for any hardcoded variable names have been updated

Finally, some information is not collected for particular cohorts.  
1. syn20852283: check if data is collected for a given cohort by investigating raw file uploads
