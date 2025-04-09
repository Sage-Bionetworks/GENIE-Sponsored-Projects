# GENIE-Sponsored-Projects
This repository will contain processing code for GENIE sponsored projects used to for creating cBioPortal files


## Installation

### Using Python

1. Close the repository and navigate to the `READMD.md` in your local directory

3. Build the environment using `pip install -e .`

### Using Docker

1. Install Docker: https://docs.docker.com/get-docker/

2. Clone this repository and navigate to the directory
```
git clone git@github.com:Sage-Bionetworks/GENIE-Sponsored-Projects.git
cd GENIE-Sponsored-Projects
```

3. Build the container locally.
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
               {NSCLC,CRC,BrCa,PANC,Prostate,AKT1,ERRB2,FGFR4,ESOPHAGO,MELANOMA,OVARIAN,RENAL} release

Run GENIE sponsored projects

positional arguments:
  {NSCLC,CRC,BrCa,PANC,Prostate,AKT1,ERRB2,FGFR4,ESOPHAGO,MELANOMA,OVARIAN,RENAL}
                        Specify project to run
  release               Specify bpc release (e.g. 1.1-consortium)

optional arguments:
  -h, --help                        Show this help message and exit
  --upload                          Upload files into Synapse BPC staging directory. Default: false
  --cbioportal {synapseID}
                                    Optional parameter to specify cbioportal folder
                                    location
  --production                      Whether to run in production mode or not. Default: false
  --use-grs                         Whether to use grs or use dd as primary mapping. Default: false
  --cpt_seq_date_replacement_type   The replacement data type to use for the cpt_seq_date value replacement.
                                     Default: derived_variable
```

Example command line:

This runs the release pipeline for BLADDER 1.1 in non-production mode (staging) with GRS enabled using
main genie clinical data to do the cpt_seq_date replacement

```
geniesp BLADDER 1.1-consortium --upload --use-grs --cpt_seq_date_replacement_type main_genie
```

Example command using docker:

This runs the release pipeline for PANC 1.1 in non-production mode (staging) using derived variable data to
do the cpt_seq_date replacement
```
docker run --rm -e SYNAPSE_AUTH_TOKEN=$SYNAPSE_AUTH_TOKEN geniesp geniesp PANC 1.1-consortium --upload
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
usage: validate_map.py [-h] [--synapse_id SYNAPSE_ID | --file FILE] [--version VERSION] [--cohort COHORT]
                       [--release {1.1-consortium,1.2-consortium,2.0-public,2.1-consortium}] [--outfile OUTFILE] [--log {debug,info,warning,error}]

Checks validity of BPC to cBioPortal mapping file

optional arguments:
  -h, --help            show this help message and exit
  --synapse_id SYNAPSE_ID, -s SYNAPSE_ID
                        Synapse ID of mapping file
  --file FILE, -f FILE  Local path to mapping file
  --version VERSION, -v VERSION
                        Synapse entity version number (default: current)
  --cohort COHORT, -c COHORT
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
1. syn25712693: BPC REDCap to cbio mapping
1. syn22296821: Dataset labels in Data files for derived variables

On rare occassions, variable name changes may also require changes in the codebase:
1. Check for any hardcoded variable names have been updated

Finally, some information is not collected for particular cohorts.  
1. syn20852283: check if data is collected for a given cohort by investigating raw file uploads
