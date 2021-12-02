# GENIE-Sponsored-Projects
This repository will contain processing code for GENIE sponsored projects used to for creating cBioPortal files


## Installation

Build Docker container to set up your python environment

```
docker build -t geniesp .
```

Note: depending on the environment in which you are running Docker, prefixing docker commands with 'sudo' may be necessary.  

## Usage
```
docker run --rm geniesp -h
usage: __main__.py [-h] [--synconfigpath SYNCONFIGPATH] [--staging]
                   {NSCLC,CRC,BrCa,PANC,Prostate,AKT1,ERRB2,FGFR4} release
                   cbiopath

Run GENIE sponsored projects

positional arguments:
  {NSCLC,CRC,BrCa,PANC,Prostate,AKT1,ERRB2,FGFR4}
                        Specify sponsored project to run
  release               Specify bpc release
  cbiopath              Specify path to cbio: must do `git clone
                        https://github.com/cBioPortal/cbioportal.git`

optional arguments:
  -h, --help            show this help message and exit
  --synconfigpath SYNCONFIGPATH
                        Specify path to .synapseConfig file
  --staging             If true, files aren't uploaded onto synapse
  ```