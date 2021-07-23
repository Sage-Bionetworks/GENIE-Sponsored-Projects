# GENIE-Sponsored-Projects
This repository will contain processing code for GENIE sponsored projects used to for creating cBioPortal files


## Installation

Use virtual env or conda to set up your python environment

```
conda activate geniesp
pip install .
pip install -r requirements
```

## Usage
```
geniesp -h
usage: geniesp [-h] [--staging] {NSCLC,CRC,BrCa,AKT1,ERRB2,FGFR4} cBioPath release

Run GENIE sponsored projects

positional arguments:
  {NSCLC,CRC,BrCa,AKT1,ERRB2,FGFR4}
                        Specify sponsored project to run
  cBioPath              Specify path to cbio: must do `git clone https://github.com/cBioPortal/cbioportal.git`
  release               Specify bpc release

optional arguments:
  -h, --help            show this help message and exit
  --staging             If true, files aren't uploaded onto synapse
```