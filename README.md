# GENIE-Sponsored-Projects
This repository will contain processing code for GENIE sponsored projects used to for creating cBioPortal files


## Installation

1. Clone this repository and navigate to the directory
```
git clone 
cd GENIE-Sponsored-Projects
```

2. Once inside the directory, clone the cbioportal repository
```
git clone https://github.com/cBioPortal/cbioportal.git
```

3. and copy Synapse credentials from home directory
```
cp ~/.synapseConfig .
```

4. Build the Docker container to set up your python environment

```
docker build -t geniesp .
```

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

optional arguments:
  -h, --help            show this help message and exit
  --staging             If true, files aren't uploaded onto synapse
  ```