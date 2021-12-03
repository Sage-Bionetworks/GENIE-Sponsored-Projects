# GENIE-Sponsored-Projects
This repository will contain processing code for GENIE sponsored projects used to for creating cBioPortal files


## Installation

1. Clone this repository and navigate to the directory
```
git clone git@github.com:Sage-Bionetworks/GENIE-Sponsored-Projects.git
# or 
git clone https://github.com/Sage-Bionetworks/GENIE-Sponsored-Projects.git

cd GENIE-Sponsored-Projects
```

2. Once inside the directory, clone the cbioportal repository
```
git clone git@github.com:cBioPortal/cbioportal.git
# or
git clone https://github.com/cBioPortal/cbioportal.git
```

3. and copy Synapse credentials from your home directory (or wherever they are cached)
```
cp ~/.synapseConfig .
```

4. Build the Docker container to set up your python environment

```
docker build -t geniesp .
```

## Usage

To view usage details, run
```
docker run --rm geniesp -h
```

Output will be as follows

```
usage: __main__.py [-h] [--staging]
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
docker -rm run geniesp PANC 1.1-consortium --staging
```

To load the files directly to Synapse, remove the --staging parameter
```
docker run -rm geniesp PANC 1.1-consortium 
```
