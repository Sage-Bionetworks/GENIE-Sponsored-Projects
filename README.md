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
docker run --rm -e SYNAPSE_AUTH_TOKEN='<PAT>' geniesp PANC 1.1-consortium --staging
```

To load the files directly to Synapse, remove the --staging parameter
```
docker run --rm -e SYNAPSE_AUTH_TOKEN='<PAT>' geniesp PANC 1.1-consortium
```

## Scripts

To validate a cBioPortal mapping file stored on synapse:
```
python validate_map.py --log info --synapse_id {synapse_id}
```

or stored in a local file:
```
python validate_map.py --log info --synapse_id {/path/to/file.csv}
```

To view full usage details:
```
python validate_map.py -h
```

which outputs:
```
Checks validity of BPC to cBioPortal mapping file

optional arguments:
  -h, --help            show this help message and exit
  --synapse_id SYNAPSE_ID, -s SYNAPSE_ID
                        Synapse ID of mapping file
  --file FILE, -f FILE  Local path to mapping file
  --version VERSION, -v VERSION
                        Synapse entity version number (default: current)
  --outfile OUTFILE, -o OUTFILE
                        Name of output file (default: output.csv)
  --log {debug,info,warning,error}, -l {debug,info,warning,error}
                        Set logging output level (default: error)
```
