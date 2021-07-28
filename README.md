# Signal peptide prediction
Python script to run 4 signal peptide prediction software and generate a consensus

This python script will run any combination of 4 different signal peptide predictors (SignalP, WolfPsort, Phobius and Deepsig) and generate a consensus list of proteins that are signal peptides based on the results of all chosen predictors.

# Requirements
* matplotlib
* python3
* [pyvenn](https://github.com/tctianchi/pyvenn)
* [Phobius](https://phobius.sbc.su.se/data.html)
* [SignalP](https://services.healthtech.dtu.dk/software.php)
* [WolfPsort](https://github.com/fmaguire/WoLFPSort)
* [Deepsig](https://github.com/BolognaBiocomp/deepsig)

## Quick usage

This script can be run with the following, where -i is the input protein file in fasta format and -s -w -p and -d for running SignalP, WolfPsort, Phobius and Deepsig.

`./SP_prediction.py -i protein.fa -swpd`

A folder in the location where SP_prediction.py is run, will be created that contains output folders for each prediction program as well as a folder for venn diagrams if more than one predictor is run.

## Notes
Currently, this is written for using the four predictors in my virtual environment. Only SignalP and Deepsig are in my path. Lines 33, 47, 63 and 81 contain the code for running the predictors and these will likely need to be changed for running in your own environment. I may update this in the future.

Included in this repository is also a script to extract the final proteins from the consensus list and write them to a new file. I may also incorporate this in a future update. It's usage is as follows `./protein_extract.py <consensus file> <protein file> <output file>` and will require Biopython.


## Full options and usage
```
usage: Signal_peptide_runner.py -i INPUT [-p] [-s] [-d] [-w] [-h]

Signal Peptide Prediction

Required Arguments:
  -i INPUT, --input INPUT

Optional Arguments:
  -p, --phobius         Run Phobius
  -s, --signalp         Run SignalP
  -d, --deepsig         Run Deepsig
  -w, --wolfpsort       Run WolfPSort
  -h, --help            show this help message and exit
```
