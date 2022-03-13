#!/bin/bash
### Workflow script of the data analysis with wget link downloading the original data files
### optimized to run on MacOS High Sierra 10.13.6, using: R 4.0.3
### created by Marketa Vlkova in 19-September-2021
### example to run: ./Run.sh

### STEP 1:
### get genome files of environmental strains from Figshare
wget --content-disposition https://figshare.com/ndownloader/files/30817357
tar -xzvf PhenAnalysis.tar.gz
rm PhenAnalysis.tar.gz

### STEP 2:
### produce Figures (3 - 7 and Supp 3 - 4)
./PlotPhenotype.R
