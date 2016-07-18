CirQuestDSV
(C) Copyright 2014-2016. All Rights Reserved.
ptdolan@stanford.edu
===========
12-13-14 

This is a repository for code related to the analysis of Circular Resequencing data. Code may be written in python  and/or R.
All code in this repository is intended for use downstream of Q-threshold count table generation by the CirSeq software package.

7-13-16

This package has been updated with a base program for analyzing the output of CirSeq ('"Q" files') and plotting general information about them. It also generate a list of high frequency mutations and agglomerates q files into useful R data structures for further analysis. Updated to python 3.0.

Please fork this project to improve annotation in the python portion ("Q20Analysis.py") or to develop analyses based on this data structure in R.

Current usage instructions:
1. Open a terminal in the same directory as the scripts and run: 

    > python Q20Analysis.py <directory -- all file names must end in "Q20.txt" > <translation start> <translation end> <next ORF start> <next ORF stop> and so on...

2. Edit 'annotAnalysis.R' to show the same directory as the python argument above.

3. Run script in R. 

3a. Open and run in RStudio or:

3b. In command line:
    > r 
then:
    > source("AnnotAnalysis.R")




