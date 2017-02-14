CirQuestDSV

ptdolan@stanford.edu
===========
Description:

This is a repository for code related to the analysis of Circular Resequencing data. Code may be written in python and/or R.
All code in this repository is intended for use downstream of Q-threshold count table generation by the CirSeq software package.

This package has been updated with a base program for analyzing the output of CirSeq ('"Q" files') and plotting general information about them. It also generate a list of high frequency mutations and agglomerates q files into useful R data structures for further analysis. Updated to python 3.0. 

Please fork this project to improve annotation in the python portion ("Q20Analysis.py") or to develop analyses based on this data structure in R.


Requirements:
Python (3.X compliant)
R

Usage:

File names:
 1. All input file names must end in "Q20.txt" or "Q20threshold.txt"!
 	This makes sure that overwriting does not occur when output files and input files are in the same directory. 

 2. Each input file requires specific naming for grouping and plotting features. these will also be parsed for data table annotation.
	e.g.: "~/pathto/myData/MouseA_plusTreatment_1-Q20.txt"

 3. For grouped passages:
 	In order for files to be grouped correctly in the plotting and analysis steps, files in directory must be labeled "MYHEADERTEXT_1-Q20.txt","HEADERTEXT_2-Q20.txt","HEADERTEXT_3-Q20.txt", those that do not match criteria are excluded from trajectory plots. Including MYHEADERTEXT is optional but useful for labeling if no header needed, "1-Q20.txt","2-Q20.txt","3-Q20.txt" is sufficient for grouping. 

Current usage instructions:

1. Open a terminal in the same directory as the scripts and run: 

    > python Q20Analysis.py \<directory with q20 files> \<translation start> \<translation end> \<next ORF start> \<next ORF stop> and so on...

    This generates an annotated q20 ("-q20annot.txt") file. You can stop here and use this annotated file for interpretation, or you can use this enhanced q20 to generate plots with the associated R script, "AnnotAnalysis.R". 
    
2. "AnnotAnalysis.R". This program reads a directory of annotated files and combines them to generate frequency plots, trajectory plots, it also performs dimension reduction analysis and outputs a number of useful output files. 

with R installed: 
	in command line: 
		Rscript ./path/to/AnnotAnalysis.R ./pathto/q20annots/

	or: 
		R
	then, in R: 
		source("./path/to/AnnotAnalysis.R","./pathto/q20annots/)

_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-


MIT License

Copyright (c) [2016-2017] [Patrick T. Dolan]

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
