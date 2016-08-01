CirQuestDSV

ptdolan@stanford.edu
===========
Description:

This is a repository for code related to the analysis of Circular Resequencing data. Code may be written in python  and/or R.
All code in this repository is intended for use downstream of Q-threshold count table generation by the CirSeq software package.

This package has been updated with a base program for analyzing the output of CirSeq ('"Q" files') and plotting general information about them. It also generate a list of high frequency mutations and agglomerates q files into useful R data structures for further analysis. Updated to python 3.0. 

Please fork this project to improve annotation in the python portion ("Q20Analysis.py") or to develop analyses based on this data structure in R.


Requirements:
Python (3.X compliant)


/!5. Usage:

File names:
 1. All file names must end in "Q20.txt" !!!
 2. For grouped passages: Files in directory must be labeled "1-Q20.txt","2-Q20.txt","3-Q20.txt", those that do not match criteria are excluded from trajectory plots.

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


4. MIT License

Copyright (c) [2016] [Patrick T. Dolan]

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
