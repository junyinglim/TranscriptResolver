TranscriptResolver
=================
# Author
Jun Ying Lim

Department of Integrative Biology and UC Museum of Paleontology, University of California, Berkeley

# Description

### TranscriptResolver
* Reconciles output from citizen science transcription efforts from Notes from Nature (http://notesfromnature.org) for the Calbug project (http://calbug.berkeley.edu/)
* Code should work for any type of data where there are replicates (e.g., replicate verbatim transcriptions of the same specimen label)
* Uses sequence alignment algorithms to find consensus strings (http://blog.notesfromnature.org/2014/01/14/checking-notes-from-nature-data/)

Different consensus methods:
* `vote_count` - Chooses the most frequently occurring value. Recommended for fields where choice of values is constrained (e.g., drop down lists)
* `consensus` - Implements a character sequence alignment on replicate strings and produces a consensus string. Recommended for fields where input is more free-style (e.g., verbatim transcription of fields)
* `metadata` - Does not perform any consensus method per se. Instead combines all values into a single string, delimited by "|"

### TranscriptPrepare
* Custom script to prepare raw Notes from Nature output for resolving using TranscriptResolver
* Main steps:
  * Normalizes or prepare fields 
  * Cleans up column headers
  * Creates a file for bulk upload into the Essig Database


### TranscriptClean
* Custom script to prepare raw Notes from Nature output for resolving using TranscriptResolver
* Main steps:
  * Normalizes or prepare fields 
  * Cleans up column headers
  * Creates a file for bulk upload into the Essig Database


# Installation
### Python packages
You will require a number of python packages to run this. You can do this by using the pip installer. If you're running Python 3, pip should be by default installed. Input the following in Terminal

`pip3 install fuzzywuzzy numpy pandas nltk python-Levenshtein pymysql biopython`

### Other binaries
The string consensus functions use the sequence alignment algorithms provided in MAFFT. You will have to install mafft into the directory `/usr/local/bin`

If you have homebrew installed, you can install it as so:
`brew install mafft`


# Usage
To run TranscriptResolver in Terminal, you simply need to <code> cd </code> to the folder containing TranscriptResolver and:

`python3 transcriptResolver.py`

The program will prompt you on the necessary parameters. You can also supply them directly in terminal.

`python3 transcriptResolver.py -wd <yourworkingdir> -file <yourfile> -stem <yourstemname> -col_id UNIQUE_ID -col_target [<field1>,<field2>] -col_method [<method1>,<method2>]`


# Contact
Feel free to contact me (junyinglim<at>berkeley.edu) if you have any questions or need help with installation!
