Notes-from-Nature
=================
## Author
Jun Ying Lim
Department of Integrative Biology and UC Museum of Paleontology, University of California, Berkeley

## Description
* Reconciles output from citizen science transcription efforts from Notes from Nature (http://notesfromnature.org) for the Calbug project (http://calbug.berkeley.edu/)
* Focuses on specimens from the Essig Museum of Entomology, University of California, Berkeley.

## Installation
### Python packages
You will require a number of python packages to run this. You can do this by using the pip installer. If you're running Python 3, pip should be by default installed. Input the following in Terminal

`pip3 install fuzzywuzzy numpy pandas nltk python-Levenshtein`

### Other binaries
The string consensus functions use the sequence alignment algorithms provided in MAFFT. You will have to install mafft into the directory `/usr/local/bin`

If you have homebrew installed, you can install it as so:
`brew install mafft`

## Usage
To run TranscriptResolver in Terminal, you simply need to <code> cd </code> to the folder containing TranscriptResolver and:

`python3 transcriptResolver.py`

The program will prompt you on the necessary parameters. You can also supply them directly in terminal.

`python3 transcriptResolver.py -wd ~/Dropbox/Projects/2013/Calbug/TranscriptResolver/tests -file test_file.csv -stem test -col_id UNIQUE_ID -col_target [DATE,TARGET] -col_method [vote_count,metadata]`


## Contact
Feel free to contact me if you have any questions!
