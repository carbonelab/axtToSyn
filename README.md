# axtToSyn

Detect synteny blocks and synteny breakpoints from pairwise genome alignment file.

axtToSyn.py is a python script that performs synteny block detection using a pairwise genome alignment file output by the UCSC pairwise genome alignment pipeline. 

## Install

Dependencies: 

 - Python >= 3.6

```
wget https://raw.githubusercontent.com/carbonelab/axtToSyn/master/axtToSyn.py
chmod +x axtToSyn.py
./axtToSyn.py -h
```

## Usage

```
$ ./axtToSyn.py -h
usage: axtToSyn [-h] [--min-score [MIN_SCORE]] [--min-blen [MIN_BLEN]]
                file outfile tname qname

Generates synteny blocks from pairwise genome alignment .axt file by block
elongation.

positional arguments:
  file                  Relative path to net.axt alignment file.
  outfile               Path to output synteny blocks file (default:
                        synblocks.csv)
  tname                 Target assembly, used to name breakpoints.
  qname                 Query assembly, used to name breakpoints.

optional arguments:
  -h, --help            show this help message and exit
  --min-score [MIN_SCORE]
                        Min alignment score of block to be considered for
                        elongation (defalt: 1e5)
  --min-blen [MIN_BLEN]
                        Min block len to be considered for elongation in the
                        first pass (defalt: 1e3)
``` 
