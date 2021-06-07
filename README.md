# fast5_demultiplexer
Demultiplex fast5 reads based on Guppy output

## Description
Reorganize fast5 file by pass/fail and by barcode. Duplicates Guppy output tree.

## TODO
- Start with the multi fast5 folder from MinKNOW. Add function to convert multi to single.
- List reads names form files in parallel
- Convert demultiplexed fast5 from single to multi in parallel

## Installation
1- Create and activate conda environment (must have conda installed):
```
conda create -n fast5_demultiplexer python=3 pip
conda activate fast5_demultiplexer
pip install ont-fast5-api
```
2- Clone repository and test pipeline:
```
git clone https://github.com/duceppemo/fast5_demultiplexer
cd fast5_demultiplexer
python fast5_demultiplexer.py -h
```

## Usage
```
usage: fast5_demultiplexer.py [-h] -b /basecalled/ -s /single_fast5/ -d /demultiplexed_fast5 [-t 16]

Reorganize fast5 single files according to Guppy demultiplexing results. WARNING: about times the disk footprint of the fast5 folder is required!

optional arguments:
  -h, --help            show this help message and exit
  -b /basecalled/, --basecalled /basecalled/
                        "Save_path" folder used for Guppy.
  -s /single_fast5/, --singles /single_fast5/
                        parent folder where all the single fast5 files are located.
  -d /demultiplexed_fast5, --demultiplexed /demultiplexed_fast5
                        Folder where to move sinfle fast5 files according to Guppy demultiplexing.
  -t 16, --threads 16   Number of threads to use. Default max available.
```
