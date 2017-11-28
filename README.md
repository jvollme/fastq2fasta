### fastq2fasta.py

convert fastq files into fasta (ignoring quality score information). Can read and write compressed files.

#### dependencies:
 - python 2.7+
 - biopython

#### usage: 

````
fastq2fasta.py [-h] -if INPUT_FASTQ [-of OUTPUT_FASTA]
               [-ct {gzip,gz,bzip2,bz2,zip,none}] [-V]
````

**Warning**: No quality information will be retained for the resulting fasta-file

#### optional arguments:
````
  -h, --help            show this help message and exit
  -if INPUT_FASTQ, --in_fastq INPUT_FASTQ
                        Input fastq file
  -of OUTPUT_FASTA, --out_fasta OUTPUT_FASTA
                        Output fasta file (Default=<input_fastq>.fasta
  -ct {gzip,gz,bzip2,bz2,zip,none}, --compression {gzip,gz,bzip2,bz2,zip,none}
                        type of compression ('gz', 'bz2', 'zip', 'none').
                        default= guess compression-format based on filename
                        extension
  -V, --version         show program's version number and exit
````
