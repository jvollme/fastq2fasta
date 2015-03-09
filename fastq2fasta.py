#!/usr/bin/python
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os.path, gzip, bz2, zipfile, sys
from Bio import SeqIO
import argparse

version = "0.1 (14.12.14)"
name = "fastq2fasta.py"
#created 08.12.2014 by John Vollmers
parser=argparse.ArgumentParser(description = "convert fastq format into fasta. version = " + version)
parser.add_argument('-if', '--in_fastq', action = "store", dest = "input_fastq", required = True, help = "Input fastq file")
parser.add_argument('-of', '--out_fasta', action = "store", dest = "output_fasta", default = None, help = "Output fasta file (Default=<input_fastq>.fasta")
parser.add_argument('-ct','--compression', action = "store", dest = "compression_type", choices = ["gzip","gz", "bzip2","bz2", "zip", "none"], help = "type of compression ('gz', 'bz2', 'zip', 'none'). default= guess compression-format based on filename extension")
parser.add_argument('-V', '--version', action="version", version=name + " version " + version) 
args = parser.parse_args()

if args.compression_type != None:
	compression_type = args.compression_type
else: #if no compression format explicitely provided, try to guess it from filename
	if args.input_fastq.endswith(".gz"):
		compression_type = "gz"
	elif args.input_fastq.endswith(".bz2"):
		compression_type = "bzip2"
	elif args.input_fastq.endswith(".zip"):
		compression_type = "zip"
	else:
		compression_type = "none"

def readwrite_fasta(infilename, outfilename):
	#	print "reading fasta"
	try:
		if compression_type in ["gzip", "gz"]:
			infile = gzip.open(infilename, 'r')
		elif compression_type in ["bzip2", "bz2"]:
			infile = bz2.BZ2File(infilename, 'r')
		elif compression_type == "zip":
			myzipfile = zipfile.ZipFile(infilename, 'r')
			if len(myzipfile.namelist()) > 1:
				raise IOError, "TOO MANY FILES IN ZIPFILE"
			else:
				infile = myzipfile.open(myzipfile.namelist()[0])
		else:
			infile = open(infilename, 'r')
		outfile=open(outfilename, "w")
		input_seq_iterator=SeqIO.parse(infile, "fastq")
		recordcounter = 0
		for record in input_seq_iterator:
			SeqIO.write(record, outfile, "fasta")
			recordcounter += 1
			sys.stdout.write("\rprocessed " + str(recordcounter) + " reads")
			sys.stdout.flush()
		outfile.close()
		infile.close()
	except Exception, ex:
		print ex.__class__.__name__ + " : " + str(ex)
		return None

input_fastq = args.input_fastq
if args.output_fasta == None:
	output_fasta = input_fastq.rstrip(".gz") + ".fasta"
else:
	output_fasta = args.output_fasta

readwrite_fasta(input_fastq, output_fasta)
print "Wrote sequences in fasta-format to : " + output_fasta
