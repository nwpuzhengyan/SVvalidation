import os
import pysam
import argparse
from pyfaidx import Fasta
import time

parser = argparse.ArgumentParser()
parser.add_argument("test", help='test file')
parser.add_argument("reference", help='benchmark file')
parser.add_argument("bam", help='bam file')
args = parser.parse_args()

os.system('python2 ../INS_pos.py '+args.test)
os.system('python2 ../create_INS.py '+args.test+' '+args.reference)
os.system('python2 ../extract_primary.py '+args.bam+' SV.txt')
os.system('minimap2 -ax map-pb INS_hg19.fa read.txt > hg19_INS.sam')
os.system('samtools view -bS hg19_INS.sam > hg19_INS.bam')
os.system('samtools sort hg19_INS.bam > hg19_INS.sort.bam')
os.system('samtools index hg19_INS.sort.bam')

os.system('python2 ../divide_read_INS.py '+args.bam+' '+args.reference+' hg19_INS.sort.bam INS_hg19.fa SV.txt')
os.system('minimap2 -ax map-pb INS_hg19.fa divide_read.txt > M.sam')
os.system('samtools view -bS M.sam > M.bam')
os.system('samtools sort M.bam > M.sort.bam')
os.system('samtools index M.sort.bam')
os.system('python2 ../get_breakpoint_INS.py M.sort.bam divide_read.txt SV.txt INS_hg19.fa')
