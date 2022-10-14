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

os.system('python2 ../DEL_pos.py '+args.test+' SV.txt')
os.system('python2 ../create_DEL.py SV.txt '+args.reference)
os.system('python2 ../extract_primary.py '+args.bam+' SV.txt')
os.system('minimap2 -ax map-pb DEL_hg19.fa read.txt > hg19_del.sam')
os.system('samtools view -bS hg19_del.sam > hg19_del.bam')
os.system('samtools sort hg19_del.bam > hg19_del.sort.bam')
os.system('samtools index hg19_del.sort.bam')

os.system('python2 ../divide_read8_2.py '+args.bam+' '+args.reference+' hg19_del.sort.bam DEL_hg19.fa SV.txt')
os.system('minimap2 -ax map-pb '+args.reference+' divide_read.txt > M.sam')
os.system('samtools view -bS M.sam > M.bam')
os.system('samtools sort M.bam > M.sort.bam')
os.system('samtools index M.sort.bam')
os.system('python2 ../get_breakpoint8_3.py M.sort.bam divide_read.txt SV.txt '+args.reference)
