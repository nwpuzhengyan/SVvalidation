import pysam
import ssw
import argparse
from pyfaidx import Fasta
import os

def through_bp(start_pos,end_pos,pos):
    if start_pos<pos-20 and end_pos>pos+20:
        return True
    else:
        return False

def write_read(chr,total_read,read_list,bp,line):
    for read_name in sorted(total_read):
        res.write('>'+str(read_name)+'\n')
        res.write(read_list[read_name][0] + '\n')

def read_bam(samfile,chr,pos1,pos2):
    read_num=0
    for read in samfile.fetch(chr, pos1-200, pos2+201,until_eof=True):
        if read.query_name in read_list:
            continue
        else:
            if read.is_supplementary == False and read.mapq>30:
                read_list[read.query_name] = [read.query_sequence]
                read_num+=1
                if read_num>50:
                    break

parser = argparse.ArgumentParser()
parser.add_argument("hg19_del_bam", help='bam file')
parser.add_argument("DEL", help='DEL file')
args = parser.parse_args()

MAPQ=0
coverage=40

samfile=pysam.AlignmentFile(args.hg19_del_bam, 'rb')
res= open('read.txt', 'w')
read_list = {}
f = open(args.DEL, "r")
line = f.readline()
while line:
    line = line.split()
    chr = str(line[0])
    m_pos=(int(line[1])+int(line[2]))/2
    read_bam(samfile,chr, int(line[1]), int(line[2]))
    line = f.readline()
total_read = list(set(sorted(read_list.keys())))
write_read(chr, total_read, read_list,m_pos,line)
f.close()
res.close()