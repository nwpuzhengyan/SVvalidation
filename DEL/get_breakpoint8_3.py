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

def soft_len(cigar):
    i=0
    S_clipping = [0, 0]
    sign = 0
    while i < len(cigar):
        num = ''
        while cigar[i].isdigit():
            num = num + cigar[i]
            i += 1
        if sign == 0 and cigar[i] == 'S':
            S_clipping[0] = int(num)
        if sign == 1 and cigar[i] == 'S':
            S_clipping[1] = int(num)
        sign = 1
        i += 1
    return S_clipping

def calculate_read_num(chr,total_read, read_list,pos1,pos2):
    single_read = 0
    useless_read = 0
    support_read = 0
    reject_read = 0
    dis_len=(pos2-pos1)*0.2
    for read_name in sorted(total_read):
        if len(read_list[read_name])!=2:
            direct= list(read_list[read_name].keys())
            if direct[0]=='l':
                single_read+=1
            else:
                single_read += 1
        else:
            end_pos=read_list[read_name]['l'][3]
            start_pos=read_list[read_name]['r'][2]
            gap_len=start_pos-end_pos
            DEL_len=abs(pos2-pos1)
            if abs(DEL_len-gap_len)<dis_len:
                support_read+=1
            else:
                reject_read+=1
    if (support_read+reject_read)==0:
        rate=0
    else:
        rate=float(support_read)/float(support_read+reject_read)
    useful_read=support_read+reject_read
    if support_read>=coverage1:
        if rate>0.7:
            SV_type = 'homo'
        elif rate>=0.3 and rate<=0.7:
            SV_type = 'heter'
        elif rate<0.3:
            SV_type = 'false'
    else:
        SV_type='false'
    res = str(useful_read) + ' ' + str(support_read) + ' ' + str(reject_read)+ ' ' + str(single_read) + ' ' + str(useless_read) + ' ' + SV_type
    return res

def read_bam(samfile,chr,pos1,pos2):
    read_list = {}
    for read in samfile.fetch(chr, pos1-200, pos2+200):
        read_name = read.query_name
        read_name = read_name.split(':')
        DEL_pos=str(chr)+':'+str(pos1)+'-'+str(pos2)
        DEL=str(read_name[0])+':'+str(read_name[1])
        if DEL_pos!=DEL:
            continue
        if read_name[2] in read_list:
            if read_name[3] in read_list[read_name[2]]:
                if read.tags[2][1]>read_list[read_name[2]][read_name[3]][1]:
                    read_list[read_name[2]][read_name[3]] = [read.mapq, read.tags[2][1], read.reference_start, \
                                                             read.reference_end, read.cigarstring, str(read_name[3])]
            else:
                read_list[read_name[2]][read_name[3]] = [read.mapq, read.tags[2][1], read.reference_start, \
                                                     read.reference_end, read.cigarstring, str(read_name[3])]
        else:
            read_list[read_name[2]] = {}
            read_list[read_name[2]][read_name[3]] = [read.mapq, read.tags[2][1], read.reference_start, \
                                                     read.reference_end, read.cigarstring, str(read_name[3])]
    return read_list

parser = argparse.ArgumentParser()
parser.add_argument("hg19_bam", help='bam file')
parser.add_argument("read", help='read file')
parser.add_argument("DEL", help='DEL file')
parser.add_argument("fasta1", help='fasta file')
args = parser.parse_args()
aligner = ssw.Aligner()
#aligner.gap_extend=2
#aligner.gap_open=3


MAPQ=30
coverage1=3

genes1 = Fasta(args.fasta1)
samfile=pysam.AlignmentFile(args.hg19_bam, 'rb')
res= open('SV_type.txt', 'w')
single_rate=open('single.txt', 'w')
f = open(args.DEL, "r")
line = f.readline()
while line:
    line = line.split()
    print(line)
    res.write(str(line[0]) + ':' + str(line[1]) + '-' + str(line[2]) + ' ')
    chr = str(line[0])
    read_list = read_bam(samfile,chr,int(line[1]), int(line[2]))
    total_read = list(set(sorted(read_list.keys())))
    read_num = calculate_read_num(chr, total_read, read_list,int(line[1]), int(line[2]))
    res.write(read_num + '\n')
    line = f.readline()
f.close()
res.close()
single_rate.close()