import pysam
import ssw
import argparse
from pyfaidx import Fasta
import os

def through_bp(start_pos,end_pos,pos1,pos2):
    if start_pos<pos1-20 and end_pos>pos1+20:
        return pos1
    elif start_pos<pos2-20 and end_pos>pos2+20:
        return pos2
    else:
        return 0

def soft_clip(chr,start_pos,cigar,read_str):
    S_clipping = [0, 0]
    i = 0
    sign = 0
    cigar_str = ''
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
        if (cigar[i] != 'H'):
            cigar_str = cigar_str + cigar[i] * int(num)
        i += 1
    res=[S_clipping[0],S_clipping[1],cigar_str]
    return res


def divide_read(chr,start_pos,cigar,read_str,bp,genes,newpos):
    S_clipping = [newpos[0], newpos[1]]
    cigar_str=newpos[2]
    str1=''
    str2=''
    bp_pos=bp-start_pos+newpos[0]
    ref = genes[chr][start_pos-S_clipping[0]:start_pos + len(cigar_str)+S_clipping[1]]
    read = read_str[0:len(read_str)]
    pos = 0
    pos1 = 0
    pos2 = 0
    while pos < len(cigar_str) and pos1 < len(read) and pos2 < len(ref):
        if cigar_str[pos] == 'M' or cigar_str[pos] == 'S':
            if pos2 < bp_pos and pos2 > (bp_pos - 5000):
                str1 += read[pos1]
            elif pos2 >= bp_pos and pos2 < (bp_pos + 5000):
                str2 += read[pos1]
            pos1 += 1
            pos2 += 1
        elif cigar_str[pos] == 'D':
            pos2 += 1
        else:
            if pos2 < bp_pos and pos2 > (bp_pos - 5000):
                str1 += read[pos1]
            elif pos2 >= bp_pos and pos2 < (bp_pos + 5000):
                str2 += read[pos1]
            pos1 += 1
        pos += 1
    return[str1,str2]

def write_read(chr,total_read,read_list,bp1,bp2,line,genes):
    for read_name in sorted(total_read):
        newpos=soft_clip(chr,read_list[read_name][2],read_list[read_name][5],read_list[read_name][6])
        newstart=read_list[read_name][2]-newpos[0]
        newend=read_list[read_name][3]+newpos[1]
        bp=through_bp(newstart,newend,bp1,bp2)
        if bp>0:
            read_seq=divide_read(chr,read_list[read_name][2],read_list[read_name][5],read_list[read_name][6],bp,genes,newpos)
            res.write('>'+str(line[0])+':'+str(line[1])+'-'+str(line[1])+':'+str(read_name)+':l'+'\n')
            res.write(read_seq[0]+'\n')
            res.write('>' + str(line[0]) + ':' + str(line[1]) + '-' + str(line[1]) + ':' + str(read_name) + ':r'+'\n')
            res.write(read_seq[1] + '\n')

def read_bam(samfile,chr,pos1,pos2):
    read_list = {}
    for read in samfile.fetch(chr, pos1-100, pos2+100):
        if read.query_name in read_list:
            if read.mapq > read_list[read.query_name][0] and read.is_secondary==False and\
                read.reference_start < pos1 - 500 and read.reference_end > pos2 + 500:
                read_list[read.query_name] = [read.mapq, read.tags[2][1], read.reference_start,read.reference_end,\
                                              read.is_secondary,read.cigarstring,read.query_sequence,read.query_length]
        else:
            if read.is_secondary==False and read.reference_start<pos1-500 and read.reference_end>pos2+500:
                read_list[read.query_name] = [read.mapq, read.tags[2][1], read.reference_start,read.reference_end,\
                                          read.is_secondary,read.cigarstring,read.query_sequence,read.query_length]
    return read_list

parser = argparse.ArgumentParser()
parser.add_argument("hg19_bam", help='bam file')
parser.add_argument("fasta1", help='fasta file')
parser.add_argument("hg19_INS_bam", help='bam file')
parser.add_argument("fasta2", help='fasta file')
parser.add_argument("INS", help='INS file')
args = parser.parse_args()

MAPQ=30
coverage=40

genes1 = Fasta(args.fasta1)
genes2 = Fasta(args.fasta2)
samfile1=pysam.AlignmentFile(args.hg19_bam, 'rb')
samfile2=pysam.AlignmentFile(args.hg19_INS_bam, 'rb')
res= open('divide_read.txt', 'w')
f = open(args.INS, "r")
line = f.readline()
while line:
    line = line.split()
    print(line)
    chr = str(line[0])
    read_list = read_bam(samfile1,chr, int(line[1]), int(line[1]))
    total_read = list(set(sorted(read_list.keys())))
    INS_read_list = read_bam(samfile2, chr, int(line[2]), int(line[3]))
    for k in INS_read_list.keys():
        if k in total_read:
            if  INS_read_list[k][1]<read_list[k][1]:
                del INS_read_list[k]
            else:
                del read_list[k]
    total_read = list(set(sorted(read_list.keys())))
    write_read(chr, total_read, read_list,int(line[1]), int(line[1]), line, genes1)
    total_read = list(set(sorted(INS_read_list.keys())))
    m_pos = (int(line[3]) + int(line[2])) / 2
    write_read(chr, total_read,INS_read_list,m_pos, m_pos, line, genes2)
    line = f.readline()
f.close()
res.close()