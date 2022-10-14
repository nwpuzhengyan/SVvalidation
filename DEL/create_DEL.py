import argparse
import os
from pyfaidx import Fasta

parser = argparse.ArgumentParser()
parser.add_argument("txt", help='inversion position txt')
parser.add_argument("fasta", help='hg19 fasta file')
args = parser.parse_args()

f = open(args.txt)
chr_pos={}
chr_num=0
for line in f.readlines():
    v = line.split()
    if int(v[2])-int(v[1])!=0:
       chr_pos.setdefault(chr_num, []).append(v[0])
       chr_pos.setdefault(chr_num, []).append(v[1])
       chr_pos.setdefault(chr_num, []).append(v[2])
       chr_num=chr_num+1
f.close()
genes = Fasta(args.fasta)
print (genes.keys())


Write_inv = open(r'DEL_hg19.fa','w')
for chr in genes.keys():
    Write_inv.write('>' + chr + '\n')
    start=0
    chr_use=0
    key_num=0
    for key in chr_pos.keys():
        key_num=key_num+1
        if chr==chr_pos[key][0]:
            print (chr+' '+str(chr_pos[key][1])+' '+str(chr_pos[key][2])+'\n')
            chr_use=1
            Write_inv.write(str(genes[chr][start:int(chr_pos[key][1])]))
            start=int(chr_pos[key][2])
            if (key_num==len(chr_pos.keys())):
                Write_inv.write(str(genes[chr])[start:len(genes[chr])] + '\n')
                break
        elif chr!=chr_pos[key][0] and chr_use==1:
            Write_inv.write(str(genes[chr])[start:len(genes[chr])] + '\n')
            break
    if chr_use==0:
        Write_inv.write(str(genes[chr]) + '\n')
