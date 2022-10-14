import argparse

parser = argparse.ArgumentParser()
parser.add_argument("txt", help='inversion position txt')
args = parser.parse_args()

f = open(args.txt)
SV = open(r'SV.txt','w')
chrname=''

line=f.readline()

while line:
    line = line.split()
    if line[0]!=chrname:
        DIS=0
    INS_len=len(line[3])
    pos1=int(line[1])+DIS
    pos2 = int(line[1])+INS_len-(int(line[2])-int(line[1])) + DIS
    SV.write(str(line[0]) + ' '+str(line[1])+' '+ str(pos1) + ' ' + str(pos2)+' '+str(INS_len) + '\n')
    DIS=DIS+INS_len
    chrname=str(line[0])
    line = f.readline()
f.close()
SV.close()