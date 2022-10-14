import argparse

parser = argparse.ArgumentParser()
parser.add_argument("txt", help='inversion position txt')
parser.add_argument("output", help='output name')
args = parser.parse_args()

f = open(args.txt)
SV = open(args.output,'w')
chrname=''

line=f.readline()

while line:
    line = line.split()
    if line[0]!=chrname:
        DIS=0
    pos1=int(line[1])+DIS
    pos2 = int(line[1]) + DIS
    SV.write(str(line[0]) + ' ' + str(line[1])+' '+str(line[2]) + ' ' + str(pos2) +'\n')
    DIS=DIS+int(line[1])-int(line[2])
    chrname=str(line[0])
    line = f.readline()

f.close()
SV.close()