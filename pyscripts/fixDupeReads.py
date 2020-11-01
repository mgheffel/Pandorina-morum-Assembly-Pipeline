#prevents duplicate read IDs by appending indexes to reads with equal IDs
import os
import sys
import shutil


#sys.argv[1] should be input fastq file
inFile=sys.argv[1]
outFile=inFile[:-6]+'-noDupe.fastq'

readID_dict={}
count=0
outf=open(outFile,'w')
with open(inFile,'r') as f:
	line=f.readline()
	while line:
		readID=line[1:-8]
		while readID in readID_dict:
			try:
				splitRead=readID.split(' ')
				readID=splitRead[0]+'- '+splitRead[1]
			except:
				print(line)
				readID+='-'
			#readID+='-'
		readID_dict[readID]=True
		outf.write('@'+readID+line[-8:])
		line=f.readline()
		outf.write(line)
		line=f.readline()
		outf.write(line)
		line=f.readline()
		outf.write(line)
		line=f.readline()
outf.close()

for key in readID_dict.keys():
	print(key)
