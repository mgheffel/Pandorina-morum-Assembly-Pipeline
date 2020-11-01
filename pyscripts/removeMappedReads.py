#removes reads from forward and reverse files that match secies declared in blas step of pipeline (Bacteria superkingdom)
import sys
import os


samFile=sys.argv[1]
r1File=sys.argv[2]
r2File=sys.argv[3]
#fReadFile=sys.argv[2]
#rReadFile=sys.argv[3]


#nonBacReadsDict={}
BacReadsDict={}
f=open(samFile,'r')
line=f.readline()
while line:
	if line[0]=='@':
		line=f.readline()
		continue
	splitLine=line.split('\t')
	mappedRef=splitLine[2]
	#if mappedRef=='*':
	#	rName=splitLine[0]
	#	nonBacReadsDict[rName]=1
	if mappedRef!='*':
		rName=splitLine[0]
		BacReadsDict[rName]=1
	line=f.readline()

#take reads from forward read file
f=open(r1File,'r')
outf=open(r1File[:-6]+'_noBac.fastq','w')
line=f.readline()
count=0
while line:
	if count==0 and line[0]=='@':
		nameline=line
		rName=line.split(' ')[0][1:]
		count=1
	elif count==1:
		seq=line
		count=2
	elif count==2:
		count=3
	else:
		qual=line
		count=0
		#if rName in nonBacReadsDict:
		#	outf.write(nameline)
		#	outf.write(seq)
		#	outf.write('+\n')
		#	outf.write(qual)
		#else:
		#	print('Removed Read: '+rName)
		#	print(nameline)
		if rName in BacReadsDict:
			print('Removed Read: '+rName)
			print(nameline)
		else:
			outf.write(nameline)
			outf.write(seq)
			outf.write('+\n')
			outf.write(qual)
	line=f.readline()
f.close()
outf.close()

#take reads from reverse reads file
f=open(r2File,'r')
outf=open(r2File[:-6]+'_noBac.fastq','w')
line=f.readline()
count=0
while line:
        if count==0 and line[0]=='@':
                nameline=line
                rName=line.split(' ')[0][1:]
                count=1
        elif count==1:
                seq=line
                count=2
        elif count==2:
                count=3
        else:
                qual=line
                count=0
                if rName in BacReadsDict:
                        print('Removed Read: '+rName)
                        print(nameline)
                else:
                        outf.write(nameline)
                        outf.write(seq)
                        outf.write('+\n')
                        outf.write(qual)
        line=f.readline()
f.close()
outf.close()

print("Removed "+str(len(BacReadsDict.keys()))+" bacterial reads")
