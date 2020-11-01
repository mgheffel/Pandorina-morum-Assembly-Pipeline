import os
import sys
import shutil
#converts fasta files to fastq files with quality line represented as highest possible quality
here=sys.argv[1]
if here=='-f':
    here='here'

#make output directory
if 'raw_fastqContigs' in os.listdir(here):
    shutil.rmtree(here+'/raw_fastqContigs')
os.mkdir(here+'/raw_fastqContigs')

for filename in os.listdir(here+'/raw_trimmedContigs'):
    if filename[-6:]!='.fasta':
        continue
    outFilename=here+'/raw_fastqContigs/'+filename[:-1]+'q'
    print(outFilename)
    filename=here+'/raw_trimmedContigs/'+filename
    print(filename)
    inf=open(filename,'r')
    outf=open(outFilename,'w')
    
    line=inf.readline()
    while line:
        if line[0]=='>':
            outf.write('@'+line[1:])
        else:
            outf.write(line)
            outf.write('+\n')
            qualLine='~'*(len(line)-1)
            outf.write(qualLine+'\n')
        line=inf.readline()
    inf.close()
    outf.close()
