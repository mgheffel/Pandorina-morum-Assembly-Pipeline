#reads files in from a directory called scaffolds_blast within the running directory (here) and downloads reference
#--sequences whose superkingdom match the regex r'[Bb]acteri[au]' into a file called 
#-- '/bacterialBlastSequences/bacRefSeqs.fasta'
import os
import sys
from Bio import Entrez
from Bio import SeqIO
import time
import re
Entrez.email = "mgheffel@ksu.edu"
#sys.argv[1]=dataset directory
here=sys.argv[1]
#here=os.getcwd()
#os.mkdir(here+'/bacterialBlastSequences')
outputFile=here+'/bacterialBlastSequences/bacRefSeqs.fasta'
regex=re.compile(r'[Bb]acteri[au]')
bacterialRefs=[]
for filename in os.listdir(here+'/scaffolds_blast'):
    print(filename)
    with open('scaffolds_blast/'+filename,'r') as f:
        line=f.readline()
        while line:
            splitline=line.split('\t')
            superkingdom=splitline[7]
            print(superkingdom)
            print(regex.search(superkingdom))
            if not regex.search(superkingdom):
                line=f.readline()
                continue
            refID=splitline[2]
            if refID not in bacterialRefs:
                bacterialRefs.append(refID)
            #scaffoldName=splitline[0]
            #eScore=splitline[1]
            #speciesName=splitline[7]
            line=f.readline()

print(len(bacterialRefs))

outf=open(outputFile,'w')
for refID in bacterialRefs:
    print(refID)
    flag=True
    while flag:
        try:
            handle = Entrez.efetch(db="nucleotide", id=refID, rettype="gb", retmode="text")
            record = SeqIO.read(handle,"genbank")
            handle.close()
            seq=record.seq
            try:
                name=record.description
            except:
                name='name-not-found'
                print('nnn')
            outf.write('>'+refID+'\n')
            outf.write(str(seq)+'\n')
            flag=False
        except:
            time.sleep(10)
print('sequences written to: '+outputFile)
outf.close()
