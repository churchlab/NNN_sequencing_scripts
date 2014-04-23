# script to find the NNN codon
# bases 31, 32, 33
# run it like this:
# python findNNNCodon.py inputfileName outputfileName

import optparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import gzip

from sys import exit

parser = optparse.OptionParser()
(options,args) = parser.parse_args()

#print args[0][-3:]

if (args[0][-3:] == ".gz"):
    readsFile = gzip.open(args[0]).readlines()
else: readsFile = open(args[0]).readlines()

outfile = open(args[1],'w')
constSeq = 'GTGATGCTTCCGCTCAGG' #pgk
first_offset = 32
wobble_offset = first_offset + 3

# initialize the codon and aa dictionaries
codons = {
'AAA': 0, 'AAC': 0, 'AAG': 0, 'AAT': 0, 'ACA': 0, 'ACC': 0, 'ACG': 0, 'ACT': 0, 'AGA': 0, 'AGC': 0, 'AGG': 0, 'AGT': 0, 'ATA': 0, 'ATC': 0, 'ATG': 0, 'ATT': 0, 'CAA': 0, 'CAC': 0, 'CAG': 0, 'CAT': 0, 'CCA': 0, 'CCC': 0, 'CCG': 0, 'CCT': 0, 'CGA': 0, 'CGC': 0, 'CGG': 0, 'CGT': 0, 'CTA': 0, 'CTC': 0, 'CTG': 0, 'CTT': 0, 'GAA': 0, 'GAC': 0, 'GAG': 0, 'GAT': 0, 'GCA': 0, 'GCC': 0, 'GCG': 0, 'GCT': 0, 'GGA': 0, 'GGC': 0, 'GGG': 0, 'GGT': 0, 'GTA': 0, 'GTC': 0, 'GTG': 0, 'GTT': 0, 'TAA': 0, 'TAC': 0, 'TAG': 0, 'TAT': 0, 'TCA': 0, 'TCC': 0, 'TCG': 0, 'TCT': 0, 'TGA': 0, 'TGC': 0, 'TGG': 0, 'TGT': 0, 'TTA': 0, 'TTC': 0, 'TTG': 0, 'TTT': 0
}

aas = {
'*': 0, 'A': 0, 'C': 0, 'D': 0, 'E': 0, 'F': 0, 'G': 0, 'H': 0, 'I': 0, 'K': 0, 'L': 0, 'M': 0, 'N': 0, 'P': 0, 'Q': 0, 'R': 0, 'S': 0, 'T': 0, 'V': 0, 'W': 0, 'Y': 0
}

for i in range(len(readsFile)):
    if i == 0 or i%4 == 0:
        rSeq = readsFile[i+1]
        #print(rSeq)
        index = rSeq.find(constSeq)
        if len(rSeq) <= index+wobble_offset: # read too short
            print "%s is too short" %rSeq
            continue
	#print index
	#print rSeq
        codon = rSeq[index+first_offset:index+wobble_offset]
        #outfile.write(codon + '\n')
        if codon not in codons:
            codons[codon] = 1
        else:
            codons[codon] = codons[codon] + 1
        coding_dna = Seq(codon,IUPAC.unambiguous_dna)
        translation = str(coding_dna.translate())
        if translation not in aas:
            aas[translation] = 1
        else:
            aas[translation] = aas[translation] + 1

#for key in codons.keys():
for key in sorted(codons.iterkeys()):
    outfile.write("%s\t%d\n" %(key,codons[key]))
    print "%s\t%d" %(key,codons[key])

#for aa in aas.keys():
for aa in sorted(aas.iterkeys()):
    outfile.write("%s\t%d\n" %(aa,aas[aa]))
    print "%s\t%d" %(aa,aas[aa])

outfile.close()
