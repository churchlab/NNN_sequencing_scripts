# script to find the NNN codon
# bases 31, 32, 33
# run it like this:
# python findNNNCodon_conditional.py inputfileName outputfileName

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

outfile = open(args[1]+".out",'w')
gca_outfile = open(args[1]+"_gca.out", "w")
gta_outfile = open(args[1]+"_gta.out", "w")
cta_outfile = open(args[1]+"_cta.out", "w")
other_outfile = open(args[1]+"_other.out", "w")

constSeq = 'GCGCGCCCAGTATGT' #tyrS
first_offset = 16
wobble_offset = first_offset + 3
check_codon_first = first_offset + 12
check_codon_last = check_codon_first + 3

# initialize the codon and aa dictionaries
codons = {
'AAA': 0, 'AAC': 0, 'AAG': 0, 'AAT': 0, 'ACA': 0, 'ACC': 0, 'ACG': 0, 'ACT': 0, 'AGA': 0, 'AGC': 0, 'AGG': 0, 'AGT': 0, 'ATA': 0, 'ATC': 0, 'ATG': 0, 'ATT': 0, 'CAA': 0, 'CAC': 0, 'CAG': 0, 'CAT': 0, 'CCA': 0, 'CCC': 0, 'CCG': 0, 'CCT': 0, 'CGA': 0, 'CGC': 0, 'CGG': 0, 'CGT': 0, 'CTA': 0, 'CTC': 0, 'CTG': 0, 'CTT': 0, 'GAA': 0, 'GAC': 0, 'GAG': 0, 'GAT': 0, 'GCA': 0, 'GCC': 0, 'GCG': 0, 'GCT': 0, 'GGA': 0, 'GGC': 0, 'GGG': 0, 'GGT': 0, 'GTA': 0, 'GTC': 0, 'GTG': 0, 'GTT': 0, 'TAA': 0, 'TAC': 0, 'TAG': 0, 'TAT': 0, 'TCA': 0, 'TCC': 0, 'TCG': 0, 'TCT': 0, 'TGA': 0, 'TGC': 0, 'TGG': 0, 'TGT': 0, 'TTA': 0, 'TTC': 0, 'TTG': 0, 'TTT': 0
}

aas = {
'*': 0, 'A': 0, 'C': 0, 'D': 0, 'E': 0, 'F': 0, 'G': 0, 'H': 0, 'I': 0, 'K': 0, 'L': 0, 'M': 0, 'N': 0, 'P': 0, 'Q': 0, 'R': 0, 'S': 0, 'T': 0, 'V': 0, 'W': 0, 'Y': 0
}

gca_table = {
'AAA': 0, 'AAC': 0, 'AAG': 0, 'AAT': 0, 'ACA': 0, 'ACC': 0, 'ACG': 0, 'ACT': 0, 'AGA': 0, 'AGC': 0, 'AGG': 0, 'AGT': 0, 'ATA': 0, 'ATC': 0, 'ATG': 0, 'ATT': 0, 'CAA': 0, 'CAC': 0, 'CAG': 0, 'CAT': 0, 'CCA': 0, 'CCC': 0, 'CCG': 0, 'CCT': 0, 'CGA': 0, 'CGC': 0, 'CGG': 0, 'CGT': 0, 'CTA': 0, 'CTC': 0, 'CTG': 0, 'CTT': 0, 'GAA': 0, 'GAC': 0, 'GAG': 0, 'GAT': 0, 'GCA': 0, 'GCC': 0, 'GCG': 0, 'GCT': 0, 'GGA': 0, 'GGC': 0, 'GGG': 0, 'GGT': 0, 'GTA': 0, 'GTC': 0, 'GTG': 0, 'GTT': 0, 'TAA': 0, 'TAC': 0, 'TAG': 0, 'TAT': 0, 'TCA': 0, 'TCC': 0, 'TCG': 0, 'TCT': 0, 'TGA': 0, 'TGC': 0, 'TGG': 0, 'TGT': 0, 'TTA': 0, 'TTC': 0, 'TTG': 0, 'TTT': 0
}

gta_table = {
'AAA': 0, 'AAC': 0, 'AAG': 0, 'AAT': 0, 'ACA': 0, 'ACC': 0, 'ACG': 0, 'ACT': 0, 'AGA': 0, 'AGC': 0, 'AGG': 0, 'AGT': 0, 'ATA': 0, 'ATC': 0, 'ATG': 0, 'ATT': 0, 'CAA': 0, 'CAC': 0, 'CAG': 0, 'CAT': 0, 'CCA': 0, 'CCC': 0, 'CCG': 0, 'CCT': 0, 'CGA': 0, 'CGC': 0, 'CGG': 0, 'CGT': 0, 'CTA': 0, 'CTC': 0, 'CTG': 0, 'CTT': 0, 'GAA': 0, 'GAC': 0, 'GAG': 0, 'GAT': 0, 'GCA': 0, 'GCC': 0, 'GCG': 0, 'GCT': 0, 'GGA': 0, 'GGC': 0, 'GGG': 0, 'GGT': 0, 'GTA': 0, 'GTC': 0, 'GTG': 0, 'GTT': 0, 'TAA': 0, 'TAC': 0, 'TAG': 0, 'TAT': 0, 'TCA': 0, 'TCC': 0, 'TCG': 0, 'TCT': 0, 'TGA': 0, 'TGC': 0, 'TGG': 0, 'TGT': 0, 'TTA': 0, 'TTC': 0, 'TTG': 0, 'TTT': 0
}

cta_table = {
'AAA': 0, 'AAC': 0, 'AAG': 0, 'AAT': 0, 'ACA': 0, 'ACC': 0, 'ACG': 0, 'ACT': 0, 'AGA': 0, 'AGC': 0, 'AGG': 0, 'AGT': 0, 'ATA': 0, 'ATC': 0, 'ATG': 0, 'ATT': 0, 'CAA': 0, 'CAC': 0, 'CAG': 0, 'CAT': 0, 'CCA': 0, 'CCC': 0, 'CCG': 0, 'CCT': 0, 'CGA': 0, 'CGC': 0, 'CGG': 0, 'CGT': 0, 'CTA': 0, 'CTC': 0, 'CTG': 0, 'CTT': 0, 'GAA': 0, 'GAC': 0, 'GAG': 0, 'GAT': 0, 'GCA': 0, 'GCC': 0, 'GCG': 0, 'GCT': 0, 'GGA': 0, 'GGC': 0, 'GGG': 0, 'GGT': 0, 'GTA': 0, 'GTC': 0, 'GTG': 0, 'GTT': 0, 'TAA': 0, 'TAC': 0, 'TAG': 0, 'TAT': 0, 'TCA': 0, 'TCC': 0, 'TCG': 0, 'TCT': 0, 'TGA': 0, 'TGC': 0, 'TGG': 0, 'TGT': 0, 'TTA': 0, 'TTC': 0, 'TTG': 0, 'TTT': 0
}

other_table = {
'AAA': 0, 'AAC': 0, 'AAG': 0, 'AAT': 0, 'ACA': 0, 'ACC': 0, 'ACG': 0, 'ACT': 0, 'AGA': 0, 'AGC': 0, 'AGG': 0, 'AGT': 0, 'ATA': 0, 'ATC': 0, 'ATG': 0, 'ATT': 0, 'CAA': 0, 'CAC': 0, 'CAG': 0, 'CAT': 0, 'CCA': 0, 'CCC': 0, 'CCG': 0, 'CCT': 0, 'CGA': 0, 'CGC': 0, 'CGG': 0, 'CGT': 0, 'CTA': 0, 'CTC': 0, 'CTG': 0, 'CTT': 0, 'GAA': 0, 'GAC': 0, 'GAG': 0, 'GAT': 0, 'GCA': 0, 'GCC': 0, 'GCG': 0, 'GCT': 0, 'GGA': 0, 'GGC': 0, 'GGG': 0, 'GGT': 0, 'GTA': 0, 'GTC': 0, 'GTG': 0, 'GTT': 0, 'TAA': 0, 'TAC': 0, 'TAG': 0, 'TAT': 0, 'TCA': 0, 'TCC': 0, 'TCG': 0, 'TCT': 0, 'TGA': 0, 'TGC': 0, 'TGG': 0, 'TGT': 0, 'TTA': 0, 'TTC': 0, 'TTG': 0, 'TTT': 0
}

gca_aas = {
'*': 0, 'A': 0, 'C': 0, 'D': 0, 'E': 0, 'F': 0, 'G': 0, 'H': 0, 'I': 0, 'K': 0, 'L': 0, 'M': 0, 'N': 0, 'P': 0, 'Q': 0, 'R': 0, 'S': 0, 'T': 0, 'V': 0, 'W': 0, 'Y': 0
}

gta_aas = {
'*': 0, 'A': 0, 'C': 0, 'D': 0, 'E': 0, 'F': 0, 'G': 0, 'H': 0, 'I': 0, 'K': 0, 'L': 0, 'M': 0, 'N': 0, 'P': 0, 'Q': 0, 'R': 0, 'S': 0, 'T': 0, 'V': 0, 'W': 0, 'Y': 0
}

cta_aas = {
'*': 0, 'A': 0, 'C': 0, 'D': 0, 'E': 0, 'F': 0, 'G': 0, 'H': 0, 'I': 0, 'K': 0, 'L': 0, 'M': 0, 'N': 0, 'P': 0, 'Q': 0, 'R': 0, 'S': 0, 'T': 0, 'V': 0, 'W': 0, 'Y': 0
}

other_aas = {
'*': 0, 'A': 0, 'C': 0, 'D': 0, 'E': 0, 'F': 0, 'G': 0, 'H': 0, 'I': 0, 'K': 0, 'L': 0, 'M': 0, 'N': 0, 'P': 0, 'Q': 0, 'R': 0, 'S': 0, 'T': 0, 'V': 0, 'W': 0, 'Y': 0
}


#codons = {}
#aas = {}

for i in range(len(readsFile)):
    if i == 0 or i%4 == 0:
        rSeq = readsFile[i+1]
        #print(rSeq)
        index = rSeq.find(constSeq)
        if len(rSeq) <= index+check_codon_last: # read too short
            print "%s is too short" %rSeq
            continue
	#print index
	#print rSeq
        codon = rSeq[index+first_offset:index+wobble_offset]
        check_codon = rSeq[ index+check_codon_first:index+check_codon_last ]
	if check_codon == "GCA":
	    if codon not in gca_table:
                gca_table[ codon ] = 1
            else:
                gca_table[ codon ] = gca_table[ codon ] + 1
            coding_dna = Seq(codon,IUPAC.unambiguous_dna)
	    translation = str(coding_dna.translate())
	    if translation not in gca_aas:
	        gca_aas[translation] = 1
	    else:
	        gca_aas[translation] = gca_aas[translation] + 1

	elif check_codon == "GTA":
	    if codon not in gta_table:
                gta_table[ codon ] = 1
            else:
                gta_table[ codon ] = gta_table[ codon ] + 1
            coding_dna = Seq(codon,IUPAC.unambiguous_dna)
	    translation = str(coding_dna.translate())
	    if translation not in gta_aas:
	        gta_aas[translation] = 1
	    else:
	        gta_aas[translation] = gta_aas[translation] + 1

	elif check_codon == "CTA":
	    if codon not in cta_table:
                cta_table[ codon ] = 1
            else:
                cta_table[ codon ] = cta_table[ codon ] + 1
            coding_dna = Seq(codon,IUPAC.unambiguous_dna)
	    translation = str(coding_dna.translate())
	    if translation not in cta_aas:
	        cta_aas[translation] = 1
	    else:
	        cta_aas[translation] = cta_aas[translation] + 1
            coding_dna = Seq(codon,IUPAC.unambiguous_dna)
	else:
	    if codon not in other_table:
                other_table[ codon ] = 1
            else:
                other_table[ codon ] = other_table[ codon ] + 1
            coding_dna = Seq(codon,IUPAC.unambiguous_dna)
	    translation = str(coding_dna.translate())
	    if translation not in other_aas:
	        other_aas[translation] = 1
	    else:
	        other_aas[translation] = other_aas[translation] + 1


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

## output all
##
##for key in codons.keys():
##    outfile.write("%s\t%d\n" %(key,codons[key]))
##    print "%s\t%d" %(key,codons[key])
##
##for aa in aas.keys():
##    outfile.write("%s\t%d\n" %(aa,aas[aa]))
##    print "%s\t%d" %(aa,aas[aa])
##outfile.close()

## output gca conditional only

#for key in gca_table.keys():
for key in sorted(gca_table.iterkeys()):
    gca_outfile.write("%s\t%d\n" %(key,gca_table[key]))
    print "%s\t%d" %(key,gca_table[key])

#for aa in gca_aas.keys():
for aa in sorted(gca_aas.iterkeys()):
    gca_outfile.write("%s\t%d\n" %(aa,gca_aas[aa]))
    print "%s\t%d" %(aa,gca_aas[aa])
gca_outfile.close()

## output gta conditional only

#for key in gta_table.keys():
for key in sorted(gta_table.iterkeys()):
    gta_outfile.write("%s\t%d\n" %(key,gta_table[key]))
    print "%s\t%d" %(key,gta_table[key])

#for aa in gta_aas.keys():
for aa in sorted(gta_aas.iterkeys()):
    gta_outfile.write("%s\t%d\n" %(aa,gta_aas[aa]))
    print "%s\t%d" %(aa,gta_aas[aa])
gta_outfile.close()

## output cta conditional only

#for key in cta_table.keys():
for key in sorted(cta_table.iterkeys()):
    cta_outfile.write("%s\t%d\n" %(key,cta_table[key]))
    print "%s\t%d" %(key,cta_table[key])

#for aa in cta_aas.keys():
for aa in sorted(cta_aas.iterkeys()):
    cta_outfile.write("%s\t%d\n" %(aa,cta_aas[aa]))
    print "%s\t%d" %(aa,cta_aas[aa])
cta_outfile.close()

## output other conditional only

#for key in other_table.keys():
for key in sorted(other_table.iterkeys()):
    other_outfile.write("%s\t%d\n" %(key,other_table[key]))
    print "%s\t%d" %(key,other_table[key])

#for aa in other_aas.keys():
for aa in sorted(other_aas.iterkeys()):
    other_outfile.write("%s\t%d\n" %(aa,other_aas[aa]))
    print "%s\t%d" %(aa,other_aas[aa])
other_outfile.close()

