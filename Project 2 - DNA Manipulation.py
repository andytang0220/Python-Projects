    # -*- coding: utf-8 -*-
codons = {"UUU":"Phe", "UUC":"Phe", "UUA":"Leu", "UUG":"Leu",
          "UCU":"Ser", "UCC":"Ser", "UCA":"Ser", "UCG":"Ser",
          "UAU":"Tyr", "UAC":"Tyr", "UAA":"STOP", "UAG":"STOP",
          "UGU":"Cys", "UGC":"Cys", "UGA":"STOP", "UGG":"Trp",
          "CUU":"Leu", "CUC":"Leu", "CUA":"Leu", "CUG":"Leu",
          "CCU":"Pro", "CCC":"Pro", "CCA":"Pro", "CCG":"Pro",
          "CAU":"His", "CAC":"His", "CAA":"Gln", "CAG":"Gln",
          "CGU":"Arg", "CGC":"Arg", "CGA":"Arg", "CGG":"Arg",
          "AUU":"Ile", "AUC":"Ile", "AUA":"Ile", "AUG":"Met",
          "ACU":"Thr", "ACC":"Thr", "ACA":"Thr", "ACG":"Thr",
          "AAU":"Asn", "AAC":"Asn", "AAA":"Lys", "AAG":"Lys",
          "AGU":"Ser", "AGC":"Ser", "AGA":"Arg", "AGG":"Arg",
          "GUU":"Val", "GUC":"Val", "GUA":"Val", "GUG":"Val",
          "GCU":"Ala", "GCC":"Ala", "GCA":"Ala", "GCG":"Ala",
          "GAU":"Asp", "GAC":"Asp", "GAA":"Glu", "GAG":"Glu",
          "GGU":"Gly", "GGC":"Gly", "GGA":"Gly", "GGG":"Gly"}
        
def main():
    """
    dna = readFile('cftr1.txt')
    comp_dna = replication(dna)
    writeFile('replication.txt',comp_dna)
    """
    #read each dna file and determine if it is normal or abnormal
    dna1 = cftr(readFile('cftr1.txt'))
    dna2 = cftr(readFile('cftr2.txt'))
    dna3 = cftr(readFile('cftr3.txt'))
    dna4 = cftr(readFile('cftr4.txt'))
    dna5 = cftr(readFile('cftr5.txt'))
    dna = [dna1, dna2, dna3, dna4, dna5]
    for x in range(5):
        if dna[x] == True:
            print('Sequence', (x+1), 'codes for normal CFTR')
        elif dna[x] == False:
            print('Sequence', (x+1), 'codes for abnormal CFTR')
    
def readFile(fileName):
    """
    input: fileName to read
    output: text in file as string
    """
    with open(fileName,'r') as textFile:
        text = textFile.read()
        
    return text.strip() 
    
def writeFile(fileName,text):
    """
    input: fileName to write
           text to write as string
    """
    with open(fileName,'w') as textFile:
        textFile.write(text.strip())

def replication(dna):
    comp_dna = '' #initialize the complementary dna string
    for x in range(len(dna)): #build up the complementary string based on the original
        if dna[x] == 'A':
            comp_dna += 'T'
        elif dna[x] == 'T':
            comp_dna += 'A'
        elif dna[x] == 'G':
            comp_dna += 'C'
        elif dna[x] == 'C':
            comp_dna += 'G'
    return comp_dna

def transcription(dna):
    mRNA_r = [] #set up the mRNA list that will later be joined into a string
    for x in dna: #replace each T with a U
        if x == 'T':
            mRNA_r.append('U')
        else:
            mRNA_r.append(x)
    mRNA = '' #initialize mRNA string
    return mRNA.join(mRNA_r)
    
    
def translation(mrna):
    mrna_list = []
    amino_acids = []
    for x in mrna: #add the individual bases into the list
        mrna_list.append(x)
    for x in range(0, len(mrna), 3):
        if (len(mrna)-3) < x: #if the next group of bases is less than 3, no codon will be generated for that group
            break
        key = '{}{}{}'.format(mrna[x], mrna[x+1], mrna[x+2]) #dict key for amino acids
        amino_acids += codons[key], ' '
    residue_seq = ''
    return residue_seq.join(amino_acids)
    

def restriction(dna, seq):    
    frags = dna.split(seq) #generate fragments of the dna strand based on where the desired seq is
    for x in range(len(frags)):
        if x == 0: #first fragment will only have the first half of seq on the end
            frags[x] += seq[0:len(seq) // 2]
        elif x == len(frags) - 1: #last fragment will only have the second half of seq at the beginning
            frags[x] = seq[len(seq) // 2: len(seq)] + frags[x]
        else: #all other fragments will have the second half at the front and the first half on the end
            frags[x] = seq[len(seq) // 2: len(seq)] + frags[x] + seq[0:len(seq) // 2]
    return frags
        
            
def isolate(dna, start, end):
    if start and end in dna: #ensure that start and end are actually in the dna sequence
        frag_set1 = dna.split(start)
        for x in frag_set1:
           if end in x: #only continue to split the fragments that have the end sequence
               frag1 = x
        frag_set2 = frag1.split(end)
        isolate = frag_set2[0] #only interested in the first split sequence
        return isolate
    else: #if neither start nor end in dna sequence, return False
       return False
   
def cftr(dna):
    #establish start and end sequences
    start = 'ATTAAAGAAAATATC'
    end = 'GGTGTTTCCTATGAT'
    seq = isolate(dna, start, end)
    if seq == 'ATCTTT': #normal CFTR sequence
        return True
    else:
        return False
 
if __name__ == "__main__": main()
