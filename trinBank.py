'''Program for processing weedBank files and allowing for calling individual 
attributes of the file.
'''
import re

codons = {'atg': 'M', 'ttt': 'F', 'ttc': 'F', 'tta': 'L', 'ttg': 'L', 'ctt': 'L',
'ctc': 'L', 'cta': 'L', 'ctg': 'L','att': 'I', 'atc': 'I', 'ata': 'I', 'gtt': 'V',
'gtc': 'V', 'gta': 'V', 'gtg': 'V', 'tct': 'S', 'tcc': 'S', 'tca': 'S', 'tcg': 'S',
'cct': 'P', 'ccc': 'P', 'cca': 'P', 'ccg': 'P', 'act': 'T', 'acc': 'T', 'aca': 'T',
'acg': 'T', 'gct': 'A', 'gcc': 'A', 'gca': 'A', 'gcg': 'A', 'tat': 'Y', 'tac': 'Y', 
'taa': '*', 'tag': '*', 'cat': 'H', 'cac': 'H', 'caa': 'Q', 'cag': 'Q', 'aat': 'N',
'aac': 'N', 'aaa': 'K', 'aag': 'K', 'gat': 'D', 'gac': 'D', 'gaa': 'E', 'gag': 'E',
'tgt': 'C', 'tgc': 'C', 'tga': '*', 'tgg': 'W', 'cgt': 'R', 'cgc': 'R', 'cga': 'R',
'cgg': 'R', 'agt': 'S', 'agc': 'S', 'aga': 'R', 'agg': 'R', 'ggt': 'G', 'ggc': 'G',
'gga': 'G', 'ggg': 'G'}

class trinBank(object):
    
    def __init__(self, trinBankFile):
        self.trinBankFile = trinBankFile
        openFile = open(trinBankFile, 'r')
        first = []
        second = []
        for line in openFile:
            line = line.strip('\n')
            items = line.split()
            first.append(items[0])
            second.append(items[1])
            
        self.gene_id = second[0]
        self.trancript_id = second[1]
        self.sprot_Top_BLASTX_hit = second[2]
        self.RNAMMER = second[3]
        self.prot_id = second[4]
        self.prot_coords = second[5]
        self.sprot_Top_BLASTP_hit = second[6]
        self.Pfam = second[7]
        self.SignalP = second[8]
        self.TmHMM = second[9]
        self.eggnog = second[10]
        self.Kegg = second[11]
        self.gene_ontology_blast = second[12]
        self.gene_ontology_pfam = second[13]
        self.transcript = second[14]
        self.peptide = second[15]
        
        if self.prot_coords != '.':
            self.translateProtein()

        
    def translateProtein(self):
        regex = '\[(\d+)-(\d+)([+-])'
        if self.peptide == "." and self.prot_coords != ".":
            protCoords = re.search(regex, self.prot_coords).groups()
            start = int(protCoords[0])
            end = int(protCoords[1])
            wobble = start + 3   
            codon = None
            translate = ''
            sequence = self.transcript
            if protCoords[2] == '-':
                sequence = self.reverseComplement()
            while wobble < end:
                codon = sequence[start:wobble]
                #codonList.append(codon)
                try:
                    translate+=codons[codon.lower()]
                except:
                    translate+="N"
                start +=3
                wobble +=3
            self.peptide = translate
        
    def reverseComplement(self):
        rev = ''
        seq = self.transcript
        for i in range(1, len(seq)+1):
            if seq[-i].lower() == 'a':
                n = 't'
            elif seq[-i].lower() == 't':
                n = 'a'
            elif seq[-i].lower() == 'g':
                n = 'c'
            elif seq[-i].lower() == 'c':
                n = 'g'
            else:
                n = 'n'
            rev += n
        return rev

newObject = trinBank('/home/scott/Dropbox/Programming/trinBank/trinBank')

print newObject.gene_id
print newObject.peptide