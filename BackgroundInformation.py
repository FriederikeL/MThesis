from Bio import SeqIO
from Bio.SeqUtils import GC123
from Bio.Align import substitution_matrices
from Bio import Align
from collections import Counter
import copy
import pickle
import os
from datetime import datetime
import math




SynonymousCodons = {
    "C": {"TGT":0, "TGC":0},
    "D": {"GAT":0, "GAC":0},
    "S": {"TCT":0, "TCG":0, "TCA":0, "TCC":0, "AGC":0, "AGT":0},
    "Q": {"CAA":0, "CAG":0},
    "M": {"ATG":0},
    "N": {"AAC":0, "AAT":0},
    "P": {"CCT":0, "CCG":0, "CCA":0, "CCC":0},
    "K": {"AAG":0, "AAA":0},
    "STOP": {"TAG":0, "TGA":0, "TAA":0},
    "T": {"ACC":0, "ACA":0, "ACG":0, "ACT":0},
    "F": {"TTT":0, "TTC":0},
    "A": {"GCA":0, "GCC":0, "GCG":0, "GCT":0},
    "G": {"GGT":0, "GGG":0, "GGA":0, "GGC":0},
    "I": {"ATC":0, "ATA":0, "ATT":0},
    "L": {"TTA":0, "TTG":0, "CTC":0, "CTT":0, "CTG":0, "CTA":0},
    "H": {"CAT":0, "CAC":0},
    "R": {"CGA":0, "CGC":0, "CGG":0, "CGT":0, "AGG":0, "AGA":0},
    "W": {"TGG":0},
    "V": {"GTA":0, "GTC":0, "GTG":0, "GTT":0},
    "E": {"GAG":0, "GAA":0},
    "Y": {"TAT":0, "TAC":0},
}
CodonsDict = {
    "TTT": 0, "TTC": 0, "TTA": 0, "TTG": 0,
    "CTT": 0, "CTC": 0, "CTA": 0, "CTG": 0,
    "ATT": 0, "ATC": 0, "ATA": 0, "ATG": 0,
    "GTT": 0, "GTC": 0, "GTA": 0, "GTG": 0,
    "TAT": 0, "TAC": 0, "TAA": 0, "TAG": 0,
    "CAT": 0, "CAC": 0, "CAA": 0, "CAG": 0,
    "AAT": 0, "AAC": 0, "AAA": 0, "AAG": 0,
    "GAT": 0, "GAC": 0, "GAA": 0, "GAG": 0,
    "TCT": 0, "TCC": 0, "TCA": 0, "TCG": 0,
    "CCT": 0, "CCC": 0, "CCA": 0, "CCG": 0,
    "ACT": 0, "ACC": 0, "ACA": 0, "ACG": 0,
    "GCT": 0, "GCC": 0, "GCA": 0, "GCG": 0,
    "TGT": 0, "TGC": 0, "TGA": 0, "TGG": 0,
    "CGT": 0, "CGC": 0, "CGA": 0, "CGG": 0,
    "AGT": 0, "AGC": 0, "AGA": 0, "AGG": 0,
    "GGT": 0, "GGC": 0, "GGA": 0, "GGG": 0}
AminoAcids = {
    "C": 0, "D": 0, "S": 0,"Q": 0, "M": 0, "N": 0, "P": 0, "K": 0, "STOP": 0, "T": 0,
    "F": 0, "A": 0, "G": 0, "I": 0, "L": 0, "H": 0, "R": 0, "W": 0, "V": 0, "E": 0, "Y": 0,
}


# dict: {<fragment:str>: {<codon:str>: <codonCount:int>}}
Fragment2Middlecodon_5 = {}
Fragment2Middlecodon_7 = {}
CodonCounts = {}

ICU = {}
CUV_5 = {}
CUV_7 = {}
RSCU = {}
SCUO = {}
GCcontent = {}


# Quick check if len(genomic Sequence) == 3* len(protein sequence)
# Remove proteins with len < 100
def qualitycheck(tag2CDS, record):
    for entry in tag2CDS:
        cds = (tag2CDS[entry])
        frac = len(cds.extract(record.seq)) / (len(cds.qualifiers['translation'][0]))
        problems = 0
        if not (round(frac) == 3):
            print("Doesnt fit")
            print(frac)
            print(cds.extract(record.seq))
            print(cds.qualifiers['translation'])
            problems += 1
    if problems == 0:
        print("Passed Quality check without any problems")
    else:
        print("Some quality problems occured in {} cases. See problematic sequences above".format(str(problems)))


def save_obj(obj, name):
    if not os.path.exists("obj"):
        os.mkdir("obj")
    with open('venv/obj/' + name + '.pkl', 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)


# Input: gb record that holds the information of one whole gb file.
# Output: dictionary that only contains gene and CDS features, of those that share the same location
#        { <locusTag> : <CDS>}
def get_cds(record):

    tag2CDS = {}

    for feature in record.features:
        if feature.type == "CDS":
            if "translation" in feature.qualifiers.keys():
                if len(feature.qualifiers['translation'][0]) > 100:
                    tag2CDS[feature.qualifiers['locus_tag'][0]] = feature
            else:
                continue
    return tag2CDS


# Calculation of the Codon Usage Vector for the middle amino acid in similar fragments
# Creates pickled Dict that holds the Codon Usage vector for each Fragment, clustered by miffe acid.
# dict: { mAcid{<fragment>:{<codon>:<codonUsage>}}}
# The Vector only contains information about the codons of the middle amino acid
def createMiddleCUVfile_middlecodon(fragment2codonlist: dict, windowsize:int):

    print("Calculating codon usage pattern")
    #codonUsageAll = {}
    legalCodons = CodonsDict.keys()
    acid2Fragments = {}

    for fragment in fragment2codonlist:
        nrFragments = sum(fragment2codonlist[fragment].values())    # collect number of Fragments
        acid = fragment[int(windowsize/2)]                          # Middle acid in the fragment
        codons4acid = copy.deepcopy(SynonymousCodons[acid])         # slice SynonymousCodons by acid key

        for middlecodon in fragment2codonlist[fragment]:

            if(middlecodon in legalCodons):
                codons4acid[middlecodon] = codons4acid[middlecodon] + (1 / nrFragments)     # save codon usage for specific Acid
            else:
                print("Middle Codon is illegal Codon. Skip Codon {}".format(middlecodon))   # Remove if codon is not one of the 64 natural codons (61 Aa,3 stop)
                nrFragments = nrFragments-1
                continue

        #codonUsageAll[fragment] = codons4acid  # save created vector for the fragment

        if acid in acid2Fragments.keys():
            acid2Fragments[acid][fragment] = codons4acid        # save created vectors for the fragment, clustered by middle acid
        else:
            acid2Fragments[acid] = {fragment: codons4acid}

    # save created dictionary to file using pickle
    save_obj(acid2Fragments, "MiddleCUV_"+ str(windowsize))


# Calculation of the average codon usage for each amino acid in the whole Sequence
# Input: genomic sequence
# Creates pickled Dict that holds average codon usage per amino acid {acid: {codon: fraction}}
def createAcidCUV_file(sequences : list):
    global CodonCounts

    acidCUV = copy.deepcopy(SynonymousCodons)
    codonlist = []

    for seq in sequences:
        codonlist +=  [str(seq[j:j + 3]) for j in range(0, len(seq), 3)]  # split sequence into codons

    CodonCounts = Counter(codonlist)                                    # calculate occurences per codon

    # calculate average codon usage for each amino acid
    for acid in SynonymousCodons:
        acidCodes = list(SynonymousCodons[acid].keys())                 # get all codes for the specific amino acid
        acidCodeCounts = {key:CodonCounts[key] for key in CodonCounts.keys() & acidCodes}   # get the occurences of these codes
        acidOccurence = sum(acidCodeCounts.values())

        # Identify STOP Codons
        if acid == "STOP":
            acidCodeCounts = {key: CodonCounts[key] for key in
                              CodonCounts.keys() & ["TGA", "TAA", "TAG"]}  # get the occuences of these codes
            acidOccurence = sum(acidCodeCounts.values())

        # get number of total occurences of the amino acid
        for code in acidCUV[acid]:
            acidCUV[acid][code] += acidCodeCounts[code]/acidOccurence

    # save created dictionary to file using pickle
    save_obj(acidCUV, "AminoAcidUsage")

    #reloaded_acidCUV = load_obj("AminoAcidUsage")      # reload data


# splits all CDS into fragments of windowsize w and collects the codon occurences for the middle amino acid
# Input: genomic sequence seq, windowsize w, list of CDS cds
# Result: fills the global fragment2Middlecodon_w for given windowsize w. {<fragment:str>: {<codon:str>: <codonCount:int>}}
def fillFragment2Middlecodon(seq, cds:list, w:int):

    # assign start position of middle codon for respective windowsize, assign dictionary to be filled.
    if w == 5:
        startPos = 6
        fragment2Middlecodon = Fragment2Middlecodon_5
    else:
        startPos = 9
        fragment2Middlecodon = Fragment2Middlecodon_7

    for protein in cds:
        peptide = str(protein.qualifiers['translation'][0])
        genseq = protein.extract(seq)

        for i in range(len(peptide)-w):
            fragment = peptide[i:i+w]
            codonseq = genseq[i * 3:(i + w) * 3]
            middlecodon = str(codonseq[startPos:startPos+3])

            if fragment in fragment2Middlecodon.keys():
                fragmentsContent = fragment2Middlecodon[fragment]
                if middlecodon in fragmentsContent.keys():
                    fragment2Middlecodon[fragment][middlecodon] +=1
                else:
                    fragment2Middlecodon[fragment][middlecodon] = 1
            else:
                codondict = {middlecodon: 1}

                fragment2Middlecodon[fragment] = codondict


# rscu = relative synonymous codon usage,
# defined as the ratio of the observed frequency of codons to the expected frequency,
# given that all the synonymous codons for the same amino acid are used equal
# Result: pickled dict {<acid>: {<codon>: <rscu value>}}
# WARNING: Can only be used if createAcidCUV_file() was called before
def calculateRSCU():

    if(sum(CodonCounts.values()) == 0):
        print("Warning! RSCU cant be computed before createAcidCUV_file() was called.")
        return None

    rscuDict = copy.deepcopy(SynonymousCodons)

    for acid in rscuDict:
        s = len(rscuDict[acid])     # number of codons coding for the acid
        countsum = 0
        # get overall occurences of amino acid
        for codon in rscuDict[acid]:
            countsum +=CodonCounts[codon]

        for codon in rscuDict[acid]:
            observed = CodonCounts[codon]
            rscuDict[acid][codon] = observed/(countsum/s)

    # save created dictionary to file using pickle
    save_obj(rscuDict, "RSCU")


# scuo = synonymous codon usage order.
# measurement of the nonrandomness in synonymous codon usage, using information entropy
# The greater the nonrandomness in synonymous codon usage is, the more ordered the sequence will be
# Result: pickled dict {<acid> : <scuo value>}
# WARNING: Can only be used if createAcidCUV_file() was called before
def calculateSCUO():

    if (sum(CodonCounts.values()) == 0):
        print("Warning! RSCU cant be computed before createAcidCUV_file() was called.")
        return None

    acid2SCUO = copy.deepcopy(AminoAcids)
    scuoDict = SynonymousCodons

    for acid_i in scuoDict:
        ni = len(scuoDict[acid_i])              # number of codons coding for the i-th acid
        if(ni == 1):                            # Skip acids with single codons (Met, Trp)

            continue
        countSum = 0
        for codon in scuoDict[acid_i]:
            countSum += CodonCounts[codon]
        entropys = []

        hijsum = 0
        for codon in scuoDict[acid_i]:
            pij = CodonCounts[codon] / countSum      # normalizing the occurences of the current codon
            hij = -pij*math.log(pij, 2)         # Entropy of the i-th amino acid of the j-th codon, according to information theory
            entropys.append(hij)
            hijsum += pij*math.log(pij, 2)      # sum of the entropies of the codons representing amino acid i

        hi  = hijsum*(-1)                       # entropy of the i-th amino acid
        himax = - math.log(1/ni, 2)             # maximum entropy for the i-th amino acid, assuming uniform distribution for random codon usage
        #himin = 0                              # minimum entropy. here: if only one of the syn. codons is used for the i-th amino acid
        ii = himax- hi                          # Information defined as the difference between the maximum entropy and the actual entropy as an index of orderliness
                                                # Here: measurement of the nonrandomness in synonymous codon usage
        oi = ii/himax                           # normalizes difference between maximum entropy and the observed entropy
                                                # oi == 1 if codon usage is biased to the extreme.


        acid2SCUO[acid_i] = oi                  # assign oi as scuo value to each amino acid

    save_obj(acid2SCUO, "SCUO")


# Calculation of the overall GC Content and the GC3 content
# Input: DNA sequence
# Output: gc, gc3
def getGCContents(seq:str):
    gcContent = GC123(seq)
    return gcContent[0], gcContent[3]

# Create background DB and Vectors from genomes located in dirpath and safe as pickled objects
def main():

    # Hard coded section that needs to be exchanged via arguments from the command line parser.
    # names have to remain the same
    #------------------------------------------------------------

    dirpath = "/Users/friedi/Documents/Cupriavidus necator/Databank"
    cutoff = 0.8
    #------------------------------------------------------------

    aligner = Align.PairwiseAligner()
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")

    start = datetime.now()
    print("start collecting background informations")

    sequences = []

    for path, subdirs, files in os.walk(dirpath):
        for file in files:
            filepath = os.path.join(path, file)
            if file.endswith(".txt") or file.endswith(".embl"):
                record = SeqIO.read(filepath, "embl")
            elif file.endswith("gb")or file.endswith(".gbff"):
                record = SeqIO.read(filepath, "genbank")
            else:
                print("File format can not be processed. Skip file: \n{}".format(filepath))
                continue

            print("Processing genome {}".format(file))

            genomeSeq = record.seq                  # extract genome Sequence
            sequences.append(genomeSeq)             # add genome sequence to list
            tag_2_cds = get_cds(record)                # create dictionary of relevant cds
            qualitycheck(tag_2_cds, record)          # check if length of genomic and protein seq match. Remove short proteins

            # split CDS into fragments and collect codon occurences
            fillFragment2Middlecodon(genomeSeq, tag_2_cds.values(), 5)
            fillFragment2Middlecodon(genomeSeq, tag_2_cds.values(), 7)


    createMiddleCUVfile_middlecodon(Fragment2Middlecodon_5, 5)  # windowsize = 5
    createMiddleCUVfile_middlecodon(Fragment2Middlecodon_7, 7)  # windowsize = 7


    createAcidCUV_file(sequences)   # average usage for each amino acid from the genome sequences
    calculateSCUO()
    calculateRSCU()

    # get GC Content values from all input sequences
    gcsum = 0
    gc3sum = 0
    seqnum = len(sequences)
    for seq in sequences:
        gc, gc3 = getGCContents(seq)
        gcsum += gc
        gc3sum += gc3
    gc_content = {gc: gcsum/ seqnum, gc3: gc3sum/seqnum}

    save_obj(gc_content, "GCcontent")

if __name__ == '__main__':
    print("let the magic happen")
    main()