from Bio import SeqIO
from Bio import GenBank
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import GC123
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.Align import substitution_matrices
from Bio import Align
from collections import OrderedDict
from collections import Counter
import itertools
import copy
import pickle
import os
from gor4 import GOR4
from datetime import datetime
import math
import sys
import pandas as pd

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

SynonymousCodons_3L= {
    "CYS": ["TGT", "TGC"],
    "ASP": ["GAT", "GAC"],
    "SER": ["TCT", "TCG", "TCA", "TCC", "AGC", "AGT"],
    "GLN": ["CAA", "CAG"],
    "MET": ["ATG"],
    "ASN": ["AAC", "AAT"],
    "PRO": ["CCT", "CCG", "CCA", "CCC"],
    "LYS": ["AAG", "AAA"],
    "STOP": ["TAG", "TGA", "TAA"],
    "THR": ["ACC", "ACA", "ACG", "ACT"],
    "PHE": ["TTT", "TTC"],
    "ALA": ["GCA", "GCC", "GCG", "GCT"],
    "GLY": ["GGT", "GGG", "GGA", "GGC"],
    "ILE": ["ATC", "ATA", "ATT"],
    "LEU": ["TTA", "TTG", "CTC", "CTT", "CTG", "CTA"],
    "HIS": ["CAT", "CAC"],
    "ARG": ["CGA", "CGC", "CGG", "CGT", "AGG", "AGA"],
    "TRP": ["TGG"],
    "VAL": ["GTA", "GTC", "GTG", "GTT"],
    "GLU": ["GAG", "GAA"],
    "TYR": ["TAT", "TAC"],
}

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

def cmdParser():
    print("This is where the cmd parser will be located")


def playArea(dirpath, aligner):
    print("Lets play around a little :)")

    cutoff = 0.8

    inputData = {'acid': [],
                 'cc5': [],
                 'cc7': [],
                 'ss':[],
                 'Result':[]}

    for path, subdirs, files in os.walk(dirpath):
        for file in files:
            filepath = os.path.join(path, file)
            if (file.endswith(".txt") or file.endswith(".embl")):
                record = SeqIO.read(filepath, "embl")
            elif (file.endswith("gb") or file.endswith(".gbff")):
                record = SeqIO.read(filepath, "genbank")
            else:
                print("File format can not be processed. Skip file: \n{}".format(filepath))
                continue

            print("Processing genome {}".format(file))

            tag2CDS = getCDS(record)  # create dictionary of relevant cds
            seq = record.seq


            start = datetime.now()


            for cds in tag2CDS.values():
                peptide = str(cds.qualifiers['translation'][0])
                genseq = cds.extract(seq)

                secStruc = getSSfromSeq(peptide)["predictions"]
                print(peptide)
                print(secStruc)

                for i in range(len(peptide)):
                    if not i < 2 or i > (len(peptide)-2):               # If first/last two Acids choose highest frequency codon. Exclude for now! Later None values?
                        if(i < 3 or i == (len(peptide)-3)):              # Only calculate CC5. Fragment Size = 7 not possible!
                            cc7 = None                                  # Assign None value to cc7
                        else:
                            cc7 = getCCVector(peptide[i - 3:i + 4], aligner, cutoff)    # If window size 7 possible, calculate cc7 values

                        cc5 = getCCVector(peptide[i - 2:i + 3], aligner, cutoff)        # get cc5 values for all possible fragments
                        acid = peptide[i]
                        ss = secStruc[i]                                                # Get sec Struc from specific acid
                        codonseq = genseq[i * 3:(i* 3)+3]                               # Get real Codon choice from genomic sequence

                        inputData["acid"].append(acid)  # Add Amino acid to input Vector
                        inputData["cc5"].append(cc5)
                        inputData["cc7"].append(cc7)
                        inputData["ss"].append(ss)
                        inputData["Result"].append(codonseq)

            df = pd.DataFrame(inputData)
            print(df)





            #singlefragCollection_5 = splitCDS_middlecodon(seq=genomeSeq, w=5, cds=list(tag2CDS.values()))  # windowsize 5
            #fragCodonCollection_5 += singlefragCollection_5

            #singlefragCollection_7 = splitCDS_middlecodon(seq=genomeSeq, w=7, cds=list(tag2CDS.values()))  # windowsize 7
            #fragCodonCollection_7 += singlefragCollection_7

            end = datetime.now()
            print("duration : {} seconds".format(str(round((end - start).total_seconds()))))




    ### ABOUT GB RECORDS
    # if feature.type == "CDS":
    #    print(feature.qualifiers['translation'])       # Amino Acid sequence
    # print(record.seq)                              # Holds genome sequence
    # print(feature.type)                    # holds information about the type (gene/mRNA/CDS)
    # print(feature.location)                # location information on the genome sequence
    # print(feature.location.strand)         # strand on the sequence that the feature is located on
    # print(feature.qualifiers)
    # print(feature.extract(record.seq))     # Genomic sequence of the gene


# Quick check if len(genomic Sequence) == 3* len(protein sequence)
# Remove proteins with len < 100
def qualitiycheck(tag2CDS, record):
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
    with open('obj/' + name + '.pkl', 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)

def load_obj(name ):
    with open('obj/' + name + '.pkl', 'rb') as f:
        return pickle.load(f)


# Input: gb record that holds the information of one whole gb file.
# Output: dictionary that only contains gene and CDS features, of those that share the same location
#        { <locusTag> : <CDS>}
def getCDS(record):

    tag2CDS = {}

    for feature in record.features:
        if feature.type == "CDS":
            if "translation" in feature.qualifiers.keys():
                if len(feature.qualifiers['translation'][0]) > 100:
                    tag2CDS[feature.qualifiers['locus_tag'][0]] = feature
            else:
                continue
    return tag2CDS


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


# splits all CDS into fragments of windowsize w and collects the corresponding codons
# Input: genomic sequence seq, windowsize w, list of CDS cds
# Output: List of [<proteinfragment>, [<codon>]] elements
def splitCDS_fullCodonSet(seq, w:int, cds:list):
    print("Collecting information for fragments of size: {}".format(w))

    fragmentContent= []

    for protein in cds:
        peptide = str(protein.qualifiers['translation'][0])
        genseq= protein.extract(seq)#.replace("T","U")       # Replace T with U to match common Codons

        for i in range(len(peptide)-w):
            fragment = peptide[i:i+w]
            codonseq= genseq[i*3:(i+w)*3]
            codonlist = [str(codonseq[j:j+3]) for j in range(0, len(codonseq),3)]
            fragmentContent.append((fragment,codonlist))

    return(fragmentContent)



# Calculation of the Codon Usage Vector for the middle amino acid in similar fragments
# Creates pickled Dict that holds the Codon Usage vector for each Fragment.
# The Vector contains usage values for all codons
def createMiddleCUVfile_FullCodonSet(fragment2codonlist: dict, windowsize:int):

    print("Calculating codon usage pattern")
    codonUsageSingle = copy.deepcopy(CodonsDict)
    codonUsageAll = {}

    for fragment in fragment2codonlist:
        nrFragments = len(fragment2codonlist[fragment])
        for codonlist in fragment2codonlist[fragment]:
            middlecodon = codonlist[0][int(windowsize/2)]
            if not middlecodon in codonUsageSingle.keys():
                print("Middle Codon is illegal Codon. Skip Codon {}".format(middlecodon))
                nrFragments = nrFragments-1
                continue
            codonUsageSingle[middlecodon] = codonUsageSingle[middlecodon] +(1/nrFragments)  #save codon usage for specific Acid

        codonUsageAll[fragment] = codonUsageSingle          # save created vector for the fragment
        codonUsageSingle = copy.deepcopy(CodonsDict)


    # save created dictionary to file using pickle
    save_obj(codonUsageAll, "MiddleCUV_"+ str(windowsize))

    # reloaded_CodonUsageAll = load_obj("MiddleCUV")      # reload data



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


# Creates Codon usage Vectors for the middle amino acids of size 5 and 7 from a collection of annotated genomes
# Creates Average Usage Vectors for all Codons in all sequences
# Input: Directory that contains all Genome files as genbank or embl
# Result: pickled files containing the Codon Usage information of the middle Codon with w = 5 and w = 7,
#       and a pickle file containing the average usage for each amino acid in all sequences
def generateFragmentDB(dirpath):

    # initialize Lists for all found fragments and the corresponding Codons. Same fragment can occur multiple times
    fragCodonCollection_5 = []                      # windowsize = 5
    fragCodonCollection_7 = []                      # windowsize = 7

    sequences = []

    for path, subdirs, files in os.walk(dirpath):
        for file in files:
            filepath = os.path.join(path, file)
            if (file.endswith(".txt") or file.endswith(".embl")):
                record = SeqIO.read(filepath, "embl")
            elif(file.endswith("gb")or file.endswith(".gbff")):
                record = SeqIO.read(filepath, "genbank")
            else:
                print("File format can not be processed. Skip file: \n{}".format(filepath))
                continue

            print("Processing genome {}".format(file))

            genomeSeq = record.seq                  # extract genome Sequence
            sequences.append(genomeSeq)             # add genome sequence to list
            tag2CDS = getCDS(record)                # create dictionary of relevant cds
            qualitiycheck(tag2CDS, record)          # check if length of genomic and protein seq match. Remove short proteins

            # split CDS into fragments and collect codon occurences
            fillFragment2Middlecodon(genomeSeq, tag2CDS.values(), 5)
            fillFragment2Middlecodon(genomeSeq, tag2CDS.values(), 7)


    createMiddleCUVfile_middlecodon(Fragment2Middlecodon_5, 5)  # windowsize = 5
    createMiddleCUVfile_middlecodon(Fragment2Middlecodon_7, 7)  # windowsize = 7

    # create average usage for each amino acid from the genome sequences
    createAcidCUV_file(sequences)
    calculateSCUO()
    calculateRSCU()
    gcsum = 0
    gc3sum = 0
    seqnum = len(sequences)
    for seq in sequences:
        gc, gc3 = getGCContents(seq)
        gcsum += gc
        gc3sum += gc3
    gcContent = {gc: gcsum/ seqnum, gc3: gc3sum/seqnum}

    save_obj(gcContent, "GCcontent")



# Predicts the secondary structure of a given Protein sequence using GOR
# input: Protein Sequence
# output: dict {predicted SS/probabilities : SS Seq / [prob values for each assignment]}
def getSSfromSeq(seq: str):

    gor4 = GOR4()
    result = gor4.predict(seq)
    #print('Predicted secondary structure', result['predictions'])
    #print('Prediction probabilities', result['probabilities'])
    return(result)



# Calculation of the overall GC Content and the GC3 content
# Input: DNA sequence
# Output: gc, gc3
def getGCContents(seq:str):
    gcContent = GC123(seq)
    return gcContent[0], gcContent[3]


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


# calculation of the match score and the matched percent of an input peptide and a fragment from the CUV file
# input: peptide for which the codon should be predicted, matched Fragment from the CUV File, aligner
def scoreMatch(inputPep, matchedFrag, aligner):
    alignments = aligner.align(inputPep, matchedFrag)
    s = alignments[0].score                             # matchscore of the input peptide and the matched fragment
    alignments = aligner.align(inputPep, inputPep)
    m = alignments[0].score                             # sum of corresponding diagonal scores of the input peptide
    p = s/m
    return(s,p)


# Calculation of the final input vector for the middle amino acid in a peptide, based on Codon Context
# input: fragment of the input Protein sequence of size 5 or 7. Aligner, cutoff
# Output: {<Syn Codon> : <weighted frequency>}
def getCCVector(inputPep, aligner, cutoff):
    """

    :param inputPep: Input Peptide of size 5 or 7
    :param aligner: Aligner using Blosum to get matching score of input peptide and fragments from DB
    :param cutoff: cutoff atwhich match percent value a fragment from the DB should ne considered
    :return: Codon Context Vector {<Synonymus Codon>:<weighted frequency>}
    """
    w = len(inputPep)

    if w == 5:
        middle = 2
        matchFragments = CUV_5[inputPep[middle]]

    elif w == 7:
        middle = 3
        matchFragments = CUV_7[inputPep[middle]]
    else:
        sys.exit("Wrong size of input peptide! Peptide: {}".format(inputPep))

    inputVector = copy.deepcopy(SynonymousCodons)[inputPep[middle]]
    successfullMatches = 0

    for fragment in matchFragments.keys():
        s,p = scoreMatch(inputPep,fragment,aligner)
        if p > cutoff:
            successfullMatches += 1
            codons = matchFragments[fragment]
            for codon in codons:
                inputVector[codon] += codons[codon]*s

    for codon in inputVector:
        if successfullMatches == 0:
            successfullMatches = 1          # Prevent division by zero

        inputVector[codon] = inputVector[codon]/successfullMatches
    return(inputVector)

# Calculation of the fraction of polar amino acids in a fragment
def get_polar_fraction(fragment):
    """

    :param fragment: Peptide for which the fraction of polar acids should be calculated
    :return: fraction of polar amino acids in the given fragment
    """
    nr_polar = 0
    size = len(fragment)
    for i in range(size):
        if fragment[i] in ["A","G","T","S","N","Q","D","E","H","R","K","P"]:
            nr_polar += 1

    return(nr_polar/size)


# Create a label that determines the direct environment of a specific amino acid
def get_pc5(triplet):
    """

    :param triplet: Tripple of currently considered Amino acid and both surrounding ones
    :return: Neighbors Reduced to 5 letter alphabet according to physio-chemical properties
        A (Aliphatic: IVL), R (aRomatic: FYWH), C (Charged: KRDE),T (Tiny: GACS), D (Diverse: TMQNP)
    """

    alipathic = ["I","V","L"]
    aRomatic = ["F","Y","W","H"]
    charged = ["E","D","R","K"]
    divers = ["T","M","Q","N","P"]
    tiny = ["G","A","C","S"]

    neighbours = []

    for i in [0,2]:
        if triplet[i] in alipathic:
            neighbours.append("A")
        elif triplet[i] in aRomatic:
            neighbours.append("R")
        elif triplet[i] in charged:
            neighbours.append("C")
        elif triplet[i] in divers:
            neighbours.append("D")
        elif triplet[i] in tiny:
            neighbours.append("T")
        else:
            print("Acid is none of known amino acids: {}".format(triplet[i]))
            neighbours.append("N")

    neighbours.sort()
    return("".join(neighbours))


# create solvent accessibility score of a fragment, based on clusters of the amino acids defined by Tomii and Kanehisa
def get_solvent_accessibility(fragment):
    """

    :param fragment: peptide fragment with the current amino acid in the middle
    :return: solvent accessibility score, calculated by assigning each amino acid to a class and score
    """
    buried = ["A","L","F","C","G","I","V","W"]
    exposed = ["P","K","Q","E","N","D"]
    intermediate = ["M","P","S","T","H","Y" ]

    score = 0
    for i in range(len(fragment)):
        if fragment[i] in buried:
            score += 0
        elif fragment[i] in exposed:
            score += 1
        elif fragment[i] in intermediate:
            score += 0.5

    return score


def generateInputDF(dirpath, aligner, cutoff):

    inputData = {'acid': [],
                 'cc5': [],
                 'cc7': [],
                 'ss': [],
                 'Result': [],
                 'instable': [],
                 'gravy': [],
                 'polar':[],
                 'pc5': [],
                 'access': []}
    nr = 0

    for path, subdirs, files in os.walk(dirpath):
        for file in files:
            filepath = os.path.join(path, file)
            if (file.endswith(".txt") or file.endswith(".embl")):
                record = SeqIO.read(filepath, "embl")
            elif (file.endswith("gb") or file.endswith(".gbff")):
                record = SeqIO.read(filepath, "genbank")
            else:
                print("File format can not be processed. Skip file: \n{}".format(filepath))
                continue

            print("Processing genome {}".format(file))

            tag2CDS = getCDS(record)  # create dictionary of relevant cds
            seq = record.seq




            for cds in tag2CDS.values():
                nr += 1
                peptide = str(cds.qualifiers['translation'][0])
                genseq = cds.extract(seq)

                secStruc = getSSfromSeq(peptide)["predictions"]
                prot_analysed = ProteinAnalysis(peptide)


                for i in range(len(peptide)):
                    if i< 1 or i >(len(peptide) - 2):
                        pc5 = None
                    else:
                        pc5 = get_pc5(peptide[i-1:i+2])

                    if i < 2 or i > (len(peptide) - 3):     # If first/last two Acids assign None values, since window not possible
                        cc7 = None
                        cc5 = None
                        gravy = None                        # none for gravy because window is too small
                        polarfrac = None                    # None for fraction of polar acids
                        access = None                       # None for solvent accessibilty

                    else:
                        if (i < 3 or i == (len(peptide) - 3)):  # Only calculate CC5. Fragment Size = 7 not possible!
                            cc7 = None                          # Assign None value to cc7
                            gravy = None                           # none for gravy because window is too small
                            access = None
                        else:
                            cc7 = getCCVector(peptide[i - 3:i + 4], aligner,cutoff) # If window size 7 possible, calculate cc7 values
                            gravy = prot_analysed.gravy()
                            access = get_solvent_accessibility(peptide[i - 3:i + 4])

                        cc5 = getCCVector(peptide[i - 2:i + 3], aligner,cutoff)     # get cc5 values for all possible fragments
                        polarfrac = get_polar_fraction(peptide[i - 2:i + 3])

                    acid = peptide[i]
                    ss = secStruc[i]                                    # Get sec Struc from specific acid
                    codonseq = genseq[i * 3:(i * 3) + 3]                # Get real Codon from genomic sequence


                    # Add Values to input Vector
                    inputData["acid"].append(acid)
                    inputData["cc5"].append(cc5)
                    inputData["cc7"].append(cc7)
                    inputData["ss"].append(ss)
                    inputData["Result"].append(codonseq)
                    inputData["instable"].append(prot_analysed.instability_index())
                    inputData["gravy"].append(gravy)
                    inputData["polar"].append(polarfrac)
                    inputData["pc5"].append(pc5)
                    inputData["access"].append(access)

                print("{} processed.".format(nr))




    # Convert Collected Values to InputDataframe
    # TODO: Include global Background Information
    df = pd.DataFrame(inputData)

    save_obj(df, "Input_Data")
    return(df)



def createModel(df):
    print("Let the fun begin!")

# Wrap it all up
def main():

    # Hard coded section that needs to be exchanged via arguments from the command line parser.
    # names have to remain the same
    #------------------------------------------------------------

    dirpath = "/Users/friedi/Documents/Cupriavidus necator/Databank"
    modelinput = "/Users/friedi/Documents/Cupriavidus necator/ModelInput"
    filepath = "/Users/friedi/Documents/Cupriavidus necator/Databank/CP066019.1.gb"
    protseq = "MSNFVALGDHITVQKGKAPLVTGYVGKGAEPYLSPEYLRGRAPADLAKAGPDAVRAAEGETILLWDGSNAGEFFRSKVGLVASTMTKISPSSVFRPAYFFHVAKQAERFLKAQTNGTGIPHVDRELLEGIKVFCPGSTEQQLLAEILDTLDTAIYETEAIIAKLKAVKQGLLHDLLTRGIDANGELRPPQAEAPHLYESSPLGWIPNEWGLAPTATRCHLITKGTTPAANEMWQGGAGIRFLRVDNLSFDGQLDLDASTFRVSLATHKGFLARSRCLEGDVLTNIVGPPLGKLGLVTKEIGEVNINQAIALFRPTEQLLPKFLLIWLSSSISQSWLRNRAKQTSGQVNLTLALCQELPLPRMTINEQQAIVDRVDAAQEQIWCEEELIRKMRLEKSGLMDDLLTGRVRVKPQLAETKQAGSA"
    genSeq = "TCTGTCTTGCTGCCGGCCGTGTTTATGCATGTGCAGGTCGCTAACCTGTTAGCAAGGACGGCACGAACAGGCCAGACCCGGCGATGCGCATACGTGATCCCCTGCCGCATCGCAGAACCCGGCATTTAACCGATTTTCTGCGGCAATGTCAAGGCGCAAGACTATGTGCACCGGCACAATTGCCACACCCGTCTGTGGCCCGTGGTGCGC"

    createCUV = False
    createInputDB = True
    cutoff = 0.8
    #------------------------------------------------------------

    aligner = Align.PairwiseAligner()
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")

    start = datetime.now()
    print("start collecting background informations")
    # Create background DB and Vectors from genomes located in dirpath or load background values from pickled objects
    if createCUV:
        generateFragmentDB(dirpath)
    else:
        ## Load Species specific background information, generated by species genome data

        global ICU; ICU = load_obj("AminoAcidUsage")            #Individual Codon Usage
        global CUV_5; CUV_5 = load_obj("MiddleCUV_5")           # Codon usage Vector, w = 5
        global CUV_7; CUV_7 = load_obj("MiddleCUV_7")           # Codon usage Vector, w = 7
        global RSCU; RSCU = load_obj("RSCU")                    # Relative synonymous codon usage
        global SCUO; SCUO= load_obj("SCUO")                     # Synonymous codon usage orderliness
        global GCcontent; GCcontent = load_obj("GCcontent")

    end = datetime.now()
    print("duration : {} seconds".format(str(round((end - start).total_seconds()))))

    ## Create Database from Trainings/Test data
    print(CUV_5)
    print(CUV_5.keys())

    print("start filling up the input dataframe")
    start = datetime.now()
    if createInputDB:
        input_data = generateInputDF(modelinput, aligner, cutoff)
    else:
        input_data = load_obj("Input_Data")

    end = datetime.now()
    print("duration : {} minutes".format(str(round((end - start).total_seconds())/60)))

    createModel(input_data)

    print("This is a change!")







    #start = datetime.now()

    #end = datetime.now()
    #print("Total duration  : {} minutes".format(str(round((end - start).total_seconds())/60)))


if __name__ == '__main__':
    print("let the magic happen")
    main()