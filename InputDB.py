import os
from datetime import datetime
import multiprocessing as mp
import pickle
from Bio import SeqIO
import math
from gor4 import GOR4
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import sys
import copy
from Bio.Align import substitution_matrices
from Bio import Align
import pandas as pd


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


fn = os.path.join(os.getcwd(), "FeatureExtraction.csv")

CUV_5 ={}
CUV_7 = {}
print("Have to load CUVs from file. This should not happen that often!")
#CUV_5 = pickle.load(open('venv/obj/MiddleCUV_5.pkl', 'rb'))
#CUV_7 = pickle.load(open('venv/obj/MiddleCUV_7.pkl', 'rb'))

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


# Predicts the secondary structure of a given Protein sequence using GOR
# input: Protein Sequence
# output: dict {predicted SS/probabilities : SS Seq / [prob values for each assignment]}
def get_ss_from_seq(seq: str):

    gor4 = GOR4()
    result = gor4.predict(seq)
    #print('Predicted secondary structure', result['predictions'])
    #print('Prediction probabilities', result['probabilities'])
    return(result)



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


# calculation of the match score and the matched percent of an input peptide and a fragment from the CUV file
# input: peptide for which the codon should be predicted, matched Fragment from the CUV File, aligner
def scoreMatch(inputPep, matchedFrag, aligner):
    alignments = aligner.align(inputPep, matchedFrag)
    s = alignments[0].score                             # matchscore of the input peptide and the matched fragment
    alignments = aligner.align(inputPep, inputPep)
    m = alignments[0].score                             # sum of corresponding diagonal scores of the input peptide
    p = s/m
    return(s,p)


def worker(cdslist,seq,cutoff, q):

    aligner = Align.PairwiseAligner()
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")

    inputData = {'acid': [],
                 'cc5': [],
                 'cc7': [],
                 'ss': [],
                 'Result': [],
                 'instable': [],
                 'gravy': [],
                 'polar': [],
                 'pc5': [],
                 'access': []}

    start = datetime.now()
    for cds in cdslist:

        peptide = str(cds.qualifiers['translation'][0])
        genseq = cds.extract(seq)

        secStruc = get_ss_from_seq(peptide)["predictions"]
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
                    # TEST! PUT BACK ASAP
                    #cc7 = getCCVector(peptide[i - 3:i + 4], aligner,cutoff) # If window size 7 possible, calculate cc7 values
                    cc7 = None
                    gravy = prot_analysed.gravy()
                    access = get_solvent_accessibility(peptide[i - 3:i + 4])

                # TEST! PUT BACK ASAP
                cc5 = None
                #cc5 = getCCVector(peptide[i - 2:i + 3], aligner,cutoff)     # get cc5 values for all possible fragments
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

    df = pd.DataFrame(inputData)


    # ALL THAT IS MISSING IS PUTTING DF IN QUEUE AND WRITING IT THE CORRECT WAY TO FILE IN LISTENER()
    done = str(round((datetime.now() - start).total_seconds()))

    res = '#Cds: ' + str(len(cdslist)), "seconds: " + done
    q.put(df)
    return df

def listener(q):
    '''listens for messages on the q, writes to file. '''

    with open(fn, 'a+') as f:
        while 1:
            m = q.get()
            if m == 'kill':
                break
            # TEST! REMOVE FIRST ROW ASAP
            #f.write(str(m) + '\n')
            f.write('acid','cc5','cc7','ss','Result','instable','gravy','polar','pc5','access',+ '\n')
            f.write(m)
            f.flush()


def generate_input_df(cds_list, seq, cutoff):
    start = datetime.now()

    # must use Manager queue here, or will not work
    manager = mp.Manager()
    q = manager.Queue()
    pool = mp.Pool(mp.cpu_count() + 2)
    processes = 100

    # put listener to work first
    watcher = pool.apply_async(listener, (q,))

    # fire off workers
    jobs = []
    chunks = int(math.ceil(len(cds_list)) / processes)
    for i in range(processes):
        job = pool.apply_async(worker, (cds_list[chunks*i:chunks*(i+1)], seq, cutoff, q))
        jobs.append(job)

    # collect results from the workers through the pool result queue
    for job in jobs:
        job.get()

    # now we are done, kill the listener
    q.put('kill')
    pool.close()
    pool.join()

    end = datetime.now()
    print("duration : {} seconds".format(str(round((end - start).total_seconds()))))

def main():
    subset2files = "/Users/friedi/Documents/Cupriavidus necator/Testset/small"
    subset3cds = "Users/friedi/Documents/Pichia_data/TestSets/TestGB"


    dirpath = subset2files
    cutoff = 0.8
    #global CUV_5
    #CUV_5 = pickle.load(open('venv/obj/MiddleCUV_5.pkl', 'rb'))
    #global CUV_7
    #CUV_7 = pickle.load(open('venv/obj/MiddleCUV_7.pkl', 'rb'))



    nr = 0
    for path, subdirs, files in os.walk(dirpath):
        start = datetime.now()
        for file in files:
            nr += 1
            filepath = os.path.join(path, file)
            if file.endswith(".txt") or file.endswith(".embl"):
                record = SeqIO.read(filepath, "embl")
            elif file.endswith("gb") or file.endswith(".gbff"):
                record = SeqIO.read(filepath, "genbank")
            else:
                print("File format can not be processed. Skip file: \n{}".format(filepath))
                continue

            print("Processing genome {}".format(file))
            tag2CDS = get_cds(record)  # create dictionary of relevant cds
            seq = record.seq

            generate_input_df(list(tag2CDS.values()),seq, cutoff)

            end = datetime.now()
            duration = str(round((end - start).total_seconds())/60)
            print("{} files processed! Took me {} minutes.\nHeading to the next one".format(str(nr),duration))


if __name__ == "__main__":
   main()