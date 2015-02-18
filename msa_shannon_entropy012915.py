#!/usr/bin/env python

'''

Compiled by Felix Francis (felixfrancier@gmail.com)

Description:    Parse multiple sequence alignment out put files and calculate Shannon's entropy for each column of the alignment.

                Shannon's entropy equation (latex format):
                
                H=-\sum_{i=1}^{M} P_i\,log_2\,P_i
                
                Entropy is a measure of the uncertainty of a probability distribution (p1, ..... , pM)
                https://stepic.org/lesson/Scoring-Motifs-157/step/7?course=Bioinformatics-Algorithms&unit=436
                
                Where, Pi is the fraction of nuleotide bases of nuleotide base type i, and M is the number of nuleotide base types (A, T, G or C)
                
                H ranges from 0 (only one base/residue in present at that position) to 4.322 (all 20 residues are equally represented in that position).
                Typically, positions with H >2.0 are considerered variable, whereas those with H < 2 are consider conserved.
                Highly conserved positions are those with H <1.0 (Litwin and Jores, 1992). A minimum number of sequences is however required (~100)
                for H to describe the diversity of a protein family.
               
'''

#Time to run the code: start timer
import time
t0 = time.time()

# Add an import path to BioPython
import sys
sys.path.append("~/Desktop/Softwares/biopython-1.65/Bio/")

# import pandas for data frame
import pandas as pd
from numpy.random import randn

# import seaborn for violin plot
import seaborn as sns
import matplotlib.pyplot as plt


##################################################################
# Function to calcuate the Shannon's entropy per alignment column
# H=-\sum_{i=1}^{M} P_i\,log_2\,P_i (http://imed.med.ucm.es/Tools/svs_help.html)
# Gaps and N's are included in the calculation
##################################################################
def shannon_entropy(list_input):
    import math
    unique_base = set(list_input)                           # Get only the unique bases in a column
    #unique_base = unique_base.discard("-")
    M   =  len(list_input)
    entropy_list = []
    # Number of residues in column
    for base in unique_base:
        n_i = list_input.count(base)                        # Number of residues of type i                   
        P_i = n_i/float(M)                                  # n_i(Number of residues of type i) / M(Number of residues in column)
        entropy_i = P_i*(math.log(P_i,2))
        entropy_list.append(entropy_i)
    sh_entropy = -(sum(entropy_list))
    #print sh_entropy
    return sh_entropy


##################################################################
# Parse the MSA (BioPython)
##################################################################
from Bio import AlignIO

##################################################################
## Import ".clustal" file
#align = (AlignIO.read("test1.clustal", "clustal"))                          # Test file
align_clustal = AlignIO.read("CLUSTAL_OM011215/clustalo-E20150108-022540-0026-80211746-oy.clustal", "clustal")  # An actual clustal Omega file
##################################################################
#print len(list(align[0]))                                                  # Gives the number of columns in the MSA
#print (list(align[0]))

#for alignment in align:
#    print (alignment)
#    print("")

##################################################################
## Import "MAFFT .clustal" file
align_mafft = AlignIO.read("MAFFT012515/_out150126125256989tRydvIXY3nJkf6i4gp7mK.aln", "clustal")
##################################################################


##################################################################
## Import "MUSCLE .clustal" file
align_muscle = AlignIO.read("MUSCLE012515/muscle-E20150126-035203-0739-5913603-es.clw", "clustal")
##################################################################


###################################################################
### Import "gold standard.clustal" file
align_gold = AlignIO.read("Test_AAA_smart/PHYLIP.cgi", "phylip")
align_gold_co = AlignIO.read("Test_AAA_smart/clustalo_aaa_021715.clustal", "clustal")
align_gold_mafft = AlignIO.read("Test_AAA_smart/mafft_aaa_021715.clustalw", "clustal")
align_gold_muscle = AlignIO.read("Test_AAA_smart/muscle_aaa_021715.clw", "clustal")
align_gold_kalign = AlignIO.read("Test_AAA_smart/kalign_aaa_021715.clustalw", "clustal")
###################################################################



###################################################################
## Error message if the lengths of the sequences in MSA are not the same
###################################################################
#seq_lengths_list = []
#for record in align:
#    seq_lengths_list.append(len(record))
##row_num = len(seq_lengths_list)                                            # Get number of rows in the MSA
#seq_lengths = set(seq_lengths_list)                                         # Get unique lengths of the sequences aligned in the MSA
#if len(seq_lengths) != 1:
#    print "Check you input Alignment!",                                     # Error message if the lengths of the sequences in MSA are not the same


##################################################################
# Function to calculate Shannon's entropy per alignment column for the whole MSA
##################################################################

def shannon_entropy_list_msa(alignment_file):
    shannon_entropy_list = []
    for col_no in xrange(len(list(alignment_file[0]))):
        list_input = list(alignment_file[:, col_no])
        shannon_entropy_list.append(shannon_entropy(list_input))
    return shannon_entropy_list


##################################################################
# Calculate Shannon's entropy for input alignments
##################################################################
MAFFT           = shannon_entropy_list_msa(align_mafft)
CLUSTAL_OMEGA   = shannon_entropy_list_msa(align_clustal)
MUSCLE          = shannon_entropy_list_msa(align_muscle)
GOLD_AA          = shannon_entropy_list_msa(align_gold)
GOLD_CO          = shannon_entropy_list_msa(align_gold_co)
GOLD_MAFFT         = shannon_entropy_list_msa(align_gold_mafft)
GOLD_MUSCLE        = shannon_entropy_list_msa(align_gold_muscle)
GOLD_KALIGN        = shannon_entropy_list_msa(align_gold_kalign)
#print min(float(s) for s in CLUSTAL_OMEGA)
##################################################################
# Violin plot of Shannon's entropy
# A violin plot is a boxplot combined with a kernel density estimate of the probability density function per point.
# http://stanford.edu/~mwaskom/software/seaborn/tutorial/plotting_distributions.html
# For future check http://stackoverflow.com/questions/27322650/how-to-use-seaborn-pointplot-and-violinplot-in-the-same-figure-change-xticks-a
##################################################################

# Create a dictonary of the Shannon's entropy list for each MSA
#d = dict( CLUSTAL_OMEGA = CLUSTAL_OMEGA, MAFFT = MAFFT, MUSCLE = MUSCLE, GOLD=GOLD)
d = dict( GOLD_AA=GOLD_AA, GOLD_CO=GOLD_CO, GOLD_MAFFT =GOLD_MAFFT,GOLD_MUSCLE=GOLD_MUSCLE, GOLD_KALIGN=GOLD_KALIGN)
# Create a pandas data frame of the dictionary above
data = pd.DataFrame(dict([ (k,pd.Series(v)) for k,v in d.iteritems() ]))


# Plot the data in dictionary
fig = plt.figure()
sns.set(style="whitegrid")
plot_vals = [v.dropna() for k, v in data.iteritems()] 
sns.violinplot(plot_vals, names=data.columns, color="coolwarm_r", lw=3);
fig.suptitle("Distribution of Shannon's entropy", fontsize=16)
plt.xlabel('Alignment tool', fontsize=16)
plt.ylabel("Shannon's entropy", fontsize=16)
plt.show()


# Time to run the code: end timer
t1 = time.time()

total = t1-t0
print (total)

