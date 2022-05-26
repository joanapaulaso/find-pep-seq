# import easygui
import os
from collections import Counter
from math import sqrt
import csv

threshold = 0
seq_min_diff = 0
apt = 0
directory = 'seqs_insilico'
included = []

for filename in os.scandir(directory):

    if filename.is_file():
        insilicoseq_file = open(filename, "r")
        insilicoseq_filter = insilicoseq_file.read()
        insilicoseq = insilicoseq_filter.split("\n")
        insilicoseq_file.close()

        # msseq_file = open("id_sequences_ms.txt", "r") # change the filename here
        msseq_file = open("id_sequences_ms_higher_3kDa.txt", "r") # change the filename here
        msseq_filter = msseq_file.read()
        msseq = msseq_filter.split("\n")
        msseq_file.close()

        def seq2vec(seq):
            
            cseq = Counter(seq) # count the aminoacids in sequence
            diffseq = set(cseq) # set of the different aminoacids in sequence
            seqlength = sqrt(sum(c * c for c in cseq.values())) # length of the seq vector

            return cseq, diffseq, seqlength

        
        def cosdis(seq1, seq2): # function that check intersections and calculates cosine
            
            intersect = seq1[1].intersection(seq2[1] ) # which amino acids are common to the two seqs?

            return sum(seq1[0][ch] * seq2[0][ch] for ch in intersect) / (seq1[2] * seq2[2])

        def seqs_len(seq1, seq2): # sequences length difference normalization

            seq1_len = len(seq1) / 100
            seq2_len = len(seq2) / 100
            seq_diff = seq1_len - seq2_len
            
            if seq_diff < 0:
                seq_diff = seq_diff * -1
            
            return seq_diff

        threshold = 0.9 # cosine threshold
        seq_min_diff = 0.08 # Max seq length difference of amino acids

        for insilicoseqline in insilicoseq:
            for msseqline in msseq:
                try:
                    cosine = cosdis(seq2vec(msseqline), seq2vec(insilicoseqline))
                    len_calc = seqs_len(msseqline, insilicoseqline)
                    if cosine > threshold and len_calc <= seq_min_diff:
                        filename = str(filename).replace("<DirEntry 'dig_", "")
                        filename = str(filename).replace("_seqs.fasta'>", "")
                        rslt = [msseqline, insilicoseqline, cosine, len_calc, filename]
                        included.append(rslt)
                        apt = 1
                    else:
                        apt = 0
                except IndexError:
                    pass

header = ["MS seq", "in silico seq", "cosine value", "seq difference", "protein", "Cosine threshold = {}".format(threshold), "Seq max diff = {}".format(seq_min_diff)]

# with open('seq_results_higher_3kDa.csv', 'w', encoding='UTF8', newline='') as f:
with open('seq_results_higher_3kDa.csv', 'w', encoding='UTF8', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(header)
    writer.writerows(included)