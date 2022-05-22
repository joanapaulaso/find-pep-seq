import easygui
from collections import Counter
from math import sqrt

apt = 0

insilicoseq_filename = easygui.fileopenbox()
insilicoseq_file = open(insilicoseq_filename, "r")
insilicoseq_filter = insilicoseq_file.read()
insilicoseq = insilicoseq_filter.split("\n")
# print("MS sequences = ")
# print(insilicoseq)
print(insilicoseq_filename)

insilicoseq_file.close()

msseq_file = open("id_sequences_ms.txt", "r") # change the filename here
msseq_filter = msseq_file.read()
msseq = msseq_filter.split("\n")
# print("in silico sequences = ")
# print(msseq)

msseq_file.close()

# seqcount = 0

# for seq in insilicoseq:
#     for silico in msseq:
#         if silico in seq:
#             seqcount += 1
#             print(silico)

# print(seqcount)

def seq2vec(seq):
    
    cseq = Counter(seq) # count the aminoacids in sequence
    diffseq = set(cseq) # set of the different aminoacids in sequence
    seqlength = sqrt(sum(c*c for c in cseq.values())) # length of the word vector

    # if apt == 1:
    #     print(lw)

    return cseq, diffseq, seqlength

# function that check intersection and calculates cosine
def cosdis(seq1, seq2):
    # which characters are common to the two words?
    intersect = seq1[1].intersection(seq2[1])
    # by definition of cosine distance we have

    if apt == 1:
        print(intersect)

    return sum(seq1[0][ch]*seq2[0][ch] for ch in intersect)/seq1[2]/seq2[2]

threshold = 0.95
for insilicoseqline in insilicoseq:
    for msseqline in msseq:
        try:
            cosine = cosdis(seq2vec(msseqline), seq2vec(insilicoseqline))
            if cosine > threshold:
                print("{},{},{}".format(msseqline, insilicoseqline, cosine))
                apt = 1
            else:
                apt = 0
        except IndexError:
            pass