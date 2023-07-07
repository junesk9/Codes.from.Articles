#!/usr/bin/env python3
"""
List the merged PE reads of overhangs matched to the given DNA sequences,
then find the perfect match of the doner sequence for genome editing trials
output counts of each case of Insertion or base editing

Junesk9 07.07.23
"""

import glob
import sys
import gzip


indel_doner = "GTGGGCAATCATGCTGACTTCCTCGATCTTTGTTGGAACAACAACTCTCATCATCACCATCATCACCACTAAGAATTCGATATCGGATCCACCACTGAGGGTTTTCATTCCTTTTGACCTTGGTGCCGTGGGTCCCTTTC"
be_doner = "GTGGGCAATCATGCTGACTTCCTCGATCTAAGCTTGAACAACAACTCTCAACCACTGAGGGTTTTCATTCCTTTTGACCTTGGTGCCGTGGGTCCCTTTC"


def revcomp(seq):
    seq = [i for i in seq]

    comp_dic = {"A":"T", "T":"A", "C":"G", "G":"C"}

    comp = [comp_dic[i] for i in seq]
    rc = comp[::-1]
    rc = "".join(rc)

    return rc

indel_rc = revcomp(indel_doner)
be_rc = revcomp(be_doner)


fi_ls = glob.glob("./K*/*/out.extendedFrags.fastq.gz")
fi_ls = sorted(fi_ls)

print("##INDEL")
for fi in fi_ls:
    sample = fi.split("/")[-2].split("_")[-1]
    with gzip.open(fi, mode="rt") as gf:
        data = gf.read()
    seqs = data.split("\n")

    seq_len = int((len(seqs)-1) / 4)
    indel, be = 0, 0

    for s in seqs:
        if indel_doner in s:
            indel += 1
        elif indel_rc in s:
            indel += 1
        elif be_doner in s:
            be += 1
        elif be_rc in s:
            be += 1
        else: pass
    print(sample, seq_len, indel, be)
 
print()
print("##base exchange")
fi_ls = glob.glob("./*/out.extendedFrags.fastq.gz")
fi_ls = sorted(fi_ls)


for fi in fi_ls:
    sample = fi.split("/")[-2].split("_")[-1]
    with gzip.open(fi, mode="rt") as gf:
        data = gf.read()
    seqs = data.split("\n")

    seq_len = int((len(seqs)-1) / 4)
    indel, be = 0, 0

    for s in seqs:
        if indel_doner in s:
            indel += 1
        elif indel_rc in s:
            indel += 1
        elif be_doner in s:
            be += 1
        elif be_rc in s:
            be += 1
        else: pass
    print(sample, seq_len, indel, be)
"""
##INDEL
Control-1 76787 0 0
Control-2 69050 33203 0
Control-3 79280 0 0
Control-4 81692 0 0
EgGSL2T2KI-1 70965 35322 0
EgGSL2T2KI-2 81603 0 0
EgGSL2T2KI-3 75774 33464 0
EgGSL2T2KI-4 69678 35380 0

##base exchange
Control-1 76787 0 0
Control-2 69050 33203 0
Control-3 79280 0 0
Control-4 81692 0 0
EgGSL2T2RW-1 79905 0 26635
EgGSL2T2RW-2 93403 0 43796
EgGSL2T2RW-3 77578 0 32231
EgGSL2T2RW-4 74104 0 36515
"""




