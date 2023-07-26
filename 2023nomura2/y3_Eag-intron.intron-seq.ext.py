#!/usr/bin/env python3

import sys
import glob

nons = "Non-Slid.intron.id.list"
slid = "Slidable.intron.id.list"

nonc = "Nonc.intron.id.list"
conv = "Conv.intron.id.list"

ref_fa = "Eag.draft-genome.20160925.dna.fa"


len_limit = [1500,2000]
seq_no = 10000


def RevComp(oligo):
    rev_dict = {"A":"T","T":"A","G":"C","C":"G","N":"N"}

    oligo = oligo.upper()
    rev = []
    for o in oligo:
        r = rev_dict[o]
        rev.append(r)

    rc = rev[::-1]
    rc = "".join(rc)

    return rc

tmp = {}
for r in open(ref_fa):
    r = r.strip()
    if r.startswith(">"):
        name = r[1:].split("_")[0]
        tmp[name] = []
    else:
        tmp[name].append(r)

for n in list(tmp.keys()):
    tmp[n] = "".join(tmp[n])
ref = tmp

out1_f = "Nonc.intron.seq.u100.fa"
out2_f = "Conv.intron.seq.u100.fa"

def ExtFa(ls, out_f):
    out = open(out_f, "w")
    for i in open(ls):
        int_id = i.strip()
        [ch, rng, ori] = int_id.split(".")
        [st, ed] = rng.split("-")
        st, ed = int(st), int(ed)
        int_len = ed - st + 1

        #if len_limit[0] <= int_len <= len_limit[1]:
        if int_len < 100:
            int_seq = ref[ch][st-1:ed]
            int_seq = RevComp(int_seq) if ori == "-" else int_seq

            name = ">" + int_id
            print(name, file=out)
            print(int_seq, file=out)
        else: pass
    out.close()

ExtFa(nonc, out1_f)
ExtFa(conv, out2_f)

