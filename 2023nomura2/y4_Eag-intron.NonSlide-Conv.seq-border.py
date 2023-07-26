#!/usr/bin/env python3


import sys


nons = "Non-Slid.intron.id.list"
conv = "Conv.intron.id.list"

nons = [i.strip() for i in open(nons)]
conv = [i.strip() for i in open(conv)]
nons_conv = list(set(nons).intersection(set(conv)))
print("[progress] %s Non-slidable, %s Conv shares %s introns" % (len(nons), len(conv), len(nons_conv)))


genome_fa = "Eag.draft-genome.20160925.dna.fa"
tmp = {}
for i in open(genome_fa):
    i = i.strip()
    if i.startswith(">"):
        name = i[1:]
        tmp[name] = []
    else:
        s = i
        tmp[name].append(s)
ref = {k:"".join(v) for k,v in tmp.items()}


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



out1 = "Eag.draft-genome.20160925.Non-Slidable.Conv.Intron-border.fa"
out2 = "Eag.draft-genome.20160925.Non-Slidable.Nonc.Intron-border.fa"

out1 = open(out1, "w")
out2 = open(out2, "w")
for k in nons:
    ## extract intron-seq
    i = k.split(".")
    i = i + ["0"] if len(i) == 3 else i

    scf = i[0]
    st = int(i[1].split("-")[0])
    ed = int(i[1].split("-")[1])
    slide = int(i[3])
    st = st + abs(slide) if st > 0 else st - abs(slide)
    ed = ed + abs(slide) if st > 0 else ed - abs(slide)
    ori = i[2]

    scf_fa = ref[scf]
    int_fa = scf_fa[st-1:ed]
    anchor5 = scf_fa[st-4:st-1]
    anchor3 = scf_fa[ed:ed+3]
    int_fa = anchor5 + int_fa + anchor3
    int_fa = RevComp(int_fa) if ori == "-" else int_fa
    int_fa = int_fa[:3].lower() + int_fa[3:-3] + int_fa[-3:].lower()
    int_fa = int_fa[3:-3]
    int_fa = int_fa[:20] + "nnnnnnnnnn" + int_fa[-20:]

    int_id = ">" + k
    if "N" in int_fa[:20] + int_fa[-20:]:
        pass
    else:

        if k in conv:
            print(int_id, file=out1)
            print(int_fa, file=out1)
        else:
            print(int_id, file=out2)
            print(int_fa, file=out2)
out1.close()
out2.close()



