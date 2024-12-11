#!/usr/bin/env python3


import sys



inf = "Eg-Chen.Isoseq-NGS-merged.genome.gff.intron.filtered.ambiguity_removed.2.shift.check.nonslide"
inf = "Euglena.gene.change.chr.gff.intron.filtered.ambiguity_removed.2.shift.check.nonslide"
inf = [i.strip().split("\t") for i in open(inf)][1:]



cntMotif = {"NC":[0,0,0], "CV":[0,0,0]} #3CAG-CTG5, 3CAG-CTG6, else
for i in inf:
    cls = i[1]
    cls = "NC" if cls.startswith("Non") else "CV"

    seq = i[3]
    seq = seq.split(" ")
    seq5 = seq[0]
    seq3 = seq[2]

    cag3 = seq5[14:17]
    ctg5 = seq3[2:5]
    ctg6 = seq3[1:4]
    #print(cag3, ctg5, ctg6)

    if cag3 != "CAG":
        cntMotif[cls][2] += 1
    else:
        if ctg5 == "CTG":
            cntMotif[cls][0] += 1
        elif ctg6 == "CTG":
            cntMotif[cls][1] += 1
        else: cntMotif[cls][2] += 1

print(cntMotif)

