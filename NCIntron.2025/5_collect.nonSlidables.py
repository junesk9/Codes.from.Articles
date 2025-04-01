#!/usr/bin/env python3


import sys



inf = "Euglena.gene.change.chr.gff.intron.filtered.ambiguity_removed.2.shift.check"
outf = inf + ".nonslide"


out = open(outf, "w")
inf = [i.strip() for i in open(inf)]
x, y, z = 0,0,0
for en,i in enumerate(inf):
    i = i.split("\t") 
    if en == 0:
        print(*i, sep="\t", file=out)
    else:
        x += 1
        cls = i[1]
        slide = int(i[5])
        if slide == 0:
            if cls.startswith("Non"):
                z += 1
            elif cls.startswith("Conv"):
                y += 1
            else: pass
            print(*i, sep="\t", file=out)
out.close()
print("[PROGRESS] Non-slidable %s-conv %s-nonc introns chosen from %s total " % (y,z,x))


