#!/usr/bin/env python3
"""
The intronic features was supplemented by agat 1.4.2 by singularity exec ~/data/utils/AGAT/agat.sif agat_sp_add_introns.pl -f Eg-Chen.Isoseq-NGS-merged.genome.verA1.gff
"""

import sys
import pandas as pd

inf = "Eg-Chen.Isoseq-NGS-merged.genome.verA1.gff.intron.filtered.ambiguity_removed.2.shift.check"
gff = "Eg-Chen.Isoseq-NGS-merged.genome.verA1.wINT.gff"


print("[PROGRESS] start parse GFF introns")
gff = [i.strip("\t").split() for i in open(gff)]
gff_dic = {}
inf_dic = {}
mids = []
x, y = 0, 0
for g in gff:
    if g[0].startswith("#"):
        pass
    elif g[1] == "Isoseq_tama" and g[2] == "intron":
        y += 1
        ch = g[0]
        st = g[3]
        ed = g[4]
        ori = g[6]
        int_id = ch + "." + st + "-" + ed + "." + ori 
        tc_id = g[-1].split("=")[-1]
        gff_dic[int_id] = tc_id
    elif g[1] == "Isoseq_tama" and g[2] == "mRNA":
        gid = g[8].split(";")[0].split("=")[-1]
        st = int(g[3])
        ed = int(g[4])
        inf_dic[gid] = [0,0,0,ed-st+1] # conv, nc-cag, nc-unknown, g_len
        mids.append(gid)
    else: pass
    x += 1
print("[PROGRESS] done parse GFF introns: %s/%s obtained" % (y,x))

print()
print("[PROGRESS] start add tc_id to the intron file")
inf = [i.strip().split() for i in open(inf)]
outf = "Eg-Chen.Isoseq-NGS-merged.genome.verA1.gff.intron.filtered.ambiguity_removed.2.shift.check.wGID"
out = open(outf, "w")

x, y = 0, 0
for i in inf:
    pnt = True
    if i[0].startswith("#Int"):
        new_i =  ["#GID"] + i
    else:
        int_id = i[0]
        #print(int_id)
        #print(list(gff_dic.keys())[:5])
        if int_id in gff_dic.keys():
            new_i = [gff_dic[int_id]] + i
            y += 1
        else: 
            pnt = False
        x += 1
    if pnt:
        print(*new_i, sep="\t", file=out)
out.close()
print("[PROGRESS] done add tc_id: %s/%s added/retrieved" % (y,x))

print()
print("[PROGRESS] parse GID/introns")
inf = [i.strip().split()[:4] for i in open(outf)][1:]
nlist = set([i.strip() for i in open("nnnnnn.list")])
#inf_dic = {m:[0,0] for m in mids}
for i in inf:
    gid = i[0]
    iid = i[1]
    #if not gid in inf_dic.keys():
        #inf_dic[gid] = [0,0] #conv, nonc
    #else: pass

    if i[2].startswith("Non-"):
        if iid in nlist:
            inf_dic[gid][2] += 1
        else:
            inf_dic[gid][1] += 1
    elif i[2].startswith("Conv"):
        inf_dic[gid][0] += 1
    else: pass

df = pd.DataFrame.from_dict(inf_dic)
df = df.T
df.columns = ["conv","nc_cagctg","nc_uncat","gene_len"]

print(df)
df.to_csv("Eg-Chen.Isoseq-NGS-merged.genome.verA1.gff.intron.filtered.ambiguity_removed.2.shift.check.wGID.gCount", sep="\t")


