#!/usr/bin/env python3


import sys
import pandas as pd

inf = "Eg-Chen.Isoseq-NGS-merged.genome.verA1.gff.intron.filtered.ambiguity_removed.2.shift.check.nonslide"
inf = "Eg-Chen.Isoseq-NGS-merged.genome.verA1.gff.intron.filtered"
inf = sys.argv[-1]
outf = inf + ".cagctg"
outf2 = inf + ".nnnnnn"

inf = [i.strip().split("\t") for i in open(inf)][1:]


cagctg = []
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

    if cag3 != "CAG":
        cntMotif[cls][2] += 1
    else:
        if ctg5 == "CTG":
            cntMotif[cls][0] += 1
            cagctg.append(i)
        elif ctg6 == "CTG":
            cntMotif[cls][1] += 1
            cagctg.append(i)
        else: 
            cntMotif[cls][2] += 1

cv_keys = ["GT-AG","GC-AG","else"]
nc_keys = ['CAG-CTG', 'CGG-CTG', 'CTG-CAG', 'CTG-CGG', 'CGG-CCG', 'CAG-nnn', 'nnn-CTG', 'CTG-nnn', 'CGG-nnn', 'nnn-CAG', 'nnn-CCG', 'nnn-CGG', 'nnn-nnn']
cnt_df = pd.DataFrame(0, index=nc_keys, columns=cv_keys)
ry_df = pd.DataFrame(0, index=nc_keys, columns=cv_keys)

def NCer(cag3, ctg5, ctg6):
    key5s = set([i.split("-")[0] for i in nc_keys])
    key3 = set([i.split("-")[1] for i in nc_keys])

    cag3 = cag3.upper()
    ctg5 = ctg5.upper()
    ctg6 = ctg6.upper()

    #key5 = "nnn" if not cag3 in key5 else cag3
    key5 = cag3 if cag3 in key5s else "nnn"
    #nc_key = True
    #print(cag3, key5, ctg5, ctg6)
    if key5 == "CAG":
        if "CTG" in [ctg5, ctg6]:
            nc_key = "CAG-CTG"
        else:
            nc_key = "CAG-nnn"
    elif key5 == "CGG":
        if "CTG" in [ctg5, ctg6]:
            nc_key = "CGG-CTG"
        elif "CCG" in [ctg5, ctg6]:
            nc_key = "CGG-CCG"
        else:
            nc_key = "CGG-nnn"
    elif key5 == "CTG":
        if "CAG" in [ctg5, ctg6]:
            nc_key = "CTG-CAG"
        elif "CGG" in [ctg5, ctg6]:
            nc_key = "CTG-CGG"
        else:
            nc_key = "CTG-nnn"
    elif key5 == "nnn":
        if "CAG" in [ctg5, ctg6]:
            nc_key = 'nnn-CAG'
        elif "CCG" in [ctg5, ctg6]:
            nc_key = 'nnn-CCG'
        elif "CGG" in [ctg5, ctg6]:
            nc_key = 'nnn-CGG'
        elif "CTG" in [ctg5, ctg6]:
            nc_key = "nnn-CTG"
        else:
            nc_key = "nnn-nnn"
    else:
        nc_key = key5+"-"+key3
        print("[ERROR] unclassified: %s" % nc_key)
        sys.exit()

    return nc_key

out2 = open(outf2, "w")
nlist = open("nnnnnn.list","w")
for i in inf:
    cv_key = i[2]
    cv_key = "else" if not cv_key in cv_keys[:2] else cv_key
    #print(cv)

    seq  = i[3]
    seq = seq.split(" ")
    seq5 = seq[0]
    seq3 = seq[2]

    cag3 = seq5[14:17]
    ctg5 = seq3[2:5]
    ctg6 = seq3[1:4]
    ry = seq3[10].upper()
    #ry = "R" if ry in ["A","G"] else "Y"

    nc_key = NCer(cag3, ctg5, ctg6)
    cnt_df.at[nc_key, cv_key] += 1

    if cv_key == "else" and nc_key == "nnn-nnn":
        iid = i[0]
        print(iid, file=nlist)
        print(*i, sep="\t", file=out2)

    if ry in ["A","G"]:
        ry_df.at[nc_key, cv_key] += 1
    else:
        pass
    
nlist.close()
out2.close()
print(cnt_df)
print()
print()
print("[puRine] at the 3'-border")
print(ry_df)


print()
print("Count 3CAG-CTG5, 3CAG-CTG6, else")
print(cntMotif)
out = open(outf, "w")
for o in cagctg:
    print(*o, sep="\t", file=out)
out.close()


