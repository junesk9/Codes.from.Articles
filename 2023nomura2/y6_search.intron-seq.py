#!/usr/env python3



import sys



def ProcessFa(fa):
    fa = [i.strip() for i in open(fa)]
    dic = {}
    for f in fa:
        if f.startswith(">"):
            name = f[1:]
            dic[name] = []
        else:
            s = f
            dic[name].append(s)
    dic = {k:"".join(v) for k,v in dic.items()}

    return dic

def RevComp(seq):
    dic = {"A":"T", "T":"A", "G":"C", "C":"G"}

    seq = [i.upper() for i in seq]
    comp = [dic[i] for i in seq if i in dic.keys()]
    rc = comp[::-1]
    rc = "".join(rc)

    return rc


target = "search.target.fa"
target = ProcessFa(target)

total = "Eag.draft-genome.20160925.Nonc.intron.seq.fa"
total = "Eag.draft-genome.20160925.Nonc.intron.seq.u100.fa"
total = ProcessFa(total)


print("[PROCESS] Parsing FASTA DONE")
print("")

targets = []
for k,v in target.items():
    rc = RevComp(v)
    for k2, v2 in total.items():
        if v == v2 or RevComp(v) == v2 or v in v2 or v2 in v or rc in v2 or v2 in rc:
            print(k, k2)
            targets.append(k2)
            print(v)
            print(v2)
        else: pass
    #print("[PROCESS] %s search done" % k)


t_chrs = [i.split(".")[0] for i in targets]

gff = "Eag.draft-genome.20180531.gene.id_edited.cds_edited.gt.gff3"
for g in open(gff):
    g = g.strip().split("\t")
    ch = g[0]
    cat = g[2]
    st, ed = int(g[3]), int(g[4])
    if ch in t_chrs:
        if cat == "gene":
            gid = g[-1].split(";")[0].split("=")[-1]
            for t in targets:
                t_line = t.split(".")
                t_ch = t_line[0]
                t_rng = t_line[1].split("-")
                t_st, t_ed = int(t_rng[0]), int(t_rng[1])

                if t_ch == ch and st < t_st and t_ed < ed:
                    print(t, gid)
                else: pass
        else: pass
    else: pass

    


