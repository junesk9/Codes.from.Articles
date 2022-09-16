#!/usr/bin/env python3
"""
*Search and count CAMTA motives (A/C/G)CGCG(T/C/G)
in the 1-kb upstream region for each gene of Brachypodium genome

*The output is a CSV table of relative position 
of each found motif from the gene start (semicolon separated)
and total count for rows and columns.

*The reference genome is adopt from Ensembl Plants DB.
(https://plants.ensembl.org/index.html)

2022.05.28
by June-Sik Kim (https://github.com/junesk9/)
"""

import sys, glob
import pandas as pd #v1.4.1

#################################################################
######################## Input Argument
ref = "./0_ref/Bdistachyon_314_v3.0.fa"
gff = "./0_ref/Bdistachyon_314_v3.1.gene.gff3"

prom_len = 1000

header = ["ACGCGT","ACGCGC","ACGCGG","CCGCGT","CCGCGC","CCGCGG","GCGCGT","GCGCGC","GCGCGG","Count"]
################################################################
######################## Function
def RevComp(seq):
    rev_d = {"A":"T","T":"A","C":"G","G":"C","N":"N"}

    seq = seq.upper()
    seq = [i for i in seq]
    comp = []
    for s in seq:
        try:
            c = rev_d[s]
        except KeyError:
            c = "N"
        comp.append(c)
    revcomp = comp[::-1]
    revcomp = "".join(revcomp)

    return revcomp


#################################################################
####################### Input parsing

ref = [i.strip() for i in open(ref)]
tmp = {}
for r in ref:
    if r.startswith(">"):
        ch = r[1:] #remove ">"
        tmp[ch] = []
    else:
        seq = r
        tmp[ch].append(seq)
ref = {i:"".join(tmp[i]) for i in tmp.keys()}

gff = [g.strip().split() for g in open(gff) if not g.startswith("#")]
tmp = {}
for g in gff:
    #print(g)
    cat = g[2]
    if cat == "gene":
        ch = g[0]
        st, ed = int(g[3]), int(g[4])
        ori = g[6]
        gid = g[-1].split(".")[0].split("=")[-1]
        #print(ch,st,ed,ori,gid)
        seq = ref[ch][st:ed+1]
        seq = RevComp(seq) if ori == "-" else seq
        #print(seq)
        if ori == "+":
            p_st = st - prom_len
            p_ed = st
        elif ori == "-":
            p_st = ed + 1
            p_ed = ed + 1 + prom_len

        p_seq = ref[ch][p_st:p_ed]
        p_seq = RevComp(p_seq) if ori == "-" else p_seq
        
        tmp[gid] = p_seq

        ## validate the parsing
        if False and gid == "Bradi1g00213":
            print(ch,st,ed,ori,gid)
            print(seq)
            print(p_seq)
            print(len(p_seq))
            sys.exit()
gff = tmp
tmp = []
################################################################
###################### Count & output result
out_dic = {k:["NA" for i in header] for k in gff.keys()}
out_dic["Count"] = [0 for i in header]
for gid, p_seq in gff.items():
    x = 0
    for i in range(prom_len-6):
        motif = p_seq[i:i+6]
        #print(motif)
        if motif in header:
            x += 1
            m_idx = header.index(motif)
            pos = str(i - prom_len)

            out_dic["Count"][m_idx] += 1
            if out_dic[gid][m_idx] == "NA":
                out_dic[gid][m_idx] = pos
            else:
                out_dic[gid][m_idx] += ";"+pos
        else: pass
    out_dic[gid][-1] = x ##count per gene
    out_dic["Count"][-1] += x ##sum of count per gene

df = pd.DataFrame.from_dict(out_dic)
df = df.transpose()
df.columns = header
print(df)
df.to_csv("Bradi.CAMTA-motif.PosByGene.v2.csv")

