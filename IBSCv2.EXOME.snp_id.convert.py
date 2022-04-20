#!/usr/bin/env python3
"""
Place Exome SNPs to the barley genes

2021.04.10
by June-Sik Kim (https://github.com/junesk9/)
"""

import sys
import glob
import time

stime = time.time()

snp_fi = "SNP.pos.info.tsv"
snps = [i.strip().split() for i in open(snp_fi).readlines()[1:]]
#print(snps[:5])
#sys.exit()

gene_info = "IBSCv2.gene.sted.txt"
gene_info = [i.strip().split() for i in open(gene_info)][1:]


out_f = "IBSCv2.EXOME.snp_id.convert.tsv"
out = open(out_f,"w")
y, z = 0, 0
for s in snps:
    y += 1
    snp = s[0]
    s_ch = s[1]
    s_pos = int(s[2])
    s_found = False
    while s_found == False:
        for g in gene_info:
            g_id = g[0]
            g_ch = g[1][3:]
            g_st = int(g[2])
            g_ed = int(g[3])
            g_rng = list(range(g_st, g_ed+1))
            if g_ch == s_ch and s_found == False:
                if s_pos in g_rng:
                    new_id = g_id+"."+str(s_pos)
                    print(snp, g_ch, new_id, sep="\t", file=out)
                    z += 1
                    if z%10 == 0:
                        rtime = int(time.time()-stime)
                        rtime = "%sm %ss" % (int(rtime/60), rtime%60)
                        print(z, y, rtime)
                    s_found = True
                else:
                    pass
            else: pass
        s_found = True
out.close()
