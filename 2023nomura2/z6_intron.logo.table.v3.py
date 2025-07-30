#!/usr/bin/env python3

"""
Generate SVG-formatted SeqLogo from intron list
(https://github.com/betteridiot/seqlogo)

pre-requisition: pdf2svg for [seqlogo] SVG output
(https://github.com/dawbarton/pdf2svg)

appearing range (pandas table size) adjustment implemented
======================================
21.0901 implement "seqLogo" directly
21.08.31 Junesk9
"""

import sys
import numpy as np
import pandas as pd
import seqlogo

#all n=703,976
#int_f1 = "0_input/Eag.draft-genome.20180531.gene.id_edited.cds_edited.gff3.intron.filtered"
#filtered n=86,939
#int_f2 = "0_input/Eag.draft-genome.20180531.gene.id_edited.cds_edited.gff3.intron.filtered.ambiguity_removed.2.noshift.new"

#int_f1 = "../2_extract.intron/Eag_consensusScaffold.1kb.rename.gff3.intron.filtered"
#int_f2 = "Eag_consensusScaffold.1kb.rename.gff3.allclear.intron"

############## INPUT
#int_f1 = "../1_extract.intron/Eag.ref.gff3.intron.filtered"
#int_f2 = "Eag.ref.gff3.allclear.intron"
#int_f3 = "5_Eug.intron.excp-CAGCTG.txt"
#int_f4 = "6_Eug.intron.excp-CAGCTG.slided.txt"
#inf = sys.argv[-1]


#inf = "2_1kb.rename.final.intron.list.add.context.v3.tsv"
#inf = "2_Eag.ref_final.intron.list.add.context.v2.tsv"
#ints = [i.strip().split("\t") for i in open(int_f2)]
#ints = [i.strip().split("\t") for i in open(int_f1)]
#ints = [i.strip().split("\t") for i in open(int_f3)]
#ints = [i.strip().split("\t") for i in open(int_f4)]

inf = "Both.intron_short.fasta"
inf = "Both.intron_overhang.fasta"
ints = [i.strip() for i in open(inf) if not i.startswith(">")]

conv5, conv3 = {}, {}
nc5, nc3 = {},{}
all5, all3 = {}, {}
cagnnn5, cagnnn3 = {}, {}
nnnctg5, nnnctg3 = {}, {}
unknown5, unknown3 = {}, {}

############## Functions
def NTcount(seq, dic):
    for en,s in enumerate(seq):
        if s != "N":
            if not s in dic.keys():
                dic[s] = [0 for i in range(len(seq))]
            else: pass
            dic[s][en] += 1
        else: pass

    return dic

############### Running body
for i in ints:
    #cls = i[1]
    #nc_cls = i[8]
    seq = i.split("nnnnnnnnnn")
    #print(seq)
    #print(nc_cls)
    #sys.exit()
    seq5, seq3 = seq[0], seq[1]

    NTcount(seq5, all5)
    NTcount(seq3, all3)




for en, x in enumerate([all5, all3]):
    outname = ["SeqLogo.BOTH-intron_oh.5prime.svg","SeqLogo.BOTH-intron_oh.3prime.svg"]
    print(x.keys())

    if len(x.keys()) < 4:
        pass
    else:
        print(x)
        #sys.exit()
        outf = outname[en]
        df = pd.DataFrame.from_dict(x, orient="index")
        df = df.sort_index()
        print(en)
        print(df)
        df = df/df.sum() ## convert count to proportion
        df = df.transpose()
        ppm = seqlogo.Ppm(df)
        seqlogo.seqlogo(ppm, ic_scale = True, format = 'svg', size = 'medium', filename = outf)

        ##5-nt reduction to express
        outf = "short."+outf
        df = df.iloc[:-5,] if en%2 == 0 else df.iloc[5:,] ##different cut for 5-prime & 3-prime
        df.index = [i for i in range(len(df.index))] ## new rownames (0-15) are required.
        ppm = seqlogo.Ppm(df)
        seqlogo.seqlogo(ppm, ic_scale = True, format = 'svg', size = 'medium', filename = outf)

        #flushing memory, for sure
        df, ppm = 0,0


