#!/usr/bin/env python3

import sys


inp_f = "Eag.draft-genome.20180531.gene.id_edited.cds_edited.gt.gff3.intron_boundary.slided"
inp_ls = [i.strip() for i in open(inp_f)]

out_f = "Eag.draft-genome.20180531.gene.id_edited.cds_edited.gt.gff3.intron_boundary.slided.class"
out = open(out_f, "w")

for i in inp_ls:
    if i.startswith("Sca"):
        print(i, file=out)
    elif i.startswith("/"):
        print(i, file=out)
    else:
        line = i.split("\t")
        #print(line)
        #['0', 'GGCACGATCCG/AGCACGTCGTCGGCTCGTGG .. TGGCCGAGAACGTGCCGCCG/CCGCGCCCGCC', 'AG CG', '_']
        #sys.exit()

        conv, nonc = "NA", "unknown"
        border = line[1].split()[-1].split("/")
        border = border[0][-1] + border[1][0] ###3-end of intron and the next basepare i.e. "GC"
        #print(border)
        #sys.exit()

        ##classify conventional intron
        c_border = "-".join(line[2].split())
        if c_border in ["GT-AG", "GC-AG"]:
            conv = c_border
        else: pass

        ##classify nonconventional intron
        seq5 = line[1].split()[0].split("/")[-1]
        seq3 = line[1].split()[-1].split("/")[0]
        #print(seq5, seq3)
        #AGCACGTCGTCGGCTCGTGG TGGCCGAGAACGTGCCGCCG

        accP3 = seq5[3:6]
        dnnM5 = seq3[-8:-5]
        dnnM6 = seq3[-9:-6]
        if accP3 == "CAG":
            if dnnM5 == "CTG":
                nonc = "CAG-CTG-5"
            elif dnnM6 == "CTG":
                nonc = "CAG-CTG-6"
            elif "N" in dnnM5 + dnnM6:
                nonc = "missing"
            else:
                nonc = "CAG-nnn"
        elif dnnM5 == "CTG":
            if accP3 == "CGG":
                nonc = "CGG-CTG-5"
            elif "N" in accP3:
                nonc = "missing"
            else:
                nonc = "nnn-CTG-5"
        elif dnnM6 == "CTG":
            if accP3 == "CGG":
                nonc = "CGG-CTG-6"
            elif "N" in accP3:
                nonc = "missing"
            else:
                nonc = "nnn-CTG-6"

        elif accP3 == "CGG":
            if dnnM5 == "CCG":
                nonc = "CGG-CCG-5"
            elif dnnM6 == "CCG":
                nonc = "CGG-CCG-6"
        elif accP3 == "CTG":
            if dnnM5 == "CAG":
                nonc = "CTG-CAG-5"
            elif dnnM6 == "CAG":
                nonc = "CTG-CAG-6"
            elif dnnM5 == "CGG":
                nonc = "CTG-CGG-5"
            elif dnnM6 == "CGG":
                nonc = "CTG-CGG-6"
        elif accP3 == "CCG":
            if dnnM5 == "CGG":
                nonc = "CCG-CGG-5"
            elif dnnM6 == "CGG":
                nonc = "CCG-CGG-6"
        elif "N" in dnnM5 + dnnM6:
            nonc = "missing"
        #sys.exit()
        
        if line[-1] == "0": ##the original and slidable introns only
            pass
        else:
            #line[2] = c_border
            new_line = line + [border, conv, nonc]
            print(*new_line,sep="\t",file=out)

out.close()

