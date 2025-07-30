#!/usr/bin/env python3


import sys
import pandas as pd

inf = "Eag.draft-genome.20180531.gene.id_edited.cds_edited.gt.gff3.intron_boundary.slided.class2"
inf = [i.strip() for i in open(inf)]


mode_noN = True ### discard all "N"-included introns

tmp = {}
slide_count = {}
for i in inf:
    if i.startswith("Scaf"):
        name = i
        tmp[name] = []
    elif i != "//":
        #0	GGCACGATCCG/AGCACGTCGTCGGCTCGTGG .. TGGCCGAGAACGTGCCGCCG/CCGCGCCCGCC	AG CG	_	GC	NA	unknown
        i = i.split("\t")
        slide = i[0]
        line = i[-3:] + [slide]
        seq = i[1]
        if mode_noN and "N" in seq:
            pass
        else:
            tmp[name].append(line)
            ###count slides
            if slide in slide_count.keys():
                slide_count[slide] += 1
            else:
                slide_count[slide] = 1

slidable, nonslidable =  {}, {}
sld_fi, nons_fi = "Slidable.intron.id.list", "Non-Slid.intron.id.list"
conv_fi, nonc_fi = "Conv.intron.id.list", "Nonc.intron.id.list"
both_fi, both_ns_fi = "Both.intron.id.list", "Both_ns.intron.id.list"
out1 = open(sld_fi, "w")
out2 = open(nons_fi, "w")
out3 = open(conv_fi, "w")
out4 = open(nonc_fi, "w")
out5 = open(both_fi, "w")
out6 = open(both_ns_fi, "w")

missing = 0
for k,v in tmp.items():
    if len(v) == 1:
        nonslidable[k] = v
        print(k, file=out2)
    elif len(v) > 1:
        slidable[k] = v
        print(k, file=out1)
    elif len(v) == 0:
        missing += 1
print()
print("Slidable: %s, Non-Slidable: %s, N-included (discard): %s ,Total: %s" % (len(slidable.keys()), len(nonslidable.keys()), missing, len(tmp.keys())))
out1.close()
out2.close()

print("### Slidable introns frequency by sliding pos.")
s_keys = sorted(list(slide_count.keys()))
for k in s_keys:
    if k != "0":
        print(k, slide_count[k], sep="\t")
    else:
        print("total", slide_count[k], sep="\t")
#end of for
print()

tmp = {**slidable, **nonslidable}
conv_keys = sorted(list(set([i[0][1] for i in tmp.values()]))) 
nonc_keys = sorted(list(set([i[0][2] for i in tmp.values()])))
border_keys = sorted(list(set([i[0][0] for i in tmp.values()])))
nonc_keys_short = sorted(list(set(["-".join(i.split("-")[:2]) for i in nonc_keys])))


#(original)
nonc_keys = ['CAG-CTG-5', 'CAG-CTG-6', 'CAG-nnn', 'CCG-CGG-5', 'CCG-CGG-6', 'CCG-nnn', 'CGG-CCG-5', 'CGG-CCG-6', 'CGG-CTG-5', 'CGG-CTG-6', 'CGG-nnn', 'CTG-CAG-5', 'CTG-CAG-6', 'CTG-CGG-5', 'CTG-CGG-6', 'CTG-nnn', 'nnn-CAG', 'nnn-CCG', 'nnn-CGG', 'nnn-CTG-5', 'nnn-CTG-6', 'unknown']
nonc_keys_short = ['CAG-CTG', 'CAG-nnn', 'CCG-CGG', 'CCG-nnn', 'CGG-CCG', 'CGG-CTG', 'CGG-nnn', 'CTG-CAG', 'CTG-CGG', 'CTG-nnn', 'nnn-CAG', 'nnn-CCG', 'nnn-CGG', 'nnn-CTG', 'unknown']

#(manually given)
nonc_keys = ['CAG-CTG-5','CAG-CTG-6', 'CGG-CTG-5', 'CGG-CTG-6', 'CTG-CAG-5','CTG-CAG-6', 'CTG-CGG-5','CTG-CGG-6', 'CGG-CCG-5','CGG-CCG-6', 'CAG-nnn', 'nnn-CTG-5','nnn-CTG-6', 'CTG-nnn', 'CGG-nnn', 'nnn-CAG', 'nnn-CCG', 'nnn-CGG', 'unknown']
nonc_keys_short = ['CAG-CTG', 'CGG-CTG', 'CTG-CAG', 'CTG-CGG', 'CGG-CCG', 'CAG-nnn', 'nnn-CTG', 'CTG-nnn', 'CGG-nnn', 'nnn-CAG', 'nnn-CCG', 'nnn-CGG', 'unknown']

## collecting intron ids here
nnn_keys = ['CAG-nnn', 'nnn-CTG', 'CTG-nnn', 'CGG-nnn', 'nnn-CAG', 'nnn-CCG', 'nnn-CGG']
nnn_dic = {i:[] for i in nnn_keys}

## for 2-exon start nt
cv_df = pd.DataFrame(0, index=conv_keys, columns=["A","G","T","C"])
nc_df = pd.DataFrame(0, index=nonc_keys_short, columns=["A","G","T","C"])


print("### Counting Non-Slidable introns ...")
ns_df1 = pd.DataFrame(0, index=nonc_keys, columns=conv_keys)
ns_df2 = pd.DataFrame(0, index=nonc_keys_short, columns=conv_keys)
ns_df3 = pd.DataFrame(0, index=border_keys, columns=conv_keys[:2]+["NonConv", "unknown"])

for k,v in nonslidable.items():
    [border, conv, nonc, slide] = v[0]
    nonc_short = "-".join(nonc.split("-")[:2])
    # if nonc_key was not matched...
    nonc = "unknown" if not nonc in nonc_keys else nonc
    nonc_short = "unknown" if not nonc_short in nonc_keys_short else nonc_short

    #nnn_dic generation
    if nonc_short in nnn_keys and conv == "NA":
        nnn_dic[nonc_short].append(k)
    else: pass

    #counting intron classes
    ns_df1.at[nonc, conv] += 1
    ns_df2.at[nonc_short, conv] += 1
    if conv == "NA": ##split-down "NonConv" to matching any of NC or not
        if nonc in ["unknown", "missing"]:
            conv = "unknown"
        else:
            conv = "NonConv"
        print(k, file=out4) #supplementary output
    else: 
        print(k, file=out3) ##supplementary output
    ns_df3.at[border, conv] += 1

    # for extract the BOTH intron ids
    if conv in ["GC-AG", "GT-AG"] and nonc_short in  ['CAG-CTG', 'CGG-CTG', 'CTG-CAG', 'CTG-CGG', 'CGG-CCG']:
        print(*[k, "nonslide", conv, nonc_short], sep="\t", file=out5)
        print(*[k, "nonslide", conv, nonc_short], sep="\t", file=out6)
    else: pass

    # the 2-exon start nt
    ex2_st = border[-1]
    if v[0][1] == "NA":
        nc_df.at[nonc_short, ex2_st] += 1
    else:
        cv_df.at[conv, ex2_st] += 1
print("### Non-slidable conv, non-conv")
print(ns_df1)
print()
print("### Non-slidable conv, non-conv-short")
print(ns_df2)
print()
print("### Non-slidable 3'-border dimer")
print(ns_df3)



print()
print()
print("### Counting Slidable introns ...")
print("### Priority order CAG-CTG > CXG-CX*G > CAG-nnn & nnn-CTG > other nnn >  missing (N-incl) > unknown")

ys_df1 = pd.DataFrame(0, index=nonc_keys, columns=conv_keys[:2])
ys_df2 = pd.DataFrame(0, index=nonc_keys_short, columns=conv_keys[:2])
ys_df3 = pd.DataFrame(0, index=nonc_keys_short, columns=conv_keys)



for k,v in slidable.items():
    border = [i[0] for i in v]
    conv = [i[1] for i in v]
    nonc = [i[2] for i in v]
    slide = [i[3] for i in v]


    if "GT-AG" in conv:
        conv_idx = conv.index("GT-AG")
        conv = conv[conv_idx]
        nonc_long = nonc[conv_idx]
        nonc_short = "-".join(nonc_long.split("-")[:2])
        # if nonc_key was not matched...
        nonc_long = "unknown" if not nonc_long in nonc_keys else nonc_long
        nonc_short = "unknown" if not nonc_short in nonc_keys_short else nonc_short

        ys_df1.at[nonc_long, conv] += 1
        ys_df2.at[nonc_short, conv] += 1
        ex2_st = border[conv_idx][-1] #the 2-exon start nt
        cv_df.at[conv, ex2_st] += 1
        print(k, file=out3) ##supplementary output

        # for extract the BOTH intron ids in the same frame
        if nonc_short in  ['CAG-CTG', 'CGG-CTG', 'CTG-CAG', 'CTG-CGG', 'CGG-CCG']:
            print(slide, conv_idx)
            slide = slide[conv_idx]
            int_id = k + "." + slide
            print(*[int_id, "slide", conv, nonc_short], sep="\t", file=out5)
        else: pass       


    elif "GC-AG" in conv:
        conv_idx = conv.index("GC-AG")
        conv = conv[conv_idx]
        nonc_long = nonc[conv_idx]
        nonc_short = "-".join(nonc_long.split("-")[:2])
        # if nonc_key was not matched...
        nonc_long = "unknown" if not nonc_long in nonc_keys else nonc_long
        nonc_short = "unknown" if not nonc_short in nonc_keys_short else nonc_short

        ys_df1.at[nonc_long, conv] += 1
        ys_df2.at[nonc_short, conv] += 1
        ex2_st = border[conv_idx][-1] #the 2-exon start nt
        cv_df.at[conv, ex2_st] += 1
        print(k, file=out3) ##supplementary output

        # for extract the BOTH intron ids in the same frame
        if nonc_short in  ['CAG-CTG', 'CGG-CTG', 'CTG-CAG', 'CTG-CGG', 'CGG-CCG']:
            slide = slide[conv_idx]
            int_id = k + "." + slide
            print(*[int_id, "slide", conv, nonc_short], sep="\t",file=out5)
        else: pass


    else:
        conv = "NA"
        print(k, file=out4) ##supplementary output



    ##nonc-short classifying
    nonc_full = [i for i in nonc]
    nonc = ["-".join(i.split("-")[:2]) for i in nonc]
    nonc = sorted(list(set(nonc)))
    if len(nonc) == 1:
        nonc_short = nonc[0]
    elif len(nonc) > 1:
        for nc in nonc_keys_short:
            if nc in nonc:
                nonc_short = nc
                break
            else:
                nonc_short = "unknown"
    else: pass
    #nnn_dic generation
    if nonc_short in nnn_keys and conv == "NA":
        nnn_dic[nonc_short].append(k)
    else: pass
    """    
    elif len(nonc) == 2 and "unknown" in nonc:
        nonc.remove("unknown")
        nonc_short = nonc[0]
    elif len(nonc) == 2 and "missing" in nonc:
        nonc.remove("missing")
        nonc_short = nonc[0]
    else:
        nonc.remove("unknown") if "unknown" in nonc else nonc
        nonc.remove("missing") if "missing" in nonc else nonc
        nonc.remove("CAG-nnn") if "CAG-nnn" in nonc and len(nonc) > 1 else nonc
        nonc_short = nonc[0]
    """
    ##set the 2-exon start nt
    n_idx = False
    for en,n in enumerate(nonc_full):
        n_short = "-".join(n.split("-")[:2])
        if n_short == nonc_short:
            n_idx = en
            ex2_st = border[n_idx][-1]
            n_short_found = True
        elif en == len(nonc_full) and n_idx == False:
            print("ERROR")
            print(nonc_full, n, n_short, nonc_short)
            print(border)
            sys.exit()
        else: pass

    ys_df3.at[nonc_short, conv] += 1
    if conv == "NA":
        nc_df.at[nonc_short, ex2_st] += 1
    else: pass


print("### Slidable conv, non-conv-short")
print(ys_df3)
print()
print("### Slidable conv, non-conv, in the same frame")
print(ys_df2)
out3.close()
out4.close()

print()
print()
print("### 2nd-exon start NT at CV introns")
print(cv_df)
print()
print("### 2nd-exon start NT at NC introns")
print(nc_df)
out5.close()
out6.close()

for k,v in nnn_dic.items():
    out_f = k + ".intron_id.txt"
    out = open(out_f, "w")
    print(*v, sep="\n", file=out)
    out.close()


"""
####################################################################
####### Supplemental part to count additional motives 23.09.22

print()
print("#### counting  CK'G-nnn, nnn-CH'G")

ckgnnn = 0
nnnchg = 0
unclass = 0
for i in inf:
    if i.startswith("Scaf"):
        name = i
        tmp[name] = []
    elif i != "//":
        #0      GGCACGATCCG/AGCACGTCGTCGGCTCGTGG .. TGGCCGAGAACGTGCCGCCG/CCGCGCCCGCC    AG CG   _       GC      NA      unknown
        i = i.split("\t")
        line = i[-3:]
        seq = i[1]
        slide = i[0]
        cat1 = i[-2]
        cat2 = i[-1]
        if "N" in seq:
            pass
        elif cat1 == "NA" and cat2 == "unknown" and name in nonslidable.keys():
            seq5p = seq.split()[0]
            seq3p = seq.split()[-1]
            seq5p = seq5p.split("/")[-1]
            seq3p = seq3p.split("/")[0]
            cag3 = seq5p[2:5]
            ctg5 = seq3p[-8:-5]
            ctg6 = seq3p[-9:-6]
            if cag3 == "CTG" and not "CAG" in [ctg5, ctg6]:
                ckgnnn += 1
            elif cag3 == "CGG" and not "CCG" in [ctg5, ctg6]:
                ckgnnn += 1
            elif cag3 != "CTG" and "CAG" in [ctg5, ctg6]:
                nnnchg += 1
            elif cag3 != "CAG" and "CTG" in [ctg5, ctg6]:
                nnnchg += 1
            elif cag3 != "CGG" and "CCG" in [ctg5, ctg6]:
                nnnchg += 1
            else: 
                unclass += 1
        else: pass
print("[CKG-nnn] %s" % ckgnnn)
print("[nnn-CHG] %s" % nnnchg)
print("[In totla] %s" % (ckgnnn + nnnchg + unclass))
"""



