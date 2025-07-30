#!/usr/bin/env python3


"""
Extract intron sequences from the intron survey form like
"1_final.intron.list" and genome multi-fasta file

This code only refers the first column to identify the region to be eluted

21.08.30
Junesk9
==================Change log
21.09.10 Also accept the slidable intron id

"""


import sys
import glob

len_limit = [1500,2000]
seq_no = 100
mode_all = True

def RevComp(oligo):
    rev_dict = {"A":"T","T":"A","G":"C","C":"G","N":"N"}

    oligo = oligo.upper()
    rev = []
    for o in oligo:
        r = rev_dict[o]
        rev.append(r)

    rc = rev[::-1]
    rc = "".join(rc)

    return rc

def ScfNameConv(scf):
    id_cnt = int(scf[8:])
    new_id = "scaffold" + str(id_cnt)

    print(scf,id_cnt,new_id)
    #sys.exit()
    
    return new_id



#ref_fa = "Eag_consensusScaffold.1kb.rename.fa"
ref_fa = "Eag.draft-genome.20160925.dna.fa"
#ref_fa = "Eag.ref.fa"

intron_fi = "Eag_consensusScaffold.1kb.rename.gff3.allclear.intron"
intron_fi = "Eag.ref.gff3.allclear.intron"
intron_fi = "Both.intron.id.list"
#intron_fi = "Both_ns.intron.id.list"
#print(intron_ls[:10])
#sys.exit()
#inp_ls = glob.glob("3_introns*txt")
#inp_ls = [intron_fi]
intron_ls = [i.strip().split("\t")[0] for i in open(intron_fi)]
#print(intron_ls[:2])
print("A total of %s introns are listed to be extracted" % len(intron_ls))



tmp = {}
for r in open(ref_fa):
    r = r.strip()
    if r.startswith(">"):
        name = r[1:].split("_")[0]
        tmp[name] = []
    else:
        tmp[name].append(r)

for n in list(tmp.keys()):
    tmp[n] = "".join(tmp[n])
ref = tmp


out_f1 = intron_fi.split(".")[0]+".intron_long.fasta"
out_f2 = intron_fi.split(".")[0]+".intron_short.fasta"
out_f3 = intron_fi.split(".")[0]+".intron_overhang.fasta"
out1 = open(out_f1,"w")
out2 = open(out_f2,"w")
out3 = open(out_f3,"w")

x, y = 0, 0
for i in intron_ls:
    int_id = ">" + i
    i = i.split(".")
    i = i + ["0"] if len(i) == 3 else i

    scf = i[0]
    st = int(i[1].split("-")[0])
    ed = int(i[1].split("-")[1])
    ori = i[2]
    slide = int(i[3])

    scf_fa = ref[scf]
    int_fa = scf_fa[st-6+slide:ed+5+slide]
    int_fa = RevComp(int_fa) if ori == "-" else int_fa
    int_fa = int_fa[:15] + "nnnnnnnnnn" + int_fa[-15:]
    print(int_id, file=out3)
    print(int_fa, file=out3)

    int_fa = scf_fa[st-1+slide:ed+slide]
    int_fa = RevComp(int_fa) if ori == "-" else int_fa
    print(int_id, file=out1)
    print(int_fa, file=out1)

    int_fa = int_fa[:20] + "nnnnnnnnnn" + int_fa[-20:]
    print(int_id, file=out2)
    print(int_fa, file=out2)
out1.close()
out2.close()
out3.close()





#ref = {tmp[i]:"".join(tmp[i]) for i in list(tmp.keys())}
#int_ls = "2_introns.CGGasACPTR.list.txt"
#int_ls = "3_introns.CGG-CTG.list.txt"
#int_ls = "3_introns.CGG-TTG.list.txt"
#int_ls = "3_introns.CAG-TTG.list.txt"





"""
####### MODE 3
conv, nonc = {}, {}
for i in intron_ls:
    int_cl = i[1]
    int_id = i[0]
    [ch, rng, ori] = int_id.split(".")
    [st, ed] = rng.split("-")
    st, ed = int(st), int(ed)
    int_len = ed - st + 1

    if len_limit[0] <= int_len <= len_limit[1]:
        int_seq = ref[ch][st-1:ed]
        int_seq = RevComp(int_seq) if ori == "-" else int_seq
        
        if int_cl.startswith("Non-"):
            nonc[int_id] = int_seq
        elif int_cl.startswith("Conv"):
            conv[int_id] = int_seq
        else: pass


key1, key2 = list(conv.keys()), list(nonc.keys())
print("[PROGRESS] %s Conv and %s Nonc Introns extracted." % (len(key1), len(key2)))
seq_no = min(seq_no, len(key1), len(key2))
print("[PROGRESS] Output generated with %s intron sequences" % seq_no)
out_f1, out_f2 = "Eag.ref.NonSlide-Conv.1500-2000.all-intron.fa", "Eag.ref.NonSlide-nonConv.1500-2000.all-intron.fa"
out1, out2 = open(out_f1, "w"), open(out_f2, "w")


for i in range(len(key1)):
    conv_k = key1[i]
    conv_s = conv[conv_k]

    conv_k = ">" + conv_k

    print(conv_k, conv_s, file = out1, sep="\n")


for i in range(len(key2)):
    nonc_k = key2[i]
    nonc_s = nonc[nonc_k]

    nonc_k = ">" + nonc_k
    
    print(nonc_k, nonc_s, file = out2, sep="\n")

out1.close()
out2.close()
sys.exit()



######mode 2
for fi in inp_ls:
    out_f = ".".join(fi.split(".")[:2])+"intron-border.fasta"
    out = open(out_f,"w")
    int_ls = [i.strip().split("\t")[0] for i in open(fi)]

    x,y=0,0
    for i in int_ls:
        x += 1
        int_id = ">"+i
        i = i.split(".")
        i = i + ["0"] if len(i) == 3 else i

        scf = i[0]
        st = int(i[1].split("-")[0])
        ed = int(i[1].split("-")[1]) 
        slide = int(i[3])
        st = st + abs(slide) if st > 0 else st - abs(slide)
        ed = ed + abs(slide) if st > 0 else ed - abs(slide)
        ori = i[2]

        scf_fa = ref[scf]
        int_fa = scf_fa[st-1:ed]
        anchor5 = scf_fa[st-4:st-1]
        anchor3 = scf_fa[ed:ed+3]
        int_fa = anchor5 + int_fa + anchor3
        int_fa = RevComp(int_fa) if ori == "-" else int_fa
        int_fa = int_fa[:3].lower() + int_fa[3:-3] + int_fa[-3:].lower()
        int_fa = int_fa[3:-3]
        int_fa = int_fa[:20] + "nnnnnnnnnn" + int_fa[-20:]

        if len_limit[0] < len(int_fa):
            print(int_id, file=out)
            print(int_fa, file=out)
            y += 1
    out.close()
    print("[PROGRESS] %s/%s introns output to %s" % (y,x,out_f))
"""

