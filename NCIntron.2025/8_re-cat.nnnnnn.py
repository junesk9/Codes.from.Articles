#!/usr/bin/env python3



import sys
import time
 
fa = "Eg-asm.chen2024.fa"
nnn = "Eg-Chen.Isoseq-NGS-merged.genome.verA1.gff.intron.filtered.ambiguity_removed.2.shift.check.nnnnnn"
outf = nnn + ".recat"
nnn = [i.strip().split("\t") for i in open(nnn)]

out = open(outf, "w")
gff = "Eg-Chen.Isoseq-NGS-merged.genome.verA1.gff"
ogff = "Eg-Chen.Isoseq-NGS-merged.genome.verA1m.gff"

try:
    waiver = "waivier_gid.ls" ##The genes/prot modifed by the intron-adjust, not to alter as feedback
    waiver = set([i.strip().split(".")[0] for i in open(waiver)])
except:
    waiver = set()

idic = {}
for n in nnn:
    slide = n[5]
    if slide == "0":
        pass
    else:
        iid = n[0]
        sl = n[7]
        idic[iid] = sl
print("[PROGRESS] %s NNNNNN wait to be re-categorized" % len(idic.keys()))

fa = [i.strip() for i in open(fa)]
dic = {}
for f in fa:
    if f.startswith(">"):
        name = f[1:]
        dic[name] = []
    else:
        seq = f
        dic[name].append(seq)
dic = {k:"".join(v) for k,v in dic.items()}
print("[PROGRESS] load genome complete")


def RevComp(seq):
    dd = {"A":"T", "T":"A", "G":"C", "C":"G", "N":"N"}

    seq = seq.upper()
    r = seq[::-1]
    rc = [dd[s] for s in r]
    rc = "".join(rc)

    return(rc)

x, y = 1, 0
mdic = {}
for k,v in idic.items():
    [ch, pos, strnd] = k.split(".")
    [st, ed] = pos.split("-")

    loc_ls = [0] + v.split(" ")
    cv_ls = []
    nc_ls = []


    for l in loc_ls:
        nst = int(st) + int(l)
        ned = int(ed) + int(l)

        iseq = dic[ch][int(nst)-1:int(ned)]
        iseq = RevComp(iseq) if strnd == "-" else iseq
        #print(iseq)
        #print(RevComp(iseq))
        #sys.exit()

        conv = iseq[:2] + "-" + iseq[-2:]
        nc1 = iseq[3:6] + "-" + iseq[-8:-5]
        nc2 = iseq[3:6] + "-" + iseq[-9:-6]

        cv_ls.append(conv)
        nc_ls.append(nc1)
        nc_ls.append(nc2)
    
    fixed = False
    if "GT-AG" in cv_ls or "GC-AG" in cv_ls:
        if "GT-AG" in cv_ls:
            recat = "GT-AG"
            idx = cv_ls.index("GT-AG")
        else:
            recat = "GC-AG"
            idx = cv_ls.index("GC-AG")
        sl = loc_ls[idx]
        fixed = True
    else: pass

    nc_motif = ["CAG-CTG", "CGG-CTG", "CTG-CAG", "CTG-CGG", "CGG-CCG"]
    for m in nc_motif:
        if m in nc_ls and fixed == False:
            recat = m
            idx = nc_ls.index(m)
            idx = int(idx / 2)
            sl = loc_ls[idx]
            fixed = True

            #print(k, m, idx, sl)
            #print(nc_ls)
            #print(loc_ls)
            #print(iseq)
            #sys.exit()

        else: pass

    if fixed == False:
        nc5 = [nc.split("-")[0] for nc in nc_ls]
        nc3 = [nc.split("-")[1] for nc in nc_ls]

        if "CAG" in nc5:
            recat = "CAG-nnn"
            idx = nc5.index("CAG")
            idx = int(idx / 2)
            sl = loc_ls[idx]
            fixed = True
        elif "CTG" in nc3 and fixed == False:
            recat = "nnn-CTG"
            idx = nc3.index("CTG")
            idx = int(idx / 2)
            sl = loc_ls[idx]
            fixed = True
        elif set(nc5).intersection(["CTG", "CGG"]) != set() and fixed == False:
            recat = "CKG-nnn"
            if "CTG" in nc5:
                idx = nc5.index("CTG")
            else:
                idx = nc5.index("CGG")
            idx = int(idx / 2)
            sl = loc_ls[idx]
            fixed = True
        elif set(nc3).intersection(["CAG", "CCG", "CGG"]) != set() and fixed == False:
            recat = "nnn-CHG"
            if "CAG" in nc3:
                idx = nc3.index("CAG")
            elif "CCG" in nc3:
                idx = nc3.index("CCG")
            else:
                idx = nc3.index("CGG")
            idx = int(idx / 2)
            sl = loc_ls[idx]
            fixed = True
        elif fixed == False:
            recat = "nnn-nnn"
            sl = 0
            fixed = True
    
    print(k, sl, recat, sep="\t",  file=out)

    if recat != "nnn-nnn":
        y += 1
        mdic[k] = sl

    x += 1
print("[PROGRESS] %s / %s nnnnnn re-categorized" % (y, x))
out.close()


print()
x, y = 0, 0
gff = [i.strip().split("\t") for i in open(gff)]
gdic = {}
for g in gff:
    ch = g[0]
    if not ch in gdic.keys():
        gdic[ch] = [g]
    else:
        gdic[ch].append(g)
print("[%s] modyfing the GFF output to %s" % (time.ctime(), ogff))

for k,v in mdic.items():
    ch = k.split(".")[0]
    [st, ed] = k.split(".")[1].split("-")
    ist = int(st) 
    ied = int(ed)

    sl = v
    nmodi = 0
    gg = gdic[ch]
    for en, g in enumerate(gg):
        ost, oed = int(g[3]), int(g[4])
        gid = {i.split("=")[0]:i.split("=")[1] for i in  g[-1].split(";")}
        gid = gid["Parent"].split(".")[0] if "Parent" in gid.keys() else gid["ID"].split(".")[0]
        if ost - 100  > ied:
            break
        if gid in waiver:
            pass
        else:
            if oed == ist - 1:
                x += 1
                nmodi += 1
                new_oed = int(g[4]) + int(sl)
                gdic[ch][en][4] = new_oed if new_oed > ost else gdic[ch][en][4]
            if ost == ied + 1:
                x += 1
                nmodi += 1
                new_ost = int(g[3]) + int(sl)
                gdic[ch][en][3] = new_ost if new_ost < oed else gdic[ch][en][3] 
    y += 1

    if y % 1000 == 0:
        print("[%s] %s:%s / %s modification applied" % (time.ctime(), y, x,len(mdic.keys())))
print("[%s] GFF modified %s/%s lines" % (time.ctime(), x, y))

out = open(ogff, "w")
for k,v in gdic.items():
    for g in v:
        print(*g, sep="\t", file=out)
out.close()
print("[%s] modified GFF file generated." % time.ctime())






