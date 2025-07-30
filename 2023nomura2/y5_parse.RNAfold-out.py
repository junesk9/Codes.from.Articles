#!/usr/bin/env python3


import glob
import sys

subset = True

fold_fi = sys.argv[-1]
fold_fi = glob.glob("./%s/*.fold" % fold_fi)


if subset == True:
    sub_ls_f = "Both.intron.id.list"
    sub_ls = [i.strip().split("\t")[0] for i in open(sub_ls_f)]
    sub_ls = [".".join(i.split(".")[:3])+".fold" for i in sub_ls]

    sub_fi = []
    for fi in fold_fi:
        f = fi.split("/")[-1]
        if f in sub_ls:
            sub_fi.append(fi)
        else: pass
    #print(len(sub_fi))
    #print(sub_fi[:4])
fold_fi = sub_fi if subset == True else fold_fi

#print(fold_fi[:5])
fold_dic = {}
for fi in fold_fi:
    for f in open(fi):
        f = f.strip()
        if f.startswith(">"):
            name = f[1:]
            fold_dic[name] = []
        elif f.endswith(")"):
            fold = f.split(" ")
            mfe = fold[-1]
            mfe = mfe[1:] if mfe.startswith("(") else mfe
            mfe = mfe[:-1] if mfe.endswith(")") else mfe
            fold = [fold[0], mfe]
            #print(fold, len(fold[0]))
            #sys.exit()

            fold_dic[name] = fold
        else: pass

def parseInt(id_ls):
    FiveP = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
    ThreeP = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
    mfe_ls = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0] ## -42 to 0

    for i in id_ls:
        fold = fold_dic[i]
        mfe = fold[-1]
        st20 = fold[0][:20]
        ed20 = fold[0][-20:]

        for en,s in enumerate(st20):
            if s == "(":
                FiveP[en] += 1
            else: pass
        for en,e in enumerate(ed20):
            if e == ")":
                ThreeP[en] += 1

        mfe = float(mfe) * -1
        mfe_pos = int(mfe / 2)
        mfe_pos = 20 if mfe_pos > 20 else mfe_pos
        mfe_ls[mfe_pos] += 1
        mfe_ls = mfe_ls[::-1]

    FiveP = [i/len(id_ls) for i in FiveP]
    ThreeP = [i/len(id_ls) for i in ThreeP]
    mfe_ls = [i/len(id_ls) for i in mfe_ls]

    return FiveP, ThreeP, mfe_ls

FiveP, ThreeP, mfe_ls = parseInt(list(fold_dic.keys()))
print(*FiveP)
print(*ThreeP)
print(*mfe_ls)

