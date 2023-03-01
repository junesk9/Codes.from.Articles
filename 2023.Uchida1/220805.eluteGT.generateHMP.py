#!/usr/bin/env python3
"""
Subset whole SNP data to genotypes allocated for phenotype data.
and transforms to the HapMap format


2022.08.05 KIM June-Sik
"""


import sys

fi = "call_method_75_TAIR9.csv"
fi = [i.strip().split(",") for i in open(fi)][1:]
header = fi[0]

pheno = "trait.txt"
pheno = [i.strip().split("\t")[-1] for i in open(pheno)]

#print(pheno)
#print(header)

ans = []
for h in header:
    if h in pheno:
        idx = header.index(h)
        ans.append(idx)
        #print(h)
    else: pass
#print(len(pheno)-1)
#print(ans)


out_f = "NB220.tair9.gt.hmp.txt"
out = open(out_f,"w")

header=fi[0]
header = [i for i in header if i in pheno]
new_head =["rs#","alleles","chrom","pos","strand","assembly#","center","protLSID","assayLSID","panelLSID","QCcode"]
new_head += header
print(*new_head,sep="\t",file=out)
#print(len(new_head))

new = []
for f in fi[1:]:
    tmp = []
    for idx, e in enumerate(f):
        if idx in ans:
            tmp.append(e)
        else: pass

    try:
        [ch,pos] = f[:2]
    except ValueError:
        print(f)
        print(tmp)
        continue
    gt = tmp
    alleles = set(gt)
    alleles = alleles.remove("NA") if "NA" in alleles else alleles
    alleles = sorted(list(alleles))

    # filter not biallelic SNPs
    if len(alleles) != 2:
        if len(alleles) > 3:
            print("[warning] more than two alt-SNPs; %s:%s, %s" % (ch, pos, ".".join(alleles)))
        elif len(alleles) == 1:
            print("[warning] no variation detected; %s:%s, %s" % (ch, pos, ".".join(alleles)))
        elif len(alleles) == 0:
            print("[warning] no genotype data read; %s:%s, %s" % (ch, pos, ".".join(alleles)))
        else: pass
        continue
    else: pass

    # assign an allele as the ref by its abundance
    all_count = [0,0]
    for g in gt:
        if g == alleles[0]:
            all_count[0] += 1
        elif g == alleles[1]:
            all_count[1] += 1
        else: pass
    alleles = alleles[::-1] if all_count[1] > all_count[0] else alleles
    
    snp_id = ":".join([ch,  pos] + alleles)
    alleles = "/".join(alleles)
    snp_info = [snp_id, alleles, ch, pos, "+", "NA", "NA", "NA", "NA", "NA", "NA"]
    #print(snp_info)
    #sys.exit()
    tmp = snp_info + gt
    print(*tmp,sep="\t",file=out)
    #sys.exit()
out.close()

  
