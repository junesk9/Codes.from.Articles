#!/usr/bin/env python3
"""
Comparing VCFs then identify the category information of the SNP input.
1,3 for SKY, 4-6 for SAC, J247 for Haruna Nijo

2022.02.27 Junesk9
"""

import sys, subprocess

vcf = sys.argv[-1]
new_vcf = ".".join(vcf.split(".")[:-1])+".diff.vcf"
vcf = [i.strip() for i in open(vcf)]

Refer_J247 = False
if Refer_J247:
    J247_vcf = "mpileup.filt.snp.vcf"
    J247_vcf = [i.strip() for i in open(J247_vcf)]
    J247_dict = {}
    for i in J247_vcf:
        if i.startswith("#"):
            pass
        else:
            i = i.split()
            snp_id = i[2]
            gt = i[-2].split(":")[0]
            J247_dict[snp_id] = gt
print("[PROGRESS] LOAD J247 VCF done")


out = open(new_vcf,"w")
for en,v in enumerate(vcf):
    #v = v.split()[:-1] ## remove the last (#6) sample data 
    v = v.split()
    if v[0].startswith("#"):
        print(*v, file=out, sep="\t")
    else:
        vcf_info = v[:9]
        snp_id = v[2]
        diff_side = "."
        gt = v[9:]
        gt = [i.split(":")[0] for i in gt]
        acc1 = set(gt[:2])
        acc2 = set(gt[2:4])
        if len(acc1) == len(acc2) == 1:
            if "0/1" in gt:
                pass
            elif "./." in gt:
                pass
            elif gt[0] == gt[2] == "0/0":
                pass
            elif gt[0] == gt[2] == "1/1":
                diff_side = "Common"
            elif gt[0] != gt[2]:
                if gt[0] == "0/0":
                    diff_side = "SAC"
                elif gt[0] != "0/0":
                    diff_side = "SKY"
                else:
                    print(snp_id, gt[0], gt[2])

        if Refer_J247:
            if not snp_id in J247_dict.keys():
                J247 = False
            else:
                J247 = J247_dict[snp_id]
                if J247 == "1/1":
                    J247 = True
                else:
                    J247 = False
        if not diff_side == ".":
            diff_side = diff_side + ";J247" if J247 else diff_side
        else:
            diff_side = "J247"

        v[6] = diff_side  
        if not v[0].startswith("HORVU"):
            print(*v, sep="\t", file=out)
        else: pass
out.close()
