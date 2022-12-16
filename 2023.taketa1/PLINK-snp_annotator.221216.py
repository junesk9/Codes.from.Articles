#!/usr/bin/env python3
"""
=====PLINK OUTPUT SNP ANNOTATOR===

PLINK output contains SNP location and P-values with variable adjustment, 
while the effect of SNP to the gene is hard to identify.

This program annotate the SNP location relative to the annotated gene (GFF3), e.g; CDS, intergenic, intronic
and also evaluate the effect to codon changes, silence or effective mutation

The output is a CSV format file based on the input PLINK file, add a column including annotated SNP info.
At present, the SNPs categorize into seven classes,
1. SNP or INDEL in intergenic region
2. SNP in intronic region
3. SNP in coding region, but silence mutation
4. SNP in coding region, also affect codon
5. INDEL in intronic region
6. INDEL in coding region; in-frame
7. INDEL in coding region; disrupted frame
8. SNP|INDEL in intron border sites (GT-AG)
9. SNP|INDEL in 5-/3-UTR

[USAGE] python3 v4_PLINK-snp_annotator.py -plink [plink.out file, mandatory] or -vcf [VCF file, need biallelic and NAME labeled] 
			 	  	  -ref [genome multi-fasta file, mandatory]
			    		  -gff [GFF file for the ref.genome, mandatory]
			    		  -out [name for output file, default: changing extention of the input file]
			    		  -t [Number of threads use, default: 1]

SEP.09 2015	
Junesk9
"""
"""
==== Change Log ===
AUG.05.2022     Now accept a custom SNP list of single lines of "chr:pos:ref:alt"... not yet perfect
JUN.13 2022     Correct minor output mistakes (non-sense; mis-sense)
FEB.24 2022     Add annotating UPR-positioned variations
FEB.16 2022     Add printing CHR-parsing results during parsing; checking SNP_ID is deprecated, now always being generated newly.
JAN.16 2020	Correct the abrupt brake with non-ATGC codon (FreeBayes result)
JUN.04 2018	Correct the same folder name issue during generating for a loop running. 
MAY.25 2018	Correct bugs with chromosome indicators like "Chr", all removed during the process
JUL.26 2017	Add new annotation for intron borders; just indicates +1 +2 -2 -1 sites of introns
FEB.14 2017	Add instant SNP name generator (might be danger when input VCF are not diallelic)
FEB.13 2017	Add except treatment (failed to extract codon), almost of the error comes when the last exon issued. while the reason unknown
                Change the output format to TSV
SEP.29 2016	Fix bug to automade output file naming. 
           	Add ERROR msg for un-annotated VCF usage
JUN.21 2016	Changes input system. now either -plink or -vcf is acceptable, mode automatically selected by input
DEC.17 2015	Add part to remove uninformative "Not SNP" flags after proper parsing
OCT.8  2015	Change to load chromosome name in the reference genome loading step 
OCT.7  2015	Debugs in gff reading; evade process lines starting with "#"
SEP.23 2015	Add distinguish part to single Deletion, appeared as "*" by GATK
SEP.18 2015	Add a module detect INDEL, define INDEL in CDS or intron
SEP.17 2015	Add "-vcf" mode, to define whether the insert file is PLINK-output or VCF
SEP.14 2015	Bugs fixed, related to "N" base in triplet; defined as "Undefinable SNP"
		Add an error pass mode in SNP indexing; Still have no idea why
SEP.09 2015	first build

"""

import sys, time, pickle, os, shutil, glob
from multiprocessing import Process


##################### arguments

opt = {"-plink":"",
       "-ref":"Radish_V1.5.fa",
       "-gff":"Radish_V1.51.gff",
       "-out":"",
       "-t":1,
       "-vcf_mode" : False,
       "-vcf" : "",
       "-snp_mode": False,
       "-snp" : ""
       }

args = sys.argv[1:]

if "-plink" in args:
    id_p = args.index("-plink")
    opt["-plink"] = args[id_p+1]
    opt["-vcf_mode"] = False
    inf = opt["-plink"]
if "-vcf" in args:
    id_v = args.index("-vcf")
    opt["-vcf"] = args[id_v+1]
    opt["-vcf_mode"] = True
    inf = opt["-vcf"]
if "-snp" in args:
    id_s = args.index("-snp")
    opt["-snp"] = args[id_s+1]
    opt["-snp_mode"] = True
    inf = opt["-snp"]
if "-plink" in args and "-vcf" in args:
    print("[ERROR] Either VCF of PLINK output can be analyzed!! ")
    print()
    sys.exit()
  
if "-ref" in args:
    id_r = args.index("-ref")
    opt["-ref"] = args[id_r+1]
if "-gff" in args:
    id_g = args.index("-gff")
    opt["-gff"] = args[id_g+1]
if "-t" in args:
    id_t = args.index("-t")   ## Number of threads to use
    opt["-t"] = int(args[id_t+1])
if "-out" in args:
    id_out = args.index("-out")
    opt["-out"] = args[id_out+1]
elif len(args) > 0:
    if opt["-plink"] != "":
        prefix = ".".join(opt["-plink"].split(".")[:-1])+".annt.tsv"
    elif opt["-vcf"] != "":
        prefix = ".".join(opt["-vcf"].split(".")[:-1])+".annt.tsv"
    elif opt["-snp"] != "":
        prefix = ".".join(opt["-snp"].split(".")[:-1])+".annt.tsv"
    opt["-out"] = prefix
elif len(args) == 0:
    print(__doc__)
    sys.exit()
else: pass

plink, ref, gff, out, num_th, vcf_mode, snp_mode = inf, opt["-ref"], opt["-gff"], opt["-out"], opt["-t"], opt["-vcf_mode"], opt["-snp_mode"]
print("""
    ##### SNP annotator ####
  
    SNP list as PLINK output (-plink) or SNP-named VCF (-vcf) as INPUT  
    GFF3 file (-gff) and Genome fasta file (-ref)
    are mendotory inputs for this program

    This runs with -input      %s
                   -ref        %s
                   -gff        %s
                   -out        %s
                   -t          %s
                   -vcf_mode   %s 
                   -snp_mode   %s

    """ % (plink,ref,gff,out,num_th,vcf_mode,snp_mode))



####################### function definitions

def translate(triplet):
    if len(triplet) != 3:
        print ("[ERROR]codon translation only accepts three letters !!")
        sys.exit()
    else:
        trans_dic = {
              "TTT":"F|Phe","TTC":"F|Phe","TTA":"L|Leu","TTG":"L|Leu","TCT":"S|Ser","TCC":"S|Ser",
              "TCA":"S|Ser","TCG":"S|Ser", "TAT":"Y|Tyr","TAC":"Y|Tyr","TAA":"*|Stp","TAG":"*|Stp",
              "TGT":"C|Cys","TGC":"C|Cys","TGA":"*|Stp","TGG":"W|Trp", "CTT":"L|Leu","CTC":"L|Leu",
              "CTA":"L|Leu","CTG":"L|Leu","CCT":"P|Pro","CCC":"P|Pro","CCA":"P|Pro","CCG":"P|Pro",
              "CAT":"H|His","CAC":"H|His","CAA":"Q|Gln","CAG":"Q|Gln","CGT":"R|Arg","CGC":"R|Arg",
              "CGA":"R|Arg","CGG":"R|Arg", "ATT":"I|Ile","ATC":"I|Ile","ATA":"I|Ile","ATG":"M|Met",
              "ACT":"T|Thr","ACC":"T|Thr","ACA":"T|Thr","ACG":"T|Thr", "AAT":"N|Asn","AAC":"N|Asn",
              "AAA":"K|Lys","AAG":"K|Lys","AGT":"S|Ser","AGC":"S|Ser","AGA":"R|Arg","AGG":"R|Arg",
              "GTT":"V|Val","GTC":"V|Val","GTA":"V|Val","GTG":"V|Val","GCT":"A|Ala","GCC":"A|Ala",
              "GCA":"A|Ala","GCG":"A|Ala", "GAT":"D|Asp","GAC":"D|Asp","GAA":"E|Glu",
              "GAG":"E|Glu","GGT":"G|Gly","GGC":"G|Gly","GGA":"G|Gly","GGG":"G|Gly"}
        try:
            aa = trans_dic[triplet]
        except KeyError: ## In case of the triplet include "N" or other
            aa = triplet
        return aa

def revcomp(triplet):
    tmp = [i.upper() for i in triplet][::-1]
    rc_seq = ""
    for i in tmp:
        if i == "A":
            rc_seq += "T"
        elif i == "T":
            rc_seq += "A"
        elif i == "G":
            rc_seq += "C"
        elif i == "C":
            rc_seq += "G"
        else:
            rc_seq += i

    return rc_seq


def annt_SNP(inp_pickle, gff_dic, ref_dic):

    with open(inp_pickle,"rb") as f:
        snp_inp = pickle.load(f)
        
    snp_dic = {}
    for i in snp_inp:
        #print(i)
        #sys.exit()  
        snp = i[1]
        info = i[1].split(":")
        #ch = info[0][-1] ## Remove "Chr" from the chromosome indicator
        ch = info[0] ## now the CHR styles were checked during parsing (STEP1, STEP2, STEP3) Feb.16.2022
        try: ## For Debug mode
            site = int(info[1])
        except IndexError: ## For debug mode
            print(info)
            sys.exit()
        ref, alt = info[2], info[3]
        snp_dic[snp] = "intergenic INDEL|SNP"
        for key,val in gff_dic.items():#key = transcript AGI  #val = [chr, ori, [5UTR], [CDS], [3UTR_rng]]
            if ch == val[0] and snp_dic[snp] == "intergenic INDEL|SNP": ## May help the process speed, when "snp_dic[snp]" is modified, skip afterward steps:
                if len(val[3]) == 0: ##Error control 2022.02.26
                    print("[ERROR] No CDS region identified", key, val[3])
                    sys.exit()
                elif site in val[3]:
                    if len(ref) == len(alt) == 1: ## Only select single-base polymorphysm
                        if alt == "*": ## GATK assign single deletion as "*"
                            snp_dic[snp] = "%s INDEL 1-bp (*) in CDS, disrupted frame; not SNP" % (key)
                        else:
                            snp_dic[snp] = key ## Further description will be added in the next step.
                    elif (len(ref) - len(alt)) % 3 != 0:
                        indel = len(alt) - len(ref) 
                        snp_dic[snp] = "%s INDEL %s-bp in CDS, disrupted frame; not SNP" % (key, indel)
                    else:
                        indel = len(alt) - len(ref)
                        snp_dic[snp] = "%s INDEL %s-bp in CDS, in frame; not SNP" % (key, indel)
                    #end of if
                elif min(val[3]) <= site <= max(val[3]):
                    #print("val",val)
                    #print("site",site)
                    #print("ref",ref,"alt",alt)
                    #sys.exit()
                    if (site - 2) in val[3] or (site + 2) in val[3]:
                        snp_dic[snp] == "%s intron border INDEL|SNP" % key
                        #print(snp_dic[snp])
                        #sys.exit()
                    else:
                        if len(ref) == len(alt) == 1: ## Only select SNP
                            snp_dic[snp] = "%s intronic SNP" % key
                        else:
                            indel = len(alt) - len(ref)
                            snp_dic[snp] = "%s intronic INDEL %s bp; not SNP" % (key, indel)
            
                elif site in val[2]:
                    snp_dic[snp] = "%s 5-UTR INDEL|SNP" % key
                elif site in val[4]:
                    snp_dic[snp] = "%s 3-UTR INDEL|SNP" % key
                else: pass
                #end of if
            else: pass
            #end of if
        #end of for 
        #print(snp_dic)

        ## Additional annotation parts for SNPs at coding sequence
        if not snp_dic[snp].endswith("SNP"):
            gene = snp_dic[snp]
            bps = gff_dic[gene][3]
            ch = gff_dic[gene][0]
            ori = gff_dic[gene][1]
            snp_loci = bps.index(site)
            codon_loci = (snp_loci+1) % 3  ## Extraxt the codon containing SNP
            error = "off" ## Error skip switch
            if codon_loci == 0:
                try:
                    rng2 = [bps[snp_loci-2],bps[snp_loci-1],bps[snp_loci]]
                except IndexError:
                    print("[Warning!] %s failed to extract codon; snp_loci: %s; gene: %s ln: %s " % (snp, (snp_loci+1), gene, len(bps)))
                    rng2 = [3,4,5]
                    error = "on"
            elif codon_loci == 1:
                try:
                    rng2 = [bps[snp_loci],bps[snp_loci+1],bps[snp_loci+2]]
                except IndexError:
                    print("[Warning!] %s failed to extract codon; snp_loci: %s; gene: %s ln: %s " % (snp, (snp_loci+1), gene, len(bps)))
                    rng2 = [3,4,5]
                    error = "on"
            elif codon_loci == 2:
                try:
                    rng2 = [bps[snp_loci-1],bps[snp_loci],bps[snp_loci+1]]
                except IndexError:
                    print("[Warning!] %s failed to extract codon; snp_loci: %s; gene: %s ln: %s " % (snp, (snp_loci+1), gene, len(bps)))
                    rng2 = [3,4,5]
                    error = "on"
            ###Debugging try 17.07.27
            #if error == "on":
                #print(gene, bps, snp_loci, codon_loci)     


            ref_seq = ""
            for bp in rng2:
                ref_seq += ref_dic[ch][bp-1]

            ## Validation of the codon extraction, sould be match to SNP info
            if not error == "on": ## If in error-skip mode, skip the belowing validation steps.
                if ref_seq[codon_loci-1].upper() != ref.upper():
                    if ref_seq[codon_loci-1].upper() == alt.upper(): #swap ref & alt; for some erronous input SNP list
                        ref = ref + alt
                        ref = ref[-1].upper()
                        alt = ref[0].upper()
                    elif ref_seq[codon_loci-1].upper() not in ["A","T","G","C"]: #[20.01.16] non-ATGC bother the step
                        pass
                    else:
                        print("[Error!] codon indexing result is fault!")
                        print(ref_seq, ref_seq[codon_loci-1], "!=" ,ref.upper())
                        out.close()  ## emergency out
                        sys.exit() 
                else: pass                                  ## Make altered codon
                alt_seq = [i for i in ref_seq]
                alt_seq[codon_loci-1] = alt
                alt_seq = "".join(alt_seq)
                #print(ref_seq,alt_seq)
            else: pass

            if ori == "-": ## Reverse complementation
                ref_seq = revcomp(ref_seq).upper()
                alt_seq = revcomp(alt_seq).upper()
            else: 
                ref_seq = ref_seq.upper()
                alt_seq = alt_seq.upper()
            ## Translation codon to amino acid
            ref_aa, alt_aa = translate(ref_seq), translate(alt_seq)

            ## Record reports to each SNP
            if ori == "+":
                aa_site = snp_loci+1
            elif ori == "-":
                aa_site = len(bps) - snp_loci
            #aa_site = int((snp_loci+1)/3)
            ## Error skip & record mode
            if error == "on":
                ref_seq = "This annotation has failed!, now on debuging"
                alt_seq = "snp_loci: %s; gene: %s ln: %s " % ((snp_loci+1), gene, len(bps))
                ref_aa, aa_site, alt_aa = "","",""
            else: pass 
            ## Actual output making mode
            if ref_aa == alt_aa:
                snp_dic[snp] = "Same-Sense SNP: %s %s>%s; %s %s" % (gene, ref_seq, alt_seq, ref_aa, aa_site)
            else:
                if len(ref_aa) != 3: ## In case of the triplet includes "N"
                    if alt_aa == "*|Stp":
                        snp_dic[snp] = "Non-Sense SNP: %s %s>%s; %s %s %s" % (gene, ref_seq, alt_seq, ref_aa, aa_site, alt_aa)
                    else:
                        snp_dic[snp] = "Mis-Sense SNP: %s %s>%s; %s %s %s" % (gene, ref_seq, alt_seq, ref_aa, aa_site, alt_aa)
                else:
                    print("[Warning!] %s: %s Codon cannot be converted to amino acid" % (gene, ref_aa))
                    snp_dic[snp] = "Undefinable SNP: %s %s>%s" % (gene, ref_aa, alt_aa)
        #end of if
    #end of for

    ## Remove uninformative "not SNP" flag 
    for snp in snp_dic.keys():
        if snp_dic[snp].endswith("not SNP"):  
            snp_dic[snp] = snp_dic[snp].split(";")[0]
        else: pass
    #end of for

    out_pick = inp_pickle+".out"
    with open(out_pick,"wb") as f:
        pickle.dump(snp_dic,f)

    #job_idx = inp_pickle.split(".")[1]
    #print("STEP4; Annotating SNP, job %s, done: %s" % (job_idx, time.ctime()))


def multiproc(folder, gff_dic, ref_dic):
    p_ls = glob.glob(folder+"/*")
    p_done = []

    for proc in p_ls:
        p = Process(target=annt_SNP, args=(proc, gff_dic, ref_dic))
        p.start()
        p_done.append(p)
    
    count_job = 0
    for proc in p_done:
        count_job += 1
        if count_job == 1:
            print("STEP4; Annotating SNP, started: %s" % (time.ctime()))
        else:
            print("STEP4; Annotating SNP, %s/%s jobs done: %s" % (count_job, num_th, time.ctime()))
        proc.join()


############################################## Load gff, convert cds to genomic location

#ch_prefix = False
gff_dic = {}
chrs = []
for i in open(gff):
    if not i.startswith("#"):
        i = i.strip().split("\t")
        if i[2] in ["CDS", "three_prime_UTR", "five_prime_UTR"]:
            category = i[2]
            ch = i[0] ##for wheat genome 1A, 1B, 1D, 2A ...
            #ch = i[0][-1] ##Remove "Chr" from the chromosome indicator
            #ch = i[0][1] ##plink output expresses chr as numeric
            chrs.append(ch)
            ori = i[6]
            gene_info = i[-1].split(";")

            gi_found = False
            for gi in gene_info:
                gi = gi.strip()
                if gi.startswith("Parent="):
                    name = gi.split("=")[-1]
                    name = name.split(":")[-1] if ":" in name else name ##Wheat Ensembl GFF forms "Parent=transcript:TraesCS1A02G000100.1"
                    name = name.split(",")[0] if "," in name else name ##TAIR9 GFF forms "Parent=AT1G01010.1,AT1G01010.1-Protein;"
                    gi_found = True
                else: pass
            if gi_found == False:
                print("[ERROR] Cannot identify the gene ID loading GFF")
                print(*gene_info, sep="\t")
                sys.exit()

            st, ed = int(i[3]), int(i[4])
            rng = [i for i in range(st,ed+1)] ##backup, try to fix bug; failed Nov12th 2015
            #rng = [i for i in range(st,ed)]

            if not name in gff_dic.keys():
                gff_dic[name] = [ch,ori,[],[],[]] ## for 5UTR, CDS, 3UTR, respectively; -1 is meaningless
            else: pass

            if category == "five_prime_UTR":
                gff_dic[name][2] += rng
            elif category == "CDS":
                gff_dic[name][3] += rng
            elif category == "three_prime_UTR":
                gff_dic[name][4] += rng
            else:
                print("[ERROR] unexpected categories: %s %s" % (name, i[2]))
            #print(gff_dic[name])
            #sys.exit()

            """
            try:
                gff_dic[name] += rng
            except:
                gff_dic[name] = [ch,ori]
                gff_dic[name] += rng
                #print(gff_dic[name])
                #sys.exit()
            """
        else: pass


## Sorting BP-pos order
for key, val in gff_dic.items():
    sort_cds = sorted(val[3])
    gff_dic[key][3] = sort_cds


#print(gff_dic["TraesCS1A02G000100.1"][:10])
print("STEP1; Parsing GFF data, done: %s" % time.ctime())
chrs = sorted(list(set(chrs)))
print("STEP1; GFF %s CHRs parsed like below ...." % len(chrs))
print(*chrs[:6], sep="\t")


################################################# Load Genome
ref = [i.strip() for i in open(ref)]
ref_dic = {}
for i in ref:
    if i.startswith(">"):
        gene = i.split()[0][1:] ##for wheat genome 1A, 1B, 1D, 2A ...
        #gene = i.split()[0][1:][-1] ## [-1] for remove "Chr" from the chromosome indicator
        #gene = i[2:]  ## Since PLINK recognizes chromosome with numeric only [>R1 --> 1]
        ref_dic[gene] = []
    else:
        ref_dic[gene].append(i)
for i in ref_dic.keys():
    ref_dic[i] = "".join(ref_dic[i])


#print(ref_dic["1A"][:10]) #check point whether the genome properlly loaded or not.
print("STEP2; Loading genome data, done: %s" % time.ctime())
chrs = sorted(list(ref_dic.keys()))
print("STEP2; GENOME %s CHRs parsed like below ...." % len(chrs))
print(*chrs[:6], sep="\t")

################################################# Load PLINK output, find the including gene

print("STEP3; Load PLINK|VCF input file, parse for %s CPU(s) usage" % str(num_th))

pick = []
for i in range(0,num_th):
    pick.append([])
    i += 1

## Store PLINK (vcf) file to a list; 
plink = [i.strip().split() for i in open(plink)]
chrs = []
if vcf_mode == True:
    print("STEP3; Turn on the VCF reading mode")
    plink_tmp = []
    #otp_error = False
    for i in plink:
        header = ["chr","snp"]
        if not i[0].startswith("#"):
            """
            if i[2] == ".":
                if otp_error == False:
                    print("[Warning] The VCF files seems not annotated (No SNP name given)")
                    print("[Warning] Retry the code is recommanded after annotating SNP with BCFTOOLS as [chr:bp:ref:alt]")
                    print("[Warning] The progress continues with considering the VCF file contains diallelic variants only")
                    print()
                    otp_error = True
                i[2] = i[0][-1]+":"+i[1]+":"+i[3]+":"+i[4] ## instantly generate SNP name ##[-1] for removing "Chr" from the chromosome indicator
                #sys.exit(1)
            """
            ch = i[0]
            chrs.append(ch)
            i[2] = ch+":"+i[1]+":"+i[3]+":"+i[4] 
            line = [i[0],i[2]]
            plink_tmp.append(line)
        else: pass
        plink = plink_tmp
    #end of for
elif snp_mode == True:
    print("STEP3; Turn on the custom SNP-list reading mode")
    plink_tmp = []
    for i in plink:
        header = ["chr","snp"]
        i = "Chr" + i[0] if not i[0].startswith("Chr") else i
        if not ":" in i:
            pass
        else:
            [ch, pos, ref, alt] = i.split(":")
            chrs.append(ch)
            line = [ch, i]
            plink_tmp.append(line)
    plink = plink_tmp 

elif vcf_mode == snp_mode == False:
    header = plink.pop(0) ## Remove the header from the plink input and store at "header".
#end of if
plink_len = len(plink)

make_folder = False
while not make_folder:
    try:
        folder = "tmp_"+str(int(time.time()))
        os.makedirs(folder)
        make_folder = True
    except:
        folder = "tmp_"+str(int(time.time()))

for en,i in enumerate(plink):
    pick_idx = en % num_th
    pick[pick_idx].append(i)
for en,i in enumerate(pick):
    prefix = folder+"/inp_tmp."+str(en)+".pickle"
    with open(prefix, "wb") as f:
        pickle.dump(i, f)
#plink, pick = [], []

#print(line)
print("STEP3; parsing PLINK|VCF input file, %s SNPs, done: %s" % (plink_len, time.ctime()))
chrs = sorted(list(set(chrs)))
print("STEP3; VCF %s CHRs parsed like below ..." % len(chrs))
print(*chrs[:6], sep="\t")
#sys.exit()

####################################### Annotate SNP with multiprocessing

multiproc(folder, gff_dic, ref_dic)
print("STEP4; Annotate %s SNPs done! %s" % (plink_len, time.ctime()))


###################################### Make output file

out_ls = glob.glob(folder+"/*.out")
snp_dic = {}
for fi in out_ls:
    with open(fi,"rb") as f:
        fi = pickle.load(f)
        snp_dic.update(fi)

print("STEP4; RESULT annt_size: ",len(snp_dic.keys()))
print("STEP4; RESULT inp_size: ", len(plink))
print("STEP5; Make output file, %s: %s" % (out, time.ctime()))
out = open(out,"w")
header.append("SNP_annt")
print(*header,file=out,sep="\t")
#sys.exit()
for i in plink:
    try:
        i.append(snp_dic[i[1]])
        print(*i,sep="\t",file=out)
    except KeyError:
        pass
out.close()
shutil.rmtree(folder)
plink, gff_dic, ref_dic, snp_dic = [],[],[],[] #Flushing memory
print("STEP5; All progresses completed, %s \n" % time.ctime())

"""
################################# Print outputs

print("STEP4; writing output file: %s" % out)
out = open(out,"w")
keys=sorted(snp_dic.keys())
for key in keys:
    print(key,snp_dic[key],sep=",",file=out)
out.close()

gff_dic, ref_dic, snp_dic = [],[],[] #Flushing memory
print("SNP annotation done!: %s" % time.ctime())



    idx = gff_dic[gene].index(site)-2
    gene_seq = ""
    for i in gff_dic[gene][2:]:
        gene_seq += ref_dic[ch][i-1]
"""
