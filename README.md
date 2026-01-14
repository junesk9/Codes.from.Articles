### Here is a deposit of in-house programming codes applied in published(ing) researches that I commited.           

***
### "HvFUL accelerates the seasonal transition of field-grown barley in the winter"
June-Sik Kim et al. (under review)

1. ***[seurat.barley-field.leaf.R](https://github.com/junesk9/Codes.from.Articles/blob/main/2023.barley-field/seurat.barley-field.leaf.R)*** - R script for Seurat analysis of the field barley data 1/5. 
2. ***[seurat.barley-field.leaf.v2a.R](https://github.com/junesk9/Codes.from.Articles/blob/main/2023.barley-field/seurat.barley-field.leaf.v2a.R)*** - R script for Seurat analysis of the field barley data 2/5.
3. ***[seurat.barley-field.leaf.v3.R](https://github.com/junesk9/Codes.from.Articles/blob/main/2023.barley-field/seurat.barley-field.leaf.v3.R)*** - R script for Seurat analysis of the field barley data 3/5.
4. ***[seurat.barley-field.leaf.v3a.R](https://github.com/junesk9/Codes.from.Articles/blob/main/2023.barley-field/seurat.barley-field.leaf.v3a.R)*** - R script for Seurat analysis of the field barley data 4/5.
5. ***[seurat.barley-field.leaf.v3b.R](https://github.com/junesk9/Codes.from.Articles/blob/main/2023.barley-field/seurat.barley-field.leaf.v3b.R)*** - R script for Seurat analysis of the field barley data 5/5.
7. ***[Barley-Field.Fig01.Rmd](https://github.com/junesk9/Codes.from.Articles/blob/main/2023.barley-field/Barley-Field.Fig01.Rmd)*** - R markdown document for Figure preparation (Part 1) [html](https://github.com/junesk9/Codes.from.Articles/blob/main/2023.barley-field/Barley-Field.Fig01.html)
8. ***[Barley-Field.Fig02A.Rmd](https://github.com/junesk9/Codes.from.Articles/blob/main/2023.barley-field/Barley-Field.Fig02.Rmd)*** - R markdown document for Figure preparation (Part 2) 
9. ***[Hexamer-zscale-abundance.py](https://github.com/junesk9/Codes.from.Articles/blob/main/2023.barley-field/z_Hexamer-zscale-abundance.v5e.py)*** - Python3 script for Z-test of oligomer abundance in given DNA sequences.[SampleOutput](https://github.com/junesk9/Codes.from.Articles/blob/main/2023.barley-field/cl45markers.1000-0nt.promoter.6mer-0mis.435n.1000t.z-test.txt).  
A. The associated transcriptome/epigenome data are available from NCBI GEO SuperSeries [GSE273147](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE273147) and [GSE226906](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE226906)    

***
### "Genetic dissection of nonconventional introns reveals co-dominant noncanonical splicing code in ***Euglena***"
Toshihisa Nomura & June-Sik Kim et al. (2025) Proc. Natl. Acad. Sci. USA [40986342](https://pubmed.ncbi.nlm.nih.gov/40986342/) [Article](https://www.pnas.org/doi/10.1073/pnas.2509937122) 

1. ***[1_check_splice_boundary1.pl](https://github.com/junesk9/Codes.from.Articles/blob/main/2023nomura2/1_check_splice_boundary1.pl)*** - Perl script to extract the intron information from inputs of the genome fasta and the GFF file.
2. ***[2_check_splice_boundary2.2.pl](https://github.com/junesk9/Codes.from.Articles/blob/main/2023nomura2/2_check_splice_boundary2.2.pl)*** - Perl script to parse the #1 output and check the sequence integrity of the exon-intron boarder. STDOUT for next process.
3. ***[3_check_splice_boundary3.pl](https://github.com/junesk9/Codes.from.Articles/blob/main/2023nomura2/3_check_splice_boundary3.pl)*** - Perl script to parse the #2 output. STDOUT for next process.
4. ***[4_check_boundary_shift.pl](https://github.com/junesk9/Codes.from.Articles/blob/main/2023nomura2/4_check_boundary_shift.pl)*** - Perl script to analyze the #3 output to check the slidalbility of each exon-intron boarder. STDOUT for next process.
5. ***[5_collect.nonSlidables.py](https://github.com/junesk9/Codes.from.Articles/blob/main/2023nomura2/5_collect.nonSlidables.py)*** - Python3 script to collect the non-slidable introns from the #4 output.
6. ***[y2_Eag-intron.count.class.v3.py](https://github.com/junesk9/Codes.from.Articles/blob/main/2023nomura2/y2_Eag-intron.count.class.v3.py)*** - Python3 script to classify and to count the intron boundry motifs.
7. ***[z2_extr.intron.fasta.v3.py](https://github.com/junesk9/Codes.from.Articles/blob/main/2023nomura2/z2_extr.intron.fasta.v3.py)*** - Python3 script to extract the intron sequence as a multi-fasta file.
8. ***[z6_intron.logo.table.v3.py](https://github.com/junesk9/Codes.from.Articles/blob/main/2023nomura2/z6_intron.logo.table.v3.py)*** - Python3 script to generate intron boundry logos.
9. ***[y5_parse.RNAfold-out.py](https://github.com/junesk9/Codes.from.Articles/blob/main/2023nomura2/y5_parse.RNAfold-out.py)*** - Python3 script to count the mean folding availablity from the RNAfold output.  
A. The associated supplementary data are available from a [FigShare data collection](https://doi.org/10.6084/m9.figshare.c.7576343)  
B. The associated genome/transcriptome data are available from a NCBI SRA [PRJNA1310665](https://www.ncbi.nlm.nih.gov/sra/?term=PRJNA1310665) and [PRJDB4359](https://www.ncbi.nlm.nih.gov/sra/?term=PRJDB4359)   

***
### "wheat single-cell analysis"

***
### "wheat leaf RNA-seq with ABA treatments"
Mega R et al. (in prep.)
1. ***[Wheat-RNA-seq.R](2026mega1/wheat-Mega.2501.v1.R)*** - R script to conduct the RNA-seq data analysis for the study.

***
### "Multiomics-based assessment of the impact of airflow on diverse plant callus cultures"
June-Sik Kim et al. (2025) Sci.Data PMID: [39900939](https://www.nature.com/articles/s41597-025-04518-7) [Article](https://www.nature.com/articles/s41597-025-04518-7) [Preprint](https://doi.org/10.1101/2024.07.17.604000)

1. ***[R-codes.Figure.Prep.v2.R](https://github.com/junesk9/Codes.from.Articles/blob/main/2024Callus-ms.rna/R-codes.Figure.Prep.v2.R)*** - R script for RNA-seq/metabolome data processing and Figure preparation.  
a. The raw and processed multiomics data are available from NCBI SRA:  [PRJDB16707](https://www.ncbi.nlm.nih.gov/bioproject/PRJDB16707), [PRJDB16708](https://www.ncbi.nlm.nih.gov/bioproject/PRJDB16708), [PRJDB16709](https://www.ncbi.nlm.nih.gov/bioproject/PRJDB16709) and [PRJDB16736](https://www.ncbi.nlm.nih.gov/bioproject/PRJDB16736).  

***
### "Multi-omics signatures of diverse plant callus cultures"
June-Sik Kim et al. (2024) Plant Biotechnol. [Article](https://www.jstage.jst.go.jp/article/plantbiotechnology/41/3/41_24.0719a/_article)

1. ***[R-codes.PlantBiotech.Figure.Prep.v1.R](https://github.com/junesk9/Codes.from.Articles/blob/main/2024Callus-ms.rna/R-codes.PlantBiotech.Figure.Prep.v1.R)*** - R script for figure preparstion.  

***
### "High-efficiency genome editing by Cas12a ribonucleoprotein complex in ***Euglena gracilis***"
Toshihisa Nomura et al. (2024) Microb.Biotechnol. PMID: [38332568](https://pubmed.ncbi.nlm.nih.gov/38332568/)

1. ***[z_AmpSeq-count.py](https://github.com/junesk9/Codes.from.Articles/blob/main/2023nomura1/z_AmpSeq-count.py)*** - Python3 script to post-process the [CRISSPresso2 v2.2.7](https://github.com/pinellolab/CRISPResso2) outputs to count DNA sequence matching cases against given doner DNAs.


***
### "Genomic signatures of Japanese malting barley breeding reflected in two modern high-quality cultivars, ‘Sukai Golden’ and ‘Sachiho Golden’"
Shin Taketa et al. (2023) Breed. Sci. PMID: [38737917](https://pubmed.ncbi.nlm.nih.gov/38737917/)

1. ***[PLINK-snp_annotator.221216.py](https://github.com/junesk9/Codes.from.Articles/blob/main/2023.taketa1/PLINK-snp_annotator.221216.py)*** - Python3 script to analyze the SNP effect to Barley coding genes.
2. ***[SNP.categorizer.221216.py](https://github.com/junesk9/Codes.from.Articles/blob/main/2023.taketa1/SNP.categorizer.221216.py)*** - Python3 script to sort the SNPs by their occurance in related accessions, Haruna Nijo, Sukai Golden, and Sachiho Golden.
3. ***[star-rsem-deseq.221216.R](https://github.com/junesk9/Codes.from.Articles/blob/main/2023.taketa1/star-rsem-deseq.221216.R)*** - R script to load the [STAR](https://github.com/alexdobin/STAR) outputs and to evaluate the differentially expressed genes by using [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html).    
a. The NGS data is deposit to DDBJ DRA database [PRJDB14939](https://www.ebi.ac.uk/ena/browser/view/PRJDB14939).

***
### "A metabolome genome-wide association study identified histidine N-pi-methyltransferase in ***Arabidopsis thaliana***"
Kai Uchida et al. (2023) Front. Plant Sci. PMID: [37360714](https://pubmed.ncbi.nlm.nih.gov/37360714/) 
 
1. ***[220805.eluteGT.generateHMP.py](https://github.com/junesk9/Codes.from.Articles/blob/main/2023.Uchida1/220805.eluteGT.generateHMP.py)*** - Python3 script to build the cutom SNP DB as a HapMap format to feed [GAPIT3](https://github.com/jiabowang/GAPIT3) run.  
a. The original 250k array-based Arabidopsis SNP dataset is available as [call_method_75.tar.gz](https://github.com/Gregor-Mendel-Institute/atpolydb/blob/master/250k_snp_data/call_method_75.tar.gz).
2. ***[nazuna-metabolite.GAPIT3.v2.mar23.R](https://github.com/junesk9/Codes.from.Articles/blob/main/2023.Uchida1/nazuna-metabolite.GAPIT3.v2.mar23.R)*** - R script to conduct [GAPIT3](https://github.com/jiabowang/GAPIT3)-aided GWAS and followed visualization.
3. ***[v4_PLINK-snp_annotator.py](https://github.com/junesk9/Codes.from.Articles/blob/main/2023.Uchida1/v4_PLINK-snp_annotator.py)*** - Python3 script to analyze the SNP effect to coding gene alternation.  


***
### "Metabolic and transcriptomic profiling during wheat seed development under progressive drought conditions"
Ryosuke Mega et al. (2023) Sci. Rep. PMID: [37696863](https://pubmed.ncbi.nlm.nih.gov/37696863/)

1. ***[wheat.mega.DEG-PW.R](https://github.com/junesk9/Codes.from.Articles/blob/main/2022.mega/wheat.mega.DEG-PW.R)*** - R script to identify DEGs from pairwise comparisons.
2. ***[wheat.mega.DEG-TC.R](https://github.com/junesk9/Codes.from.Articles/blob/main/2022.mega/wheat.mega.DEG-TC.R)*** - R script to indeify DEGs from the whole timecourse comparison.
3. ***[wheat.mega.GO-KEGG.v2.R](https://github.com/junesk9/Codes.from.Articles/blob/main/2022.mega/wheat.mega.GO-KEGG.v2.R)*** - R script to identify GO/KEGG terms statistically accumulated in a given gene set with the wheat [GOSlim](https://github.com/junesk9/Codes.from.Articles/blob/main/2022.mega/wheat.GOSlim.Ensembl52.txt.gz) and [KEGG](https://github.com/junesk9/Codes.from.Articles/blob/main/2022.mega/wheat.KEGG.Ensembl49.txt.gz) datasets.

***
### "Time-series transcriptome of Brachypodium distachyon during bacterial flagellin-induced pattern-triggered immunity"
Tsubasa Ogasahara et al. (2022) Front. Plant Sci. PMID: [36186055](https://pubmed.ncbi.nlm.nih.gov/36186055/)  

1. ***[search.prom_motif.220528.py](2022.ogasahara/search.prom_motif.220528.py)*** - Python3 script to search & count the CAMTA and other cis-elements in the 1-kb Brachypodium promoter regions.


***  
### "Exome-wide variation in a diverse barley panel reveals genetic associations with ten agronomic traits in Eastern landraces"  
June-Sik Kim et al. (2023) J. Genet. Genom.  PMID: [36566016](https://pubmed.ncbi.nlm.nih.gov/36566016/) 
  
1. ***[211020.pheno-plots.R](https://github.com/junesk9/In-house.codes.published/blob/main/211020.pheno-plots.R)*** – R script for drawing various phenotype plots  
2. ***[211020.qqman-gapit3.R](https://github.com/junesk9/In-house.codes.published/blob/main/211020.qqman-gapit3.R)*** - R script for Manhattan plots from the [GAPIT3](https://github.com/jiabowang/GAPIT3) resultant .GWAS.Results.csv files.  
3. ***[211221.LDHeatmap.R](https://github.com/junesk9/In-house.codes.published/blob/main/211221.LDHeatmap.R)*** - R script to visualize LD blocks with [LDHeatmap](https://sfustatgen.github.io/LDheatmap/).  
4. ***[Density_plot.v2.R](https://github.com/junesk9/In-house.codes.published/blob/main/Density_plot.v2.R)*** - R script to draw transparent density plots.  
5. ***[IBSCv2.EXOME.snp_id.convert.py](https://github.com/junesk9/In-house.codes.published/blob/main/IBSCv2.EXOME.snp_id.convert.py)*** - Python3 script to identify the allocated gene locus for each barley SNP.  
a. The raw exome-seq sequence data was deposit to the NCBI SRA database [PRJNA825586](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA825586).  
b. The raw exome BeadChip data is available from [NBRP Barley database](http://earth.nig.ac.jp/~dclust/download/iSelect_50K_SV274_A2_DesignStrand.xlsx).  
c. The extracted exomic SNP dataset was deposit to the EVA database [PRJEB53039](https://www.ebi.ac.uk/eva/?eva-study=PRJEB53039).  
d. The match list of SV-OU [SV-OU.matchlist](https://github.com/junesk9/Codes.from.Articles/blob/main/2022.barley.exome/SV-OU.matchlist.tsv)
  

***  
### "Arabidopsis TBP-ASSOCIATED FACTOR 12 ortholog NOBIRO6 controls root elongation with unfolded protein response cofactor activity."
June-Sik Kim et al. (2022) Proc. Natl. Acad. Sci. USA  PMID: [35115407](https://pubmed.ncbi.nlm.nih.gov/35115407/)  

1. ***[200929.box.barplots.scatter.R](2022PNAS/200208.AT.sleuth-DEG.R)*** - R script for drawing various plots in the manuscript  
2. ***[200208.AT.sleuth-DEG.R](2022PNAS/200208.AT.sleuth-DEG.R)***       - R script for DEG identifying using [Kallsito-Sleuth pipeline](https://www.nature.com/articles/nmeth.4324)  
a. The raw and processed RNA-seq data are available via the NCBI GEO database [GSE163911](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE163911)  
b. The raw DNA-seq data is avaialbel via the NCBI SRA database [PRJNA687636](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA687636)  

***  

