***
### "Genomic signatures of Japanese malting barley breeding reflected in two modern high-quality cultivars, ‘Sukai Golden’ and ‘Sachiho Golden’"
Shin Takeda et al. (under review)

1. ***[PLINK-snp_annotator.221216.py](https://github.com/junesk9/Codes.from.Articles/blob/main/2023.taketa1/PLINK-snp_annotator.221216.py)*** - Python3 script to analyze the SNP effect to Barley coding genes.
2. ***[SNP.categorizer.221216.py](https://github.com/junesk9/Codes.from.Articles/blob/main/2023.taketa1/SNP.categorizer.221216.py)*** - Python3 script to sort the SNPs by occurance in related accessions, Haruna Nijo, Sukai Golden, and Sachiho Golden.
3. ***[star-rsem-deseq.221216.R](https://github.com/junesk9/Codes.from.Articles/blob/main/2023.taketa1/star-rsem-deseq.221216.R)*** - R script to load the [STAR](https://github.com/alexdobin/STAR) outputs and evaluate the differentially expressed genes by using [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html).

***
### "A metabolome genome-wide association study identified histidine N-pi-methyltransferase in ***Arabidopsis thaliana***"
Kai Uchida et al. (accepted) Front. Plant Sci.

1. ***[220805.eluteGT.generateHMP.py](https://github.com/junesk9/Codes.from.Articles/blob/main/2023.Uchida1/220805.eluteGT.generateHMP.py)*** - Python3 script to build the cutom SNP DB as a HapMap format to feed [GAPIT3](https://github.com/jiabowang/GAPIT3) run.  
a. The original 250k array-based Arabidopsis SNP dataset is available as [call_method_75.tar.gz](https://github.com/Gregor-Mendel-Institute/atpolydb/blob/master/250k_snp_data/call_method_75.tar.gz).
2. ***[nazuna-metabolite.GAPIT3.v2.mar23.R](https://github.com/junesk9/Codes.from.Articles/blob/main/2023.Uchida1/nazuna-metabolite.GAPIT3.v2.mar23.R)*** - R script to conduct [GAPIT3](https://github.com/jiabowang/GAPIT3)-aided GWAS and followed visualization.
3. ***[v4_PLINK-snp_annotator.py](https://github.com/junesk9/Codes.from.Articles/blob/main/2023.Uchida1/v4_PLINK-snp_annotator.py)*** - Python3 script to analyze the SNP effect to coding gene alternation.  


***
### "Metabolic and transcriptomic profiling during wheat seed development under progressive drought conditions"
Ryosuke Mega et al. (under review)

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
  

***  
### "Arabidopsis TBP-ASSOCIATED FACTOR 12 ortholog NOBIRO6 controls root elongation with unfolded protein response cofactor activity."
June-Sik Kim et al. (2022) Proc. Natl. Acad. Sci. USA  PMID: [35115407](https://pubmed.ncbi.nlm.nih.gov/35115407/)  

1. ***[200929.box.barplots.scatter.R](2022PNAS/200208.AT.sleuth-DEG.R)*** - R script for drawing various plots in the manuscript  
2. ***[200208.AT.sleuth-DEG.R](2022PNAS/200208.AT.sleuth-DEG.R)***       - R script for DEG identifying using [Kallsito-Sleuth pipeline](https://www.nature.com/articles/nmeth.4324)  
a. The raw and processed RNA-seq data are available via the NCBI GEO database [GSE163911](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE163911)  
b. The raw DNA-seq data is avaialbel via the NCBI SRA database [PRJNA687636](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA687636)  

***  

