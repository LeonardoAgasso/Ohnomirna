.DELETE_ON_ERROR:

makino_2010_ohnologs.xls:
	wget -O $@ 'https://www.pnas.org/action/downloadSupplement?doi=10.1073%2Fpnas.0914697107&file=st01.xls'

parampriya_2020_ohnologs_pairs_strict.tsv:
	wget -O $@ 'http://ohnologs.curie.fr/cgi-bin/DownloadBrowse.cgi?crit=[0]&org=hsapiens&opt=pairs&wgd=2R'
	#If not working try manually from Ohnologs V2 database: http://ohnologs.curie.fr/cgi-bin/BrowsePage.cgi?org=hsapiens
	#NR = 4,259
	
makino_2010_ohnologs.tsv: makino_2010_ohnologs.xls
	#Took sheet 7 Makino_2010_Ohnologs.xls. 
	#NR = 9,062

parampriya_2020_Ohnologs_pairs_strict_clean.tsv: parampriya_2020_ohnologs_pairs_strict.tsv
	awk -F"\t" 'NR>1 {print $$1 "\t" $$2}' $< > $@
	#NR = 4,258

makino_2010_ohnologs_clean.tsv: makino_2010_ohnologs.tsv
	awk -F"\t" 'NR>1 && $$1~/^ENSG/ {print $$1 "\t" $$4}' $< > $@
	#NR = 9,057 
	#(correctly removed first two and least three lines)

ohnologs_makino_parampriya_merged.tsv: makino_2010_ohnologs_clean.tsv parampriya_2020_ohnologs_pairs_strict_clean.tsv 
	cat $^ | sort | uniq > $@
	#NR = 11,800
	#(correctly slightly less than 4258+9057=13315)

human_paralogs_biomart.tsv:
	wget -O $@ --post-file ../../local/src/mart_paralogs_query.xml 'http://www.ensembl.org/biomart/martservice'
	#NR = 3,607,809
	
human_paralogs_biomart.halved.filtered.tsv: human_paralogs_biomart.tsv
	awk -F"\t" '$$1<$$3 && $$3!="" && $$11=="within_species_paralog" {print $$0}' $< > $@
	#Filtering and halving consists in picking paralogs and consider just A-B (since both A-B and B-A are present)
	#"within_species_paralog" because see https://www.ensembl.org/info/website/glossary.html
	#NR = 58,417
	
ohnologs_makino_parampriya_merged.filtered.tsv: ohnologs_makino_parampriya_merged.tsv human_paralogs_biomart.halved.filtered.tsv
	filter_2col 1 3 --both-orders $< < $(word 2, $^) > $@

ssd_paralogs.tsv: ohnologs_makino_parampriya_merged.tsv human_paralogs_biomart.halved.filtered.tsv
	filter_2col 1 3 -v --both-orders $< < $(word 2, $^) > $@
	
commAncestors_list.tsv: human_paralogs_biomart.halved.filtered.tsv
	awk -F"\t" '{print $$10}' $< | sort | uniq > $@
	
commAncestors_list.labeled.tsv: comm_ancestors.tsv
	#Manual labelling with an ordinal number using https://www.ncbi.nlm.nih.gov/Taxonomy
	
barplot_ssdOnly.tsv: ssd_paralogs.tsv
	awk -F"\t" '{print $$10}' $< | sort | uniq -c | awk '{print $$1 "\t" $$2}' > $@
	
barplot_wgdOnly.tsv: ohnologs_makino_parampriya_merged.filtered.tsv
	awk -F"\t" '{print $$10}' $< | sort | uniq -c | awk '{print $$1 "\t" $$2}' > $@
	
barplot_all.tsv: human_paralogs_biomart.halved.filtered.tsv
	awk -F"\t" '{print $$10}' $< | sort | uniq -c | awk '{print $$1 "\t" $$2}' > $@


#----------------------------MOCK-MAKEFILE (bedtools intersect training)-------------------------------------------

gencode.gtf.gz:
	wget -O $@ 'https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_41/gencode.v41.annotation.gtf.gz'

gencode.gtf: gencode.gtf.gz
	gzip -d $<

gencode.geneOnly.gtf: gencode.gtf
	awk -F"\t" '$3=="gene"' $< > $@
	
gencode.geneOnly.protcod.gtf: gencode.geneOnly.gtf
	awk -F"\t" '$9~/gene_type "protein_coding"/' $< > $@
	
gencode.geneOnly.miRNA.gtf: gencode.geneOnly.gtf
	awk -F"\t" '$9~/gene_type "miRNA"/' $< > $@

gencode.geneOnly.protcod.bed: gencode.geneOnly.protcod.gtf
	awk -F"\t" '{print $1 "\t" $4 "\t" $5 "\t" $9}' $< > $@

gencode.geneOnly.miRNA.bed: gencode.geneOnly.miRNA.gtf
	awk -F"\t" '{print $1 "\t" $4 "\t" $5 "\t" $9}' $< > $@
	
gencode.genOnly.miRNA_intWithProtcod.bed: gencode.geneOnly.protcod.bed gencode.geneOnly.miRNA.bed
	bedtools intersect -a $< -b $(word 2, $^) > $@

#----------------------------END MOCK-MAKEFILE---------------------------------------------------------------------



