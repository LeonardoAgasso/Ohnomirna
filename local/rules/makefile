.DELETE_ON_ERROR:

makino_2010_ohnologs.xls:
	wget -O $@ 'https://www.pnas.org/action/downloadSupplement?doi=10.1073%2Fpnas.0914697107&file=st01.xls'

parampriya_2020_ohnologs_pairs_strict.tsv:
	wget -O $@ 'http://ohnologs.curie.fr/cgi-bin/DownloadBrowse.cgi?crit=[0]&org=hsapiens&opt=pairs&wgd=2R'
	#If not working try manually from Ohnologs V2 database: http://ohnologs.curie.fr/cgi-bin/BrowsePage.cgi?org=hsapiens
	
makino_2010_ohnologs.tsv: makino_2010_ohnologs.xls
	#Took sheet 7 Makino_2010_Ohnologs.xls

parampriya_2020_Ohnologs_pairs_strict_clean.tsv: parampriya_2020_ohnologs_pairs_strict.tsv
	awk -F"\t" 'NR>1 {print $$1 "\t" $$2}' $< > $@

makino_2010_ohnologs_clean.tsv: makino_2010_ohnologs.tsv
	awk -F"\t" 'NR>1 && $$1~/^ENSG/ {print $$1 "\t" $$4}' $< > $@

ohnologs_makino_parampriya_merged.tsv: makino_2010_ohnologs_clean.tsv parampriya_2020_ohnologs_pairs_strict_clean.tsv 
	cat $^ > $@

human_paralogs_biomart.tsv:
	wget -O $@ --post-file ../../local/src/mart_paralogs_query.xml 'http://www.ensembl.org/biomart/martservice'
	
human_paralogs_biomart.halved.filtered.tsv: human_paralogs_biomart.tsv
	awk -F"\t" '$$1<$$3 && $$3!="" && $$7~/^ENSP/ && $$11=="within_species_paralog" {print $$0}' $< > $@
	#Filter consists in picking paralogs which are protein-coding, and consider just A-B (since both A-B and B-A are present)
	#"within_species_paralog" because see https://www.ensembl.org/info/website/glossary.html
	
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



