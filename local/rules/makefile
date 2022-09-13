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
	
halved_human_paralogs_biomart.tsv: human_paralogs_biomart.tsv
	awk -F"\t" '$1<$3 && $3!="" {print $0}' $< > $@
	
ohnologs_makino_parampriya_merged.filtered.tsv: ohnologs_makino_parampriya_merged.tsv halved_human_paralogs_biomart.tsv
	filter_2col 1 3 --both-orders $< < $(word 2, $^) > $@

ssd_paralogs.tsv: ohnologs_makino_parampriya_merged.tsv halved_human_paralogs_biomart.tsv
	filter_2col 1 3 -v --both-orders $< < $(word 2, $^) > $@
	
ssd_barplot.tsv: ssd_paralogs.tsv
	awk -F"\t" '{print $10}' $< | sort | uniq -c | awk -F"\t" '{print $1 "\t" $2}' > $@
	
wgd_barplot.tsv: wgd_ohnologs_makino_parampriya_merged.filtered.tsv
	awk -F"t" '{print $10}' $< | sort | uniq -c | awk -F"\t" '{print $1 "\t" $2}' > $@
	
all_barplot.tsv: halved_human_paralogs_biomart.tsv
	awk -F"\t" '{print $10}' $< | sort | uniq -c | awk -F"\t" '{print $1 "\t" $2}' > $@



