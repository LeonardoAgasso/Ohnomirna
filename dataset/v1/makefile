Makino_2010_Ohnologs.xls:
	wget -O $@ 'https://www.pnas.org/action/downloadSupplement?doi=10.1073%2Fpnas.0914697107&file=st01.xls'

ParamPriya_2020_Ohnologs_Pairs_Strict.tsv:
	wget -O $@ 'http://ohnologs.curie.fr/cgi-bin/DownloadBrowse.cgi?crit=[0]&org=hsapiens&opt=pairs&wgd=2R'
	#If not working try manually from Ohnologs V2 database: http://ohnologs.curie.fr/cgi-bin/BrowsePage.cgi?org=hsapiens
	
Makino_2010_Ohnologs.tsv: Makino_2010_Ohnologs.xls
	# Took sheet 7 Makino_2010

## add awk rules (add a $)

human_paralogs_biomart.tsv:
	wget -O $@ --post-file mart_paralogs_query.xml 'http://www.ensembl.org/biomart/martservice'



