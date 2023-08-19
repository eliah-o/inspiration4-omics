
library(tidyverse)
### LOAD IN REGRESSION OUTPUT AND MERGE

ranks = c('s','g','p','c','o','f')
lmers = c('','_NON-LMER')
lmers = c('')#,'_NON-LMER')
contam = c('decontam')#$,'nodecontam')

for(i in ranks){
  if(i == 's'){
    secondrankbecausesloppynaming = 'species'
  }
  if(i == 'g'){
    secondrankbecausesloppynaming = 'genus'
  }
  if(i == 'f'){
    secondrankbecausesloppynaming = 'family'
  }
  if(i == 'o'){
    secondrankbecausesloppynaming = 'order'
  }
  if(i == 'c'){
    secondrankbecausesloppynaming = 'class'
  }
  if(i == 'p'){
    secondrankbecausesloppynaming = 'phylum'
  }
  for(ii in lmers){
    for(iii in contam){
      print(paste(i,ii,iii))
      regdata = list()
      ### GTDB
      regdata[[1]] = read.table(paste0("~/Dropbox (Mason Lab)/i4/revisions/regressionoutput/taxonomic/regression_output_oral",ii,"_bacteria_metagenomics_",i,"_xtree_",iii,"_05-0025.tsv"),sep='\t') %>% mutate(seqtype = 'METAGENOMICS') %>% mutate(regressionclass = 'Oral')   %>% mutate(dset = 'GTDB')
      regdata[[2]] = read.table(paste0("~/Dropbox (Mason Lab)/i4/revisions/regressionoutput/taxonomic/regression_output_oral",ii,"_bacteria_metatranscriptomics_",i,"_xtree_",iii,"_005-0025.tsv"),sep='\t')  %>% mutate(seqtype = 'METATRANSCRIPTOMICS') %>% mutate(regressionclass = 'Oral') %>% mutate(dset = 'GTDB')
      regdata[[3]] = read.table(paste0("~/Dropbox (Mason Lab)/i4/revisions/regressionoutput/taxonomic/regression_output_nasal",ii,"_bacteria_metagenomics_",i,"_xtree_",iii,"_05-0025.tsv"),sep='\t') %>% mutate(seqtype = 'METAGENOMICS') %>% mutate(regressionclass = 'Nasal')   %>% mutate(dset = 'GTDB')
      regdata[[4]] = read.table(paste0("~/Dropbox (Mason Lab)/i4/revisions/regressionoutput/taxonomic/regression_output_nasal",ii,"_bacteria_metatranscriptomics_",i,"_xtree_",iii,"_005-0025.tsv"),sep='\t') %>% mutate(seqtype = 'METATRANSCRIPTOMICS') %>% mutate(regressionclass = 'Nasal') %>% mutate(dset = 'GTDB')
      regdata[[5]] = read.table(paste0("~/Dropbox (Mason Lab)/i4/revisions/regressionoutput/taxonomic/regression_output_skin",ii,"_bacteria_metagenomics_",i,"_xtree_",iii,"_05-0025.tsv"),sep='\t') %>% mutate(seqtype = 'METAGENOMICS') %>% mutate(regressionclass = 'Skin') %>% mutate(dset = 'GTDB')
      regdata[[6]] = read.table(paste0("~/Dropbox (Mason Lab)/i4/revisions/regressionoutput/taxonomic/regression_output_skin",ii,"_bacteria_metatranscriptomics_",i,"_xtree_",iii,"_005-0025.tsv"),sep='\t') %>% mutate(seqtype = 'METATRANSCRIPTOMICS') %>% mutate(regressionclass = 'Skin') %>% mutate(dset = 'GTDB')
      regdata[[7]] = read.table(paste0("~/Dropbox (Mason Lab)/i4/revisions/regressionoutput/taxonomic/regression_output_skin_site_by_site",ii,"_bacteria_metagenomics_",i,"_xtree_",iii,"_05-0025.tsv"),sep='\t') %>% mutate(seqtype = 'METAGENOMICS') %>% mutate(regressionclass = 'Skin_Components') %>% mutate(dset = 'GTDB')
      regdata[[8]] =  read.table(paste0("~/Dropbox (Mason Lab)/i4/revisions/regressionoutput/taxonomic/regression_output_skin_site_by_site",ii,"_bacteria_metatranscriptomics_",i,"_xtree_",iii,"_005-0025.tsv"),sep='\t')%>% mutate(seqtype = 'METATRANSCRIPTOMICS') %>% mutate(regressionclass = 'Skin_Components') %>% mutate(dset = 'GTDB')
      
      ### METAPHLAN4
      regdata[[9]] = read.table(paste0("~/Dropbox (Mason Lab)/i4/revisions/regressionoutput/taxonomic/regression_output_oral",ii,"_bacteria_metagenomics_",i,"_metaphlan4_",iii,"_default.tsv"),sep='\t') %>% mutate(seqtype = 'METAGENOMICS') %>% mutate(regressionclass = 'Oral')   %>% mutate(dset = 'METAPHLAN4')
      regdata[[10]] = read.table(paste0("~/Dropbox (Mason Lab)/i4/revisions/regressionoutput/taxonomic/regression_output_oral",ii,"_bacteria_metatranscriptomics_",i,"_metaphlan4_",iii,"_default.tsv"),sep='\t')  %>% mutate(seqtype = 'METATRANSCRIPTOMICS') %>% mutate(regressionclass = 'Oral') %>% mutate(dset = 'METAPHLAN4')
      regdata[[11]] = read.table(paste0("~/Dropbox (Mason Lab)/i4/revisions/regressionoutput/taxonomic/regression_output_nasal",ii,"_bacteria_metagenomics_",i,"_metaphlan4_",iii,"_default.tsv"),sep='\t') %>% mutate(seqtype = 'METAGENOMICS') %>% mutate(regressionclass = 'Nasal')   %>% mutate(dset = 'METAPHLAN4')
      regdata[[12]] = read.table(paste0("~/Dropbox (Mason Lab)/i4/revisions/regressionoutput/taxonomic/regression_output_nasal",ii,"_bacteria_metatranscriptomics_",i,"_metaphlan4_",iii,"_default.tsv"),sep='\t') %>% mutate(seqtype = 'METATRANSCRIPTOMICS') %>% mutate(regressionclass = 'Nasal') %>% mutate(dset = 'METAPHLAN4')
      regdata[[13]] = read.table(paste0("~/Dropbox (Mason Lab)/i4/revisions/regressionoutput/taxonomic/regression_output_skin",ii,"_bacteria_metagenomics_",i,"_metaphlan4_",iii,"_default.tsv"),sep='\t') %>% mutate(seqtype = 'METAGENOMICS') %>% mutate(regressionclass = 'Skin') %>% mutate(dset = 'METAPHLAN4')
      regdata[[14]] = read.table(paste0("~/Dropbox (Mason Lab)/i4/revisions/regressionoutput/taxonomic/regression_output_skin",ii,"_bacteria_metatranscriptomics_",i,"_metaphlan4_",iii,"_default.tsv"),sep='\t') %>% mutate(seqtype = 'METATRANSCRIPTOMICS') %>% mutate(regressionclass = 'Skin') %>% mutate(dset = 'METAPHLAN4')
      regdata[[15]] = read.table(paste0("~/Dropbox (Mason Lab)/i4/revisions/regressionoutput/taxonomic/regression_output_skin_site_by_site",ii,"_bacteria_metagenomics_",i,"_metaphlan4_",iii,"_default.tsv"),sep='\t') %>% mutate(seqtype = 'METAGENOMICS') %>% mutate(regressionclass = 'Skin_Components') %>% mutate(dset = 'METAPHLAN4')
      regdata[[16]] =  read.table(paste0("~/Dropbox (Mason Lab)/i4/revisions/regressionoutput/taxonomic/regression_output_skin_site_by_site",ii,"_bacteria_metatranscriptomics_",i,"_metaphlan4_",iii,"_default.tsv"),sep='\t')%>% mutate(seqtype = 'METATRANSCRIPTOMICS') %>% mutate(regressionclass = 'Skin_Components') %>% mutate(dset = 'METAPHLAN4')
      
      ### BACTERIAL ASSEMBLY
      #regdata[[9]] = read.table("i4_data_packet/differential_abundance_or_expression_analysis/regression_output_oral",ii,"_bacterial_dna_",secondrankbecausesloppynaming,"_xtree_.05_.0025_ASSEMBLIES.tsv"),sep='\t') %>% mutate(seqtype = 'METAGENOMICS') %>% mutate(regressionclass = 'Oral') %>% mutate(dset = 'ASSEMBLED-BACTERIA')
      #regdata[[10]] = read.table("i4_data_packet/differential_abundance_or_expression_analysis/regression_output_oral",ii,"_bacterial_rna_",secondrankbecausesloppynaming,"_xtree_.05_.0025_ASSEMBLIES.tsv"),sep='\t') %>% mutate(seqtype = 'METATRANSCRIPTOMICS') %>% mutate(regressionclass = 'Oral') %>% mutate(dset = 'ASSEMBLED-BACTERIA')
      #regdata[[11]] = read.table("i4_data_packet/differential_abundance_or_expression_analysis/regression_output_nasal",ii,"_bacterial_dna_",secondrankbecausesloppynaming,"_xtree_.05_.0025_ASSEMBLIES.tsv"),sep='\t') %>% mutate(seqtype = 'METAGENOMICS') %>% mutate(regressionclass = 'Nasal')   %>% mutate(dset = 'ASSEMBLED-BACTERIA')
      #regdata[[12]] = read.table("i4_data_packet/differential_abundance_or_expression_analysis/regression_output_nasal",ii,"_bacterial_rna_",secondrankbecausesloppynaming,"_xtree_.05_.0025_ASSEMBLIES.tsv"),sep='\t') %>% mutate(seqtype = 'METATRANSCRIPTOMICS') %>% mutate(regressionclass = 'Nasal') %>% mutate(dset = 'ASSEMBLED-BACTERIA')
      #regdata[[13]] = read.table("i4_data_packet/differential_abundance_or_expression_analysis/regression_output_skin_bacterial_dna_",secondrankbecausesloppynaming,"_xtree_.05_.0025_ASSEMBLIES.tsv"),sep='\t') %>% mutate(seqtype = 'METAGENOMICS') %>% mutate(regressionclass = 'Skin') %>% mutate(dset = 'ASSEMBLED-BACTERIA')
      #regdata[[14]] = read.table("i4_data_packet/differential_abundance_or_expression_analysis/regression_output_skin_bacterial_rna_",secondrankbecausesloppynaming,"_xtree_.05_.0025_ASSEMBLIES.tsv"),sep='\t') %>% mutate(seqtype = 'METATRANSCRIPTOMICS') %>% mutate(regressionclass = 'Skin')  %>% mutate(dset = 'ASSEMBLED-BACTERIA')
      #regdata[[15]] = read.table("i4_data_packet/differential_abundance_or_expression_analysis/regression_output_skin_site_by_site",ii,"_bacterial_dna_",secondrankbecausesloppynaming,"_xtree_.05_.0025_ASSEMBLIES.tsv"),sep='\t') %>% mutate(seqtype = 'METAGENOMICS') %>% mutate(regressionclass = 'Skin_Components')  %>% mutate(dset = 'ASSEMBLED-BACTERIA')
      #regdata[[16]] = read.table("i4_data_packet/differential_abundance_or_expression_analysis/regression_output_skin_site_by_site",ii,"_bacterial_rna_",secondrankbecausesloppynaming,"_xtree_.05_.0025_ASSEMBLIES.tsv"),sep='\t') %>% mutate(seqtype = 'METATRANSCRIPTOMICS') %>% mutate(regressionclass = 'Skin_Components') %>% mutate(dset = 'ASSEMBLED-BACTERIA')
      
      ### VIRAL GENBANK
      
      regdata[[17]] = read.csv(paste0("~/Dropbox (Mason Lab)/i4/revisions/regressionoutput/taxonomic/regression_output_oral",ii,"_virus_metagenomics_",i,"_xtree_",iii,"_01-005.tsv"),sep='\t') %>% mutate(seqtype = 'METAGENOMICS') %>% mutate(regressionclass = 'Oral') %>% mutate(dset = 'GENBANK-VIRUSES')
      regdata[[18]] = read.csv(paste0("~/Dropbox (Mason Lab)/i4/revisions/regressionoutput/taxonomic/regression_output_oral",ii,"_virus_metatranscriptomics_",i,"_xtree_",iii,"_01-005.tsv"),sep='\t') %>% mutate(seqtype = 'METATRANSCRIPTOMICS') %>% mutate(regressionclass = 'Oral') %>% mutate(dset = 'GENBANK-VIRUSES')
      regdata[[19]] = read.csv(paste0("~/Dropbox (Mason Lab)/i4/revisions/regressionoutput/taxonomic/regression_output_nasal",ii,"_virus_metagenomics_",i,"_xtree_",iii,"_01-005.tsv"),sep='\t') %>% mutate(seqtype = 'METAGENOMICS') %>% mutate(regressionclass = 'Nasal')   %>% mutate(dset = 'GENBANK-VIRUSES')
      regdata[[20]] = read.csv(paste0("~/Dropbox (Mason Lab)/i4/revisions/regressionoutput/taxonomic/regression_output_nasal",ii,"_virus_metatranscriptomics_",i,"_xtree_",iii,"_01-005.tsv"),sep='\t') %>% mutate(seqtype = 'METATRANSCRIPTOMICS') %>% mutate(regressionclass = 'Nasal') %>% mutate(dset = 'GENBANK-VIRUSES')
      regdata[[21]] = read.csv(paste0("~/Dropbox (Mason Lab)/i4/revisions/regressionoutput/taxonomic/regression_output_skin",ii,"_virus_metagenomics_",i,"_xtree_",iii,"_01-005.tsv"),sep='\t') %>% mutate(seqtype = 'METAGENOMICS') %>% mutate(regressionclass = 'Skin') %>% mutate(dset = 'GENBANK-VIRUSES')
      regdata[[22]] = read.csv(paste0("~/Dropbox (Mason Lab)/i4/revisions/regressionoutput/taxonomic/regression_output_skin",ii,"_virus_metatranscriptomics_",i,"_xtree_",iii,"_01-005.tsv"),sep='\t') %>% mutate(seqtype = 'METATRANSCRIPTOMICS') %>% mutate(regressionclass = 'Skin') %>% mutate(dset = 'GENBANK-VIRUSES')
      regdata[[23]] = read.csv(paste0("~/Dropbox (Mason Lab)/i4/revisions/regressionoutput/taxonomic/regression_output_skin_site_by_site",ii,"_virus_metagenomics_",i,"_xtree_",iii,"_01-005.tsv"),sep='\t') %>% mutate(seqtype = 'METAGENOMICS') %>% mutate(regressionclass = 'Skin_Components') %>% mutate(dset = 'GENBANK-VIRUSES')
      regdata[[24]] = read.csv(paste0("~/Dropbox (Mason Lab)/i4/revisions/regressionoutput/taxonomic/regression_output_skin_site_by_site",ii,"_virus_metatranscriptomics_",i,"_xtree_",iii,"_01-005.tsv"),sep='\t') %>% mutate(seqtype = 'METATRANSCRIPTOMICS') %>% mutate(regressionclass = 'Skin_Components') %>% mutate(dset = 'GENBANK-VIRUSES')
      
      ### VIRAL ASSEMBLED
      
      #regdata[[25]] = read.csv("i4_data_packet/differential_abundance_or_expression_analysis/regression_output_oral",ii,"_viral_dna_",secondrankbecausesloppynaming,"_xtree_.01_.0.005_ASSEMBLIES.tsv"),sep='\t') %>% mutate(seqtype = 'METAGENOMICS') %>% mutate(regressionclass = 'Oral') %>% mutate(dset = 'ASSEMBLED-VIRUSES')
      #regdata[[26]] = read.csv("i4_data_packet/differential_abundance_or_expression_analysis/regression_output_oral",ii,"_viral_rna_",secondrankbecausesloppynaming,"_xtree_.01_.0.005_ASSEMBLIES.tsv"),sep='\t') %>% mutate(seqtype = 'METATRANSCRIPTOMICS') %>% mutate(regressionclass = 'Oral') %>% mutate(dset = 'ASSEMBLED-VIRUSES')
      #regdata[[27]] = read.csv("i4_data_packet/differential_abundance_or_expression_analysis/regression_output_nasal",ii,"_viral_dna_",secondrankbecausesloppynaming,"_xtree_.01_.0.005_ASSEMBLIES.tsv"),sep='\t') %>% mutate(seqtype = 'METAGENOMICS') %>% mutate(regressionclass = 'Nasal')   %>% mutate(dset = 'ASSEMBLED-VIRUSES')
      #regdata[[28]] = read.csv("i4_data_packet/differential_abundance_or_expression_analysis/regression_output_nasal",ii,"_viral_rna_",secondrankbecausesloppynaming,"_xtree_.01_.0.005_ASSEMBLIES.tsv"),sep='\t') %>% mutate(seqtype = 'METATRANSCRIPTOMICS') %>% mutate(regressionclass = 'Nasal') %>% mutate(dset = 'ASSEMBLED-VIRUSES')
      #regdata[[29]] = read.csv("i4_data_packet/differential_abundance_or_expression_analysis/regression_output_skin_viral_dna_",secondrankbecausesloppynaming,"_xtree_.01_.0.005_ASSEMBLIES.tsv"),sep='\t') %>% mutate(seqtype = 'METAGENOMICS') %>% mutate(regressionclass = 'Skin') %>% mutate(dset = 'ASSEMBLED-VIRUSES')
      #regdata[[30]] = read.csv("i4_data_packet/differential_abundance_or_expression_analysis/regression_output_skin_viral_rna_",secondrankbecausesloppynaming,"_xtree_.01_.0.005_ASSEMBLIES.tsv"),sep='\t') %>% mutate(seqtype = 'METATRANSCRIPTOMICS') %>% mutate(regressionclass = 'Skin') %>% mutate(dset = 'ASSEMBLED-VIRUSES')
      #regdata[[31]] = read.csv("i4_data_packet/differential_abundance_or_expression_analysis/regression_output_skin_site_by_site",ii,"_viral_dna_",secondrankbecausesloppynaming,"_xtree_.01_.0.005_ASSEMBLIES.tsv"),sep='\t') %>% mutate(seqtype = 'METAGENOMICS') %>% mutate(regressionclass = 'Skin_Components') %>% mutate(dset = 'ASSEMBLED-VIRUSES')
      #regdata[[32]] = read.csv("i4_data_packet/differential_abundance_or_expression_analysis/regression_output_skin_site_by_site",ii,"_viral_rna_",secondrankbecausesloppynaming,"_xtree_.01_.0.005_ASSEMBLIES.tsv"),sep='\t') %>% mutate(seqtype = 'METATRANSCRIPTOMICS') %>% mutate(regressionclass = 'Skin_Components') %>% mutate(dset = 'ASSEMBLED-VIRUSES')
      
      ### KRAKEN2 -- CONF 0.2 / MASK
      regdata[[25]] = read.csv(paste0("~/Dropbox (Mason Lab)/i4/revisions/regressionoutput/taxonomic/regression_output_oral",ii,"_bac-vir-fung_metagenomics_",secondrankbecausesloppynaming,"_kraken2_",iii,"_conf-mask.tsv"),sep='\t') %>% mutate(seqtype = 'METAGENOMICS') %>% mutate(regressionclass = 'Oral') %>% mutate(dset = 'KRAKEN2-CONF0.2-MASK')
      regdata[[26]] = read.csv(paste0("~/Dropbox (Mason Lab)/i4/revisions/regressionoutput/taxonomic/regression_output_oral",ii,"_bac-vir-fung_metatranscriptomics_",secondrankbecausesloppynaming,"_kraken2_",iii,"_conf-mask.tsv"),sep='\t') %>% mutate(seqtype = 'METATRANSCRIPTOMICS') %>% mutate(regressionclass = 'Oral') %>% mutate(dset = 'KRAKEN2-CONF0.2-MASK')
      regdata[[27]] = read.csv(paste0("~/Dropbox (Mason Lab)/i4/revisions/regressionoutput/taxonomic/regression_output_nasal",ii,"_bac-vir-fung_metagenomics_",secondrankbecausesloppynaming,"_kraken2_",iii,"_conf-mask.tsv"),sep='\t') %>% mutate(seqtype = 'METAGENOMICS') %>% mutate(regressionclass = 'Nasal')   %>% mutate(dset = 'KRAKEN2-CONF0.2-MASK')
      regdata[[28]] = read.csv(paste0("~/Dropbox (Mason Lab)/i4/revisions/regressionoutput/taxonomic/regression_output_nasal",ii,"_bac-vir-fung_metatranscriptomics_",secondrankbecausesloppynaming,"_kraken2_",iii,"_conf-mask.tsv"),sep='\t') %>% mutate(seqtype = 'METATRANSCRIPTOMICS') %>% mutate(regressionclass = 'Nasal') %>% mutate(dset = 'KRAKEN2-CONF0.2-MASK')
      regdata[[29]] = read.csv(paste0("~/Dropbox (Mason Lab)/i4/revisions/regressionoutput/taxonomic/regression_output_skin",ii,"_bac-vir-fung_metagenomics_",secondrankbecausesloppynaming,"_kraken2_",iii,"_conf-mask.tsv"),sep='\t') %>% mutate(seqtype = 'METAGENOMICS') %>% mutate(regressionclass = 'Skin') %>% mutate(dset = 'KRAKEN2-CONF0.2-MASK')
      regdata[[30]] = read.csv(paste0("~/Dropbox (Mason Lab)/i4/revisions/regressionoutput/taxonomic/regression_output_skin",ii,"_bac-vir-fung_metatranscriptomics_",secondrankbecausesloppynaming,"_kraken2_",iii,"_conf-mask.tsv"),sep='\t') %>% mutate(seqtype = 'METATRANSCRIPTOMICS') %>% mutate(regressionclass = 'Skin') %>% mutate(dset = 'KRAKEN2-CONF0.2-MASK')
      regdata[[31]] = read.csv(paste0("~/Dropbox (Mason Lab)/i4/revisions/regressionoutput/taxonomic/regression_output_skin_site_by_site",ii,"_bac-vir-fung_metagenomics_",secondrankbecausesloppynaming,"_kraken2_",iii,"_conf-mask.tsv"),sep='\t') %>% mutate(seqtype = 'METAGENOMICS') %>% mutate(regressionclass = 'Skin_Components') %>% mutate(dset = 'KRAKEN2-CONF0.2-MASK')
      regdata[[32]] = read.csv(paste0("~/Dropbox (Mason Lab)/i4/revisions/regressionoutput/taxonomic/regression_output_skin_site_by_site",ii,"_bac-vir-fung_metatranscriptomics_",secondrankbecausesloppynaming,"_kraken2_",iii,"_conf-mask.tsv"),sep='\t') %>% mutate(seqtype = 'METATRANSCRIPTOMICS') %>% mutate(regressionclass = 'Skin_Components') %>% mutate(dset = 'KRAKEN2-CONF0.2-MASK')
      
      ### KRAKEN2 -- CONF 0.2 / NO MASK
      regdata[[33]] = read.csv(paste0("~/Dropbox (Mason Lab)/i4/revisions/regressionoutput/taxonomic/regression_output_oral",ii,"_bac-vir-fung_metagenomics_",secondrankbecausesloppynaming,"_kraken2_",iii,"_nomask-conf.tsv"),sep='\t') %>% mutate(seqtype = 'METAGENOMICS') %>% mutate(regressionclass = 'Oral') %>% mutate(dset = 'KRAKEN2-CONF0.2-NOMASK')
      regdata[[34]] = read.csv(paste0("~/Dropbox (Mason Lab)/i4/revisions/regressionoutput/taxonomic/regression_output_oral",ii,"_bac-vir-fung_metatranscriptomics_",secondrankbecausesloppynaming,"_kraken2_",iii,"_nomask-conf.tsv"),sep='\t') %>% mutate(seqtype = 'METATRANSCRIPTOMICS') %>% mutate(regressionclass = 'Oral') %>% mutate(dset = 'KRAKEN2-CONF0.2-NOMASK')
      regdata[[35]] = read.csv(paste0("~/Dropbox (Mason Lab)/i4/revisions/regressionoutput/taxonomic/regression_output_nasal",ii,"_bac-vir-fung_metagenomics_",secondrankbecausesloppynaming,"_kraken2_",iii,"_nomask-conf.tsv"),sep='\t') %>% mutate(seqtype = 'METAGENOMICS') %>% mutate(regressionclass = 'Nasal')   %>% mutate(dset = 'KRAKEN2-CONF0.2-NOMASK')
      regdata[[36]] = read.csv(paste0("~/Dropbox (Mason Lab)/i4/revisions/regressionoutput/taxonomic/regression_output_nasal",ii,"_bac-vir-fung_metatranscriptomics_",secondrankbecausesloppynaming,"_kraken2_",iii,"_nomask-conf.tsv"),sep='\t') %>% mutate(seqtype = 'METATRANSCRIPTOMICS') %>% mutate(regressionclass = 'Nasal') %>% mutate(dset = 'KRAKEN2-CONF0.2-NOMASK')
      regdata[[37]] = read.csv(paste0("~/Dropbox (Mason Lab)/i4/revisions/regressionoutput/taxonomic/regression_output_skin",ii,"_bac-vir-fung_metagenomics_",secondrankbecausesloppynaming,"_kraken2_",iii,"_nomask-conf.tsv"),sep='\t') %>% mutate(seqtype = 'METAGENOMICS') %>% mutate(regressionclass = 'Skin') %>% mutate(dset = 'KRAKEN2-CONF0.2-NOMASK')
      regdata[[38]] = read.csv(paste0("~/Dropbox (Mason Lab)/i4/revisions/regressionoutput/taxonomic/regression_output_skin",ii,"_bac-vir-fung_metatranscriptomics_",secondrankbecausesloppynaming,"_kraken2_",iii,"_nomask-conf.tsv"),sep='\t') %>% mutate(seqtype = 'METATRANSCRIPTOMICS') %>% mutate(regressionclass = 'Skin') %>% mutate(dset = 'KRAKEN2-CONF0.2-NOMASK')
      regdata[[39]] = read.csv(paste0("~/Dropbox (Mason Lab)/i4/revisions/regressionoutput/taxonomic/regression_output_skin_site_by_site",ii,"_bac-vir-fung_metagenomics_",secondrankbecausesloppynaming,"_kraken2_",iii,"_nomask-conf.tsv"),sep='\t') %>% mutate(seqtype = 'METAGENOMICS') %>% mutate(regressionclass = 'Skin_Components') %>% mutate(dset = 'KRAKEN2-CONF0.2-NOMASK')
      regdata[[40]] = read.csv(paste0("~/Dropbox (Mason Lab)/i4/revisions/regressionoutput/taxonomic/regression_output_skin_site_by_site",ii,"_bac-vir-fung_metatranscriptomics_",secondrankbecausesloppynaming,"_kraken2_",iii,"_nomask-conf.tsv"),sep='\t') %>% mutate(seqtype = 'METATRANSCRIPTOMICS') %>% mutate(regressionclass = 'Skin_Components') %>% mutate(dset = 'KRAKEN2-CONF0.2-NOMASK')
      
      ### KRAKEN2 -- CONF 0.0 / MASK
      regdata[[41]] = read.csv(paste0("~/Dropbox (Mason Lab)/i4/revisions/regressionoutput/taxonomic/regression_output_oral",ii,"_bac-vir-fung_metagenomics_",secondrankbecausesloppynaming,"_kraken2_",iii,"_mask-noconf.tsv"),sep='\t') %>% mutate(seqtype = 'METAGENOMICS') %>% mutate(regressionclass = 'Oral') %>% mutate(dset = 'KRAKEN2-CONF0.0-MASK')
      regdata[[42]] = read.csv(paste0("~/Dropbox (Mason Lab)/i4/revisions/regressionoutput/taxonomic/regression_output_oral",ii,"_bac-vir-fung_metatranscriptomics_",secondrankbecausesloppynaming,"_kraken2_",iii,"_mask-noconf.tsv"),sep='\t') %>% mutate(seqtype = 'METATRANSCRIPTOMICS') %>% mutate(regressionclass = 'Oral') %>% mutate(dset = 'KRAKEN2-CONF0.0-MASK')
      regdata[[43]] = read.csv(paste0("~/Dropbox (Mason Lab)/i4/revisions/regressionoutput/taxonomic/regression_output_nasal",ii,"_bac-vir-fung_metagenomics_",secondrankbecausesloppynaming,"_kraken2_",iii,"_mask-noconf.tsv"),sep='\t') %>% mutate(seqtype = 'METAGENOMICS') %>% mutate(regressionclass = 'Nasal')   %>% mutate(dset = 'KRAKEN2-CONF0.0-MASK')
      regdata[[44]] = read.csv(paste0("~/Dropbox (Mason Lab)/i4/revisions/regressionoutput/taxonomic/regression_output_nasal",ii,"_bac-vir-fung_metatranscriptomics_",secondrankbecausesloppynaming,"_kraken2_",iii,"_mask-noconf.tsv"),sep='\t') %>% mutate(seqtype = 'METATRANSCRIPTOMICS') %>% mutate(regressionclass = 'Nasal') %>% mutate(dset = 'KRAKEN2-CONF0.0-MASK')
      regdata[[45]] = read.csv(paste0("~/Dropbox (Mason Lab)/i4/revisions/regressionoutput/taxonomic/regression_output_skin",ii,"_bac-vir-fung_metagenomics_",secondrankbecausesloppynaming,"_kraken2_",iii,"_mask-noconf.tsv"),sep='\t') %>% mutate(seqtype = 'METAGENOMICS') %>% mutate(regressionclass = 'Skin') %>% mutate(dset = 'KRAKEN2-CONF0.0-MASK')
      regdata[[46]] = read.csv(paste0("~/Dropbox (Mason Lab)/i4/revisions/regressionoutput/taxonomic/regression_output_skin",ii,"_bac-vir-fung_metatranscriptomics_",secondrankbecausesloppynaming,"_kraken2_",iii,"_mask-noconf.tsv"),sep='\t') %>% mutate(seqtype = 'METATRANSCRIPTOMICS') %>% mutate(regressionclass = 'Skin') %>% mutate(dset = 'KRAKEN2-CONF0.0-MASK')
      regdata[[47]] = read.csv(paste0("~/Dropbox (Mason Lab)/i4/revisions/regressionoutput/taxonomic/regression_output_skin_site_by_site",ii,"_bac-vir-fung_metagenomics_",secondrankbecausesloppynaming,"_kraken2_",iii,"_mask-noconf.tsv"),sep='\t') %>% mutate(seqtype = 'METAGENOMICS') %>% mutate(regressionclass = 'Skin_Components') %>% mutate(dset = 'KRAKEN2-CONF0.0-MASK')
      regdata[[48]] = read.csv(paste0("~/Dropbox (Mason Lab)/i4/revisions/regressionoutput/taxonomic/regression_output_skin_site_by_site",ii,"_bac-vir-fung_metatranscriptomics_",secondrankbecausesloppynaming,"_kraken2_",iii,"_mask-noconf.tsv"),sep='\t') %>% mutate(seqtype = 'METATRANSCRIPTOMICS') %>% mutate(regressionclass = 'Skin_Components') %>% mutate(dset = 'KRAKEN2-CONF0.0-MASK')
      
      ### KRAKEN2 -- CONF 0.0 / NO MASK
      regdata[[49]] = read.csv(paste0("~/Dropbox (Mason Lab)/i4/revisions/regressionoutput/taxonomic/regression_output_oral",ii,"_bac-vir-fung_metagenomics_",secondrankbecausesloppynaming,"_kraken2_",iii,"_nomask-noconf.tsv"),sep='\t') %>% mutate(seqtype = 'METAGENOMICS') %>% mutate(regressionclass = 'Oral') %>% mutate(dset = 'KRAKEN2-CONF0.0-NOMASK')
      regdata[[50]] = read.csv(paste0("~/Dropbox (Mason Lab)/i4/revisions/regressionoutput/taxonomic/regression_output_oral",ii,"_bac-vir-fung_metatranscriptomics_",secondrankbecausesloppynaming,"_kraken2_",iii,"_nomask-noconf.tsv"),sep='\t') %>% mutate(seqtype = 'METATRANSCRIPTOMICS') %>% mutate(regressionclass = 'Oral') %>% mutate(dset = 'KRAKEN2-CONF0.0-NOMASK')
      regdata[[51]] = read.csv(paste0("~/Dropbox (Mason Lab)/i4/revisions/regressionoutput/taxonomic/regression_output_nasal",ii,"_bac-vir-fung_metagenomics_",secondrankbecausesloppynaming,"_kraken2_",iii,"_nomask-noconf.tsv"),sep='\t') %>% mutate(seqtype = 'METAGENOMICS') %>% mutate(regressionclass = 'Nasal')   %>% mutate(dset = 'KRAKEN2-CONF0.0-NOMASK')
      regdata[[52]] = read.csv(paste0("~/Dropbox (Mason Lab)/i4/revisions/regressionoutput/taxonomic/regression_output_nasal",ii,"_bac-vir-fung_metatranscriptomics_",secondrankbecausesloppynaming,"_kraken2_",iii,"_nomask-noconf.tsv"),sep='\t') %>% mutate(seqtype = 'METATRANSCRIPTOMICS') %>% mutate(regressionclass = 'Nasal') %>% mutate(dset = 'KRAKEN2-CONF0.0-NOMASK')
      regdata[[53]] = read.csv(paste0("~/Dropbox (Mason Lab)/i4/revisions/regressionoutput/taxonomic/regression_output_skin",ii,"_bac-vir-fung_metagenomics_",secondrankbecausesloppynaming,"_kraken2_",iii,"_nomask-noconf.tsv"),sep='\t') %>% mutate(seqtype = 'METAGENOMICS') %>% mutate(regressionclass = 'Skin') %>% mutate(dset = 'KRAKEN2-CONF0.0-NOMASK')
      regdata[[54]] = read.csv(paste0("~/Dropbox (Mason Lab)/i4/revisions/regressionoutput/taxonomic/regression_output_skin",ii,"_bac-vir-fung_metatranscriptomics_",secondrankbecausesloppynaming,"_kraken2_",iii,"_nomask-noconf.tsv"),sep='\t') %>% mutate(seqtype = 'METATRANSCRIPTOMICS') %>% mutate(regressionclass = 'Skin') %>% mutate(dset = 'KRAKEN2-CONF0.0-NOMASK')
      regdata[[55]] = read.csv(paste0("~/Dropbox (Mason Lab)/i4/revisions/regressionoutput/taxonomic/regression_output_skin_site_by_site",ii,"_bac-vir-fung_metagenomics_",secondrankbecausesloppynaming,"_kraken2_",iii,"_nomask-noconf.tsv"),sep='\t') %>% mutate(seqtype = 'METAGENOMICS') %>% mutate(regressionclass = 'Skin_Components') %>% mutate(dset = 'KRAKEN2-CONF0.0-NOMASK')
      regdata[[56]] = read.csv(paste0("~/Dropbox (Mason Lab)/i4/revisions/regressionoutput/taxonomic/regression_output_skin_site_by_site",ii,"_bac-vir-fung_metatranscriptomics_",secondrankbecausesloppynaming,"_kraken2_",iii,"_nomask-noconf.tsv"),sep='\t') %>% mutate(seqtype = 'METATRANSCRIPTOMICS') %>% mutate(regressionclass = 'Skin_Components') %>% mutate(dset = 'KRAKEN2-CONF0.0-NOMASK')
      


      ### PHANTA
      regdata[[57]] = read.delim(paste0("~/Dropbox (Mason Lab)/i4/revisions/regressionoutput/taxonomic/regression_output_oral",ii,"_virus_metagenomics_",i,"_phanta_",iii,"_default.tsv"),sep='\t') %>% mutate(seqtype = 'METAGENOMICS') %>% mutate(regressionclass = 'Oral') %>% mutate(dset = 'PHANTA-VIRUSES')
      regdata[[58]] = read.delim(paste0("~/Dropbox (Mason Lab)/i4/revisions/regressionoutput/taxonomic/regression_output_oral",ii,"_virus_metatranscriptomics_",i,"_phanta_",iii,"_default.tsv"),sep='\t') %>% mutate(seqtype = 'METATRANSCRIPTOMICS') %>% mutate(regressionclass = 'Oral') %>% mutate(dset = 'PHANTA-VIRUSES')
      regdata[[59]] = read.delim(paste0("~/Dropbox (Mason Lab)/i4/revisions/regressionoutput/taxonomic/regression_output_nasal",ii,"_virus_metagenomics_",i,"_phanta_",iii,"_default.tsv"),sep='\t') %>% mutate(seqtype = 'METAGENOMICS') %>% mutate(regressionclass = 'Nasal')   %>% mutate(dset = 'PHANTA-VIRUSES')
      regdata[[60]] = read.delim(paste0("~/Dropbox (Mason Lab)/i4/revisions/regressionoutput/taxonomic/regression_output_nasal",ii,"_virus_metatranscriptomics_",i,"_phanta_",iii,"_default.tsv"),sep='\t') %>% mutate(seqtype = 'METATRANSCRIPTOMICS') %>% mutate(regressionclass = 'Nasal') %>% mutate(dset = 'PHANTA-VIRUSES')
      regdata[[61]] = read.delim(paste0("~/Dropbox (Mason Lab)/i4/revisions/regressionoutput/taxonomic/regression_output_skin",ii,"_virus_metagenomics_",i,"_phanta_",iii,"_default.tsv"),sep='\t') %>% mutate(seqtype = 'METAGENOMICS') %>% mutate(regressionclass = 'Skin') %>% mutate(dset = 'PHANTA-VIRUSES')
      regdata[[62]] = read.delim(paste0("~/Dropbox (Mason Lab)/i4/revisions/regressionoutput/taxonomic/regression_output_skin",ii,"_virus_metatranscriptomics_",i,"_phanta_",iii,"_default.tsv"),sep='\t') %>% mutate(seqtype = 'METATRANSCRIPTOMICS') %>% mutate(regressionclass = 'Skin') %>% mutate(dset = 'PHANTA-VIRUSES')
      regdata[[63]] = read.delim(paste0("~/Dropbox (Mason Lab)/i4/revisions/regressionoutput/taxonomic/regression_output_skin_site_by_site",ii,"_virus_metagenomics_",i,"_phanta_",iii,"_default.tsv"),sep='\t') %>% mutate(seqtype = 'METAGENOMICS') %>% mutate(regressionclass = 'Skin_Components') %>% mutate(dset = 'PHANTA-VIRUSES')
      regdata[[64]] = read.delim(paste0("~/Dropbox (Mason Lab)/i4/revisions/regressionoutput/taxonomic/regression_output_skin_site_by_site",ii,"_virus_metatranscriptomics_",i,"_phanta_",iii,"_default.tsv"),sep='\t') %>% mutate(seqtype = 'METATRANSCRIPTOMICS') %>% mutate(regressionclass = 'Skin_Components') %>% mutate(dset = 'PHANTA-VIRUSES')
      
      if(i=='s'){
      ### GENE
      regdata[[65]] = read.csv(paste0("/Users/bradentemp/Dropbox (Mason Lab)/i4/revisions/regressionoutput/gene/regression_output_oral",ii,"_merged_gene_metag_FILTERED.tsv"),sep='\t') %>% mutate(seqtype = 'METAGENOMICS') %>% mutate(regressionclass = 'Oral') %>% mutate(dset = 'GENE-CATALOG')
      regdata[[66]] = read.csv(paste0("/Users/bradentemp/Dropbox (Mason Lab)/i4/revisions/regressionoutput/gene/regression_output_oral",ii,"_merged_gene_metat_FILTERED.tsv"),sep='\t') %>% mutate(seqtype = 'METATRANSCRIPTOMICS') %>% mutate(regressionclass = 'Oral') %>% mutate(dset = 'GENE-CATALOG')
      regdata[[67]] = read.csv(paste0("/Users/bradentemp/Dropbox (Mason Lab)/i4/revisions/regressionoutput/gene/regression_output_nasal",ii,"_merged_gene_metag_FILTERED.tsv"),sep='\t') %>% mutate(seqtype = 'METAGENOMICS') %>% mutate(regressionclass = 'Nasal')   %>% mutate(dset = 'GENE-CATALOG')
      regdata[[68]] = read.csv(paste0("/Users/bradentemp/Dropbox (Mason Lab)/i4/revisions/regressionoutput/gene/regression_output_nasal",ii,"_merged_gene_metat_FILTERED.tsv"),sep='\t') %>% mutate(seqtype = 'METATRANSCRIPTOMICS') %>% mutate(regressionclass = 'Nasal') %>% mutate(dset = 'GENE-CATALOG')
      regdata[[69]] = read.csv(paste0("/Users/bradentemp/Dropbox (Mason Lab)/i4/revisions/regressionoutput/gene/regression_output_skin",ii,"_merged_gene_metag_FILTERED.tsv"),sep='\t') %>% mutate(seqtype = 'METAGENOMICS') %>% mutate(regressionclass = 'Skin') %>% mutate(dset = 'GENE-CATALOG')
      regdata[[70]] = read.csv(paste0("/Users/bradentemp/Dropbox (Mason Lab)/i4/revisions/regressionoutput/gene/regression_output_skin",ii,"_merged_gene_metat_FILTERED.tsv"),sep='\t') %>% mutate(seqtype = 'METATRANSCRIPTOMICS') %>% mutate(regressionclass = 'Skin') %>% mutate(dset = 'GENE-CATALOG')
      regdata[[71]] = read.csv(paste0("/Users/bradentemp/Dropbox (Mason Lab)/i4/revisions/regressionoutput/gene/regression_output_skinseparates",ii,"_merged_gene_metag_FILTERED.tsv"),sep='\t') %>% mutate(seqtype = 'METAGENOMICS') %>% mutate(regressionclass = 'Skin_Components') %>% mutate(dset = 'GENE-CATALOG')
      regdata[[72]] = read.csv(paste0("/Users/bradentemp/Dropbox (Mason Lab)/i4/revisions/regressionoutput/gene/regression_output_skinseparates",ii,"_merged_gene_metat_FILTERED.tsv"),sep='\t') %>% mutate(seqtype = 'METATRANSCRIPTOMICS') %>% mutate(regressionclass = 'Skin_Components') %>% mutate(dset = 'GENE-CATALOG')
      }

      ### OVERALL
      #regdata[[49]] = read.table("i4_data_packet/differential_abundance_or_expression_analysis/regression_output_overall_bacterial_dna_",secondrankbecausesloppynaming,"_xtree_.05_.0025_GTDB_r207.tsv"),sep='\t') %>% mutate(seqtype = 'METAGENOMICS') %>% mutate(regressionclass = 'Overall')   %>% mutate(dset = 'GTDB')%>% mutate(term = if_else(term == 'Timepoint_RecodePRE-LAUNCH','PRE-LAUNCH --- OVERALL',if_else(term == 'Timepoint_RecodePOST-LAUNCH','POST-LAUNCH --- OVERALL',term)))
      #regdata[[50]] = read.table("i4_data_packet/differential_abundance_or_expression_analysis/regression_output_overall_bacterial_rna_",secondrankbecausesloppynaming,"_xtree_.05_.0025_GTDB_r207.tsv"),sep='\t') %>% mutate(seqtype = 'METATRANSCRIPTOMICS') %>% mutate(regressionclass = 'Overall') %>% mutate(dset = 'GTDB')%>% mutate(term = if_else(term == 'Timepoint_RecodePRE-LAUNCH','PRE-LAUNCH --- OVERALL',if_else(term == 'Timepoint_RecodePOST-LAUNCH','POST-LAUNCH --- OVERALL',term)))
      #regdata[[51]] = read.table("i4_data_packet/differential_abundance_or_expression_analysis/regression_output_overall_bacterial_dna_",secondrankbecausesloppynaming,"_xtree_.05_.0025_ASSEMBLIES.tsv"),sep='\t') %>% mutate(seqtype = 'METAGENOMICS') %>% mutate(regressionclass = 'Overall') %>% mutate(dset = 'ASSEMBLED-BACTERIA')%>% mutate(term = if_else(term == 'Timepoint_RecodePRE-LAUNCH','PRE-LAUNCH --- OVERALL',if_else(term == 'Timepoint_RecodePOST-LAUNCH','POST-LAUNCH --- OVERALL',term)))
      #regdata[[52]] = read.table("i4_data_packet/differential_abundance_or_expression_analysis/regression_output_overall_bacterial_rna_",secondrankbecausesloppynaming,"_xtree_.05_.0025_ASSEMBLIES.tsv"),sep='\t') %>% mutate(seqtype = 'METATRANSCRIPTOMICS') %>% mutate(regressionclass = 'Overall') %>% mutate(dset = 'ASSEMBLED-BACTERIA')%>% mutate(term = if_else(term == 'Timepoint_RecodePRE-LAUNCH','PRE-LAUNCH --- OVERALL',if_else(term == 'Timepoint_RecodePOST-LAUNCH','POST-LAUNCH --- OVERALL',term)))
      #regdata[[53]] = read.csv("i4_data_packet/differential_abundance_or_expression_analysis/regression_output_overall_viral_dna_",secondrankbecausesloppynaming,"_xtree_.01_.0.005_VIRAL_GENBANK.tsv"),sep='\t') %>% mutate(seqtype = 'METAGENOMICS') %>% mutate(regressionclass = 'Overall') %>% mutate(dset = 'GENBANK-VIRUSES')%>% mutate(term = if_else(term == 'Timepoint_RecodePRE-LAUNCH','PRE-LAUNCH --- OVERALL',if_else(term == 'Timepoint_RecodePOST-LAUNCH','POST-LAUNCH --- OVERALL',term)))
      #regdata[[54]] = read.csv("i4_data_packet/differential_abundance_or_expression_analysis/regression_output_overall_viral_rna_",secondrankbecausesloppynaming,"_xtree_.01_.0.005_VIRAL_GENBANK.tsv"),sep='\t') %>% mutate(seqtype = 'METATRANSCRIPTOMICS') %>% mutate(regressionclass = 'Overall') %>% mutate(dset = 'GENBANK-VIRUSES')%>% mutate(term = if_else(term == 'Timepoint_RecodePRE-LAUNCH','PRE-LAUNCH --- OVERALL',if_else(term == 'Timepoint_RecodePOST-LAUNCH','POST-LAUNCH --- OVERALL',term)))
      #regdata[[55]] = read.csv("i4_data_packet/differential_abundance_or_expression_analysis/regression_output_overall_viral_dna_",secondrankbecausesloppynaming,"_xtree_.01_.0.005_ASSEMBLIES.tsv"),sep='\t') %>% mutate(seqtype = 'METAGENOMICS') %>% mutate(regressionclass = 'Overall') %>% mutate(dset = 'ASSEMBLED-VIRUSES')%>% mutate(term = if_else(term == 'Timepoint_RecodePRE-LAUNCH','PRE-LAUNCH --- OVERALL',if_else(term == 'Timepoint_RecodePOST-LAUNCH','POST-LAUNCH --- OVERALL',term)))
      #regdata[[56]] = read.csv("i4_data_packet/differential_abundance_or_expression_analysis/regression_output_overall_viral_rna_",secondrankbecausesloppynaming,"_xtree_.01_.0.005_ASSEMBLIES.tsv"),sep='\t') %>% mutate(seqtype = 'METATRANSCRIPTOMICS') %>% mutate(regressionclass = 'Overall') %>% mutate(dset = 'ASSEMBLED-VIRUSES')%>% mutate(term = if_else(term == 'Timepoint_RecodePRE-LAUNCH','PRE-LAUNCH --- OVERALL',if_else(term == 'Timepoint_RecodePOST-LAUNCH','POST-LAUNCH --- OVERALL',term)))
      #regdata[[57]] = read.csv("i4_data_packet/differential_abundance_or_expression_analysis/regression_output_overall_bacterial-viral_dna_",secondrankbecausesloppynaming,"_kraken2_none_KRAKEN2_REFSEQDB.tsv"),sep='\t') %>% mutate(seqtype = 'METAGENOMICS') %>% mutate(regressionclass = 'Overall') %>% mutate(dset = 'KRAKEN2')%>% mutate(term = if_else(term == 'Timepoint_RecodePRE-LAUNCH','PRE-LAUNCH --- OVERALL',if_else(term == 'Timepoint_RecodePOST-LAUNCH','POST-LAUNCH --- OVERALL',term)))
      #regdata[[58]] = read.csv("i4_data_packet/differential_abundance_or_expression_analysis/regression_output_overall_bacterial-viral_rna_",secondrankbecausesloppynaming,"_kraken2_none_KRAKEN2_REFSEQDB.tsv"),sep='\t') %>% mutate(seqtype = 'METATRANSCRIPTOMICS') %>% mutate(regressionclass = 'Overall') %>% mutate(dset = 'KRAKEN2') %>% mutate(term = if_else(term == 'Timepoint_RecodePRE-LAUNCH','PRE-LAUNCH --- OVERALL',if_else(term == 'Timepoint_RecodePOST-LAUNCH','POST-LAUNCH --- OVERALL',term)))
      #regdata[[59]] = read.csv("i4_data_packet/differential_abundance_or_expression_analysis/regression_output_overall_merged_gene_metag.tsv"),sep=' ') %>% mutate(seqtype = 'METAGENOMICS') %>% mutate(regressionclass = 'Overall') %>% mutate(dset = 'GENE-CATALOG')
      #regdata[[60]] = read.csv("i4_data_packet/differential_abundance_or_expression_analysis/regression_output_overall_merged_gene_metat.tsv"),sep=' ') %>% mutate(seqtype = 'METATRANSCRIPTOMICS') %>% mutate(regressionclass = 'Overall') %>% mutate(dset = 'GENE-CATALOG')
      
      regression_output = bind_rows(regdata)
      
      regression_output$term = gsub('Timepoint_Recode','',regression_output$term)
      regression_output$term = gsub(':','',regression_output$term)
      regression_output$term = gsub(':Body.Location','',regression_output$term)
      regression_output$term = gsub('PRE-LAUNCH','PRE-LAUNCH --- ',regression_output$term)
      regression_output$term = gsub('POST-LAUNCH','POST-LAUNCH --- ',regression_output$term)
      regression_output$term = gsub('isskin','Skin',regression_output$term)
      regression_output$term = gsub('isoral','Oral',regression_output$term)
      regression_output$term = gsub('nape','Nape of Neck',regression_output$term)
      regression_output$term = gsub('postauric','Post-Auricular',regression_output$term)
      regression_output$term = gsub('bb','Belly Button',regression_output$term)
      regression_output$term = gsub('gc','Gluteal Crease',regression_output$term)
      regression_output$term = gsub('Tzone','T-Zone',regression_output$term)
      regression_output$term = gsub('web','Toe Web Space',regression_output$term)
      regression_output$term = gsub('fore','Forearm',regression_output$term)
      regression_output$term = gsub(' ---  --- ',' --- ',regression_output$term)
      regression_output$term = gsub('isnasal','Nasal',regression_output$term)
      regression_output = regression_output %>% mutate(term = if_else(term == 'PRE-LAUNCH --- ','PRE-LAUNCH',term,))
      regression_output = regression_output %>% mutate(term = if_else(term == 'POST-LAUNCH --- ','POST-LAUNCH',term,))
      regression_output = regression_output %>% mutate(launch = if_else(grepl('PRE',term),'PRE-LAUNCH','POST-LAUNCH')) %>% mutate(direction = if_else(estimate>0,'positive','negative')) %>% filter(term != 'POST-LAUNCH --- OVERALL')%>% filter(term != 'PRE-LAUNCH --- OVERALL')
      #temp0 =regression_output %>% filter(grepl('---',term)) %>% mutate(regclass2 = strsplit(as.character(term),' --- ') %>% map_chr(2))
      #temp1 = regression_output %>% filter(!grepl('---',term),grepl('LAUNCH',term)) %>% mutate(regclass2 = 'Main effects (pre/post launch)')
      #temp2 = regression_output %>% filter(!grepl('---',term),!grepl('LAUNCH',term)) %>% mutate(regclass2 = paste0('Main effects (',term,')'))
      #regression_output = bind_rows(temp0,temp1,temp2) %>% group_by(regclass2) %>% mutate(BH_adjusted = p.adjust(p.value,method='BH'),BY_adjusted = p.adjust(p.value,method='BY'),BONFERRONI_adjusted = p.adjust(p.value,method='bonferroni')) %>% ungroup

      saveRDS(regression_output,paste0('~/Dropbox (Mason Lab)/i4/revisions/regressionoutput/parsed/regression_data_timetrends_',secondrankbecausesloppynaming,ii,'_',iii,'.rds'))

      regression_output = regression_output  %>% filter(grepl('---',term)) 
      regression_output$term = factor(regression_output$term,levels = c("PRE-LAUNCH --- Armpit","PRE-LAUNCH --- Belly Button","PRE-LAUNCH --- Forearm", "PRE-LAUNCH --- Gluteal Crease","PRE-LAUNCH --- Nape of Neck","PRE-LAUNCH --- Nasal","PRE-LAUNCH --- Post-Auricular","PRE-LAUNCH --- T-Zone","PRE-LAUNCH --- Toe Web Space","PRE-LAUNCH --- Skin","PRE-LAUNCH --- Oral","POST-LAUNCH --- Armpit","POST-LAUNCH --- Belly Button", "POST-LAUNCH --- Forearm","POST-LAUNCH --- Gluteal Crease" ,"POST-LAUNCH --- Nape of Neck", "POST-LAUNCH --- Nasal", "POST-LAUNCH --- Post-Auricular","POST-LAUNCH --- T-Zone", "POST-LAUNCH --- Toe Web Space","POST-LAUNCH --- Skin","POST-LAUNCH --- Oral"))
      regression_output = regression_output %>% mutate(launch = if_else(grepl('PRE',term),'PRE-LAUNCH','POST-LAUNCH')) %>% mutate(direction = if_else(estimate>0,'positive','negative')) %>% filter(term != 'POST-LAUNCH --- OVERALL')%>% filter(term != 'PRE-LAUNCH --- OVERALL')
      
      a =regression_output %>% filter(launch == 'PRE-LAUNCH') %>% mutate(regclass2 = strsplit(as.character(term),' --- ') %>% map_chr(2))  %>% dplyr::rename(PRE=launch,adj_pre = BH_adjusted,est_pre = estimate) %>% select(yvar,dset,seqtype,regclass2,PRE,adj_pre,est_pre)
      #a2 =regression_output %>% filter(launch == 'PRE-LAUNCH')%>% filter(term == 'PRE-LAUNCH --- OVERALL') %>% mutate(regclass2 = 'Overall') %>% mutate(mergeterm = paste(yvar,dset,seqtype,regclass2)) %>% dplyr::rename(PRE=launch,adj_pre = BH_adjusted,est_pre = estimate) %>% select(yvar,dset,seqtype,regclass2,PRE,adj_pre,est_pre)
      #a3 =regression_output %>% filter(regressionclass == 'Overall') %>% filter(launch == 'PRE-LAUNCH')%>% filter(term == 'PRE-LAUNCH') %>% mutate(regclass2 = 'Overall') %>% mutate(mergeterm = paste(yvar,dset,seqtype,regclass2)) %>% dplyr::rename(PRE=launch,adj_pre = BH_adjusted,est_pre = estimate) %>% select(yvar,dset,seqtype,regclass2,PRE,adj_pre,est_pre)
      
     #a = bind_rows(a,a3)
      
      b =regression_output %>% filter(launch == 'POST-LAUNCH') %>% filter(grepl('---',term)) %>% mutate(regclass2 = strsplit(as.character(term),' --- ') %>% map_chr(2))%>% dplyr::rename(POST=launch,adj_post = BH_adjusted,est_post = estimate) %>% select(yvar,dset,seqtype,regclass2,POST,adj_post,est_post)
      #b2 =regression_output %>% filter(launch == 'POST-LAUNCH') %>% filter(term == 'POST-LAUNCH --- OVERALL') %>% mutate(regclass2 = 'Overall') %>% mutate(mergeterm = paste(yvar,dset,seqtype,regclass2)) %>% dplyr::rename(POST=launch,adj_post = BH_adjusted,est_post = estimate) %>% select(yvar,dset,seqtype,regclass2,POST,adj_post,est_post)
      #b3 =regression_output %>% filter(regressionclass == 'Overall')%>% filter(launch == 'POST-LAUNCH') %>% filter(term == 'POST-LAUNCH') %>% mutate(regclass2 = 'Overall') %>% mutate(mergeterm = paste(yvar,dset,seqtype,regclass2)) %>% dplyr::rename(POST=launch,adj_post = BH_adjusted,est_post = estimate) %>% select(yvar,dset,seqtype,regclass2,POST,adj_post,est_post)
      #b = bind_rows(b,b3)
      
      c = inner_join(a,b)
      
      c = c %>% filter(adj_pre<0.05 |  adj_post<0.05)
      
      ### MODERATE CONSERVATISM
      #c = c %>% mutate(timetrend = 'NONE')%>% mutate(timetrend = if_else(est_pre<0 & est_post>0 & adj_pre <0.05 |est_pre<0 & est_post>0 & adj_post <0.05,'Increased in/after flight',timetrend)) %>% mutate(timetrend = if_else(est_pre>0 & est_post<0 & adj_pre <0.05 | est_pre>0 & est_post<0 & adj_post <0.05,'Decreased in/after flight',timetrend)) %>% mutate(timetrend = if_else(est_pre<0 & est_post<0 & adj_pre<0.05 | est_pre<0 & est_post<0 & adj_post <0.05,'Increased in flight',timetrend)) %>% mutate(timetrend = if_else(est_pre>0 & est_post>0 & adj_pre<0.05 | est_pre>0 & est_post>0  & adj_post <0.05,'Decreased in flight',timetrend)) 
      
      ### HIGH CONSERVATISM
      #c = c %>% mutate(timetrend = 'NONE')%>% mutate(timetrend = if_else(est_pre<0 & est_post>0 & adj_pre <0.05 |est_pre<0 & est_post>0 & adj_post <0.05,'Increased in/after flight',timetrend)) %>% mutate(timetrend = if_else(est_pre>0 & est_post<0 & adj_pre <0.05 | est_pre>0 & est_post<0 & adj_post <0.05,'Decreased in/after flight',timetrend)) %>% mutate(timetrend = if_else(est_pre<0 & est_post<0 & adj_pre<0.05 & adj_post <0.05,'Increased in flight',timetrend)) %>% mutate(timetrend = if_else(est_pre>0 & est_post>0 & adj_pre<0.05 & adj_post <0.05,'Decreased in flight',timetrend)) 
      
      ### BOTH
      c = c %>% mutate(timetrendrank = 'NONE')%>% mutate(timetrendrank = if_else(est_pre<0 & est_post>0 & adj_pre <0.05 |est_pre<0 & est_post>0 & adj_post <0.05,'Persistent increase',timetrendrank)) %>% mutate(timetrendrank = if_else(est_pre>0 & est_post<0 & adj_pre <0.05 | est_pre>0 & est_post<0 & adj_post <0.05,'Persistent decrease',timetrendrank)) %>% mutate(timetrendrank = if_else(est_pre<0 & est_post<0 & adj_pre<0.05 | est_pre<0 & est_post<0 & adj_post <0.05,'Transient increase -- LP',timetrendrank)) %>% mutate(timetrendrank = if_else(est_pre>0 & est_post>0 & adj_pre<0.05 | est_pre>0 & est_post>0  & adj_post <0.05,'Transient decrease -- LP',timetrendrank)) %>% mutate(timetrendrank = if_else(est_pre<0 & est_post<0 & adj_pre<0.05 & adj_post <0.05,'Transient increase',timetrendrank)) %>% mutate(timetrendrank = if_else(est_pre>0 & est_post>0 & adj_pre<0.05 & adj_post <0.05,'Transient decrease',timetrendrank)) 
      c = c %>% mutate(timetrend = gsub(' -- LP','',timetrendrank))
      
      fortrendline = c %>% filter(regclass2 != 'OVERALL') %>% mutate(mergeid = paste(yvar,dset,seqtype,regclass2))
      
      write.csv(fortrendline,paste0('~/Dropbox (Mason Lab)/i4/revisions/regressionoutput/parsed/regression_data_timetrends_',secondrankbecausesloppynaming,ii,'_',iii,'.csv'))
     # write.csv(bind_rows(regdata),paste0('~/Dropbox (Mason Lab)/i4/revisions/regressionoutput/parsed/regression_data_',secondrankbecausesloppynaming,ii,'_',iii,'.csv'))
    }
  }
}
