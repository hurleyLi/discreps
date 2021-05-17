library("LOLA")
library("ggplot2")
library("dplyr")
options(warn = -1)

#######################################################################
########## Function to perform enrichment analyis using LOLA ##########
#######################################################################

enrich = function(assembly, type) {
  # Load pre-assembled loci for tested genomic elements
  regionDB = loadRegionDB(paste0("reference/LOLA/", assembly))
  
  # Overlapping DISCEPs with various genomic elements
  regionSetA = readBed(paste0("inputs/DISCREPs_",assembly,"_10kb.",type,".bed"))
  universe = readBed(paste0("inputs/total.",assembly,".bed.binCount.filtered"))
  
  locResults = runLOLA(regionSetA, universe, regionDB, cores=1)
  locResults$group = ifelse(type == "overlap", type, paste(assembly,type))
  return(locResults)
}

###########################################################
########## Function to perform Fisher exact test ##########
###########################################################

fisher_test = function(df, columns) {
  df$fisher = apply(df, 1, function(x) fisher.test( 
    matrix( c(
      as.numeric(x[columns[1]]),
      as.numeric(x[columns[2]]),
      as.numeric(x[columns[3]]),
      as.numeric(x[columns[4]])
    ), nrow = 2),
  ))

  df$pValue = apply(df, 1, function(x) x['fisher']$fisher$p.value)
  df$oddsRatio = apply(df, 1, function(x) x['fisher']$fisher$estimate)
  df$CI95_lower = apply(df, 1, function(x) x['fisher']$fisher$conf.int[1])
  df$CI95_higher = apply(df, 1, function(x) x['fisher']$fisher$conf.int[2])
  df$log_OR = log2(df$oddsRatio)
  df$log_CI_higher = log2(df$CI95_higher)
  df$log_CI_lower = log2(df$CI95_lower)
  return(df)
}

##############################################################
########## Enrichment analyses for each DISCREP set ##########
##############################################################

allResults = NULL
for (assembly in c('GRCh37','GRCh38')) {
  for (type in c('overlap','non-overlap')) {
    allResults = rbind(allResults, enrich(assembly, type))
  }
}

# Perform statistical tests
df_use = allResults[,c('support','b','c','d','description','group')] # b,c,d are the default columns names from LOLA
df_use_ready = df_use %>%
  group_by(group, description) %>% 
  summarize(support = sum(support), b = sum(b), c = sum(c), d = sum(d))
df_use_ready = fisher_test(df_use_ready, c('support','b','c','d')) # these are the column names from LOLA by default

# Order results by the odds ratio in the overlap DISCREPs set
df_use_ready_overlap = df_use_ready[df_use_ready$group == "overlap",]
df_use_ready$description = factor(df_use_ready$description, 
  levels = df_use_ready_overlap[order(df_use_ready_overlap$oddsRatio, decreasing=TRUE),]$description )
df_use_ready$group = factor(df_use_ready$group,
  levels = c('GRCh38 non-overlap','GRCh37 non-overlap','overlap'))


df_use_ready$qValue = p.adjust(df_use_ready$pValue)
df_use_ready$sig = apply(df_use_ready, 1, function(x) ifelse(as.numeric(x['qValue']) < 0.01, 'sig', "xno"))

######################################################################################
########## Plotting of the enrichment analyses results for each DISCREP set ##########
######################################################################################
p = ggplot(data = df_use_ready,
  aes(x = group,  y = log_OR, ymin = log_CI_lower, ymax = log_CI_higher)) +
  geom_pointrange(aes(col=group, ), linetype = 0)+
  geom_hline(aes(fill=group),yintercept =0, linetype=2)+
  xlab('')+ ylab("Odds Ratio and 95% CI (log scale)")+
  geom_errorbar(aes(ymin=log_CI_lower, ymax=log_CI_higher,col=group, linetype = sig, ),width=0.4,cex=1)+
  facet_wrap(~description, strip.position="left",nrow=9,  ) +
  theme(plot.title=element_text(size=16,face="bold"),
    plot.margin = margin(5, 5, 20, 5),
    axis.text.y=element_blank(),
    axis.text.x=element_text(face="bold", size = 15),
    axis.ticks.y=element_blank(),
    axis.title=element_text(size=12,face="bold"),
    axis.title.x=element_text(face="bold", size = 20, hjust=-0.4, vjust = -0.4),
        
    strip.text.y = element_text(hjust=0.5,vjust = 0.5,angle=180,face="bold"),
    strip.text.y.left = element_text(angle = 0),
    strip.text = element_text(size = 15),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    legend.text=element_text(size=22),
    legend.title=element_blank(),
    legend.position = c(0.7, 0.19))+
  scale_y_continuous(breaks=c( log2(1), log2(2), log2(5), log2(10), log2(20),log2(50), log2(100),log2(200)),
    labels = c('1', '2','5','10','20','50','100','200'))+
  scale_color_manual(values=c("#311F23", "#6666DD", "#DD8452" ),
    breaks=c("overlap", "GRCh37 non-overlap", "GRCh38 non-overlap"),
    labels=c("Overlapping DISCREPs", "Unique to GRCh37", "Unique to GRCh38"))+
  scale_linetype_manual(values=c("solid", "dashed"))+
  guides(linetype=FALSE)+
  coord_flip()
p

##################################################################
########## Compare between different groups of DISCREPs ##########
##################################################################

# construct data frame for pairwised comparisons
df_use_compare = df_use %>%
  group_by(group, description) %>% 
  summarize(support = sum(support), not_support = sum(c))

df_unique_hg19 = df_use_compare[df_use_compare$group == "GRCh37 non-overlap", ]
df_unique_hg38 = df_use_compare[df_use_compare$group == "GRCh38 non-overlap", ]
df_overlap = df_use_compare[df_use_compare$group == "overlap", ]

df_overlap_hg19 = merge(df_overlap, df_unique_hg19, by = 'description')
df_overlap_hg38 = merge(df_overlap, df_unique_hg38, by = 'description')
df_hg19_hg38 = merge(df_unique_hg19, df_unique_hg38, by = 'description')

fisherTests = lapply(list(df_overlap_hg19, df_overlap_hg38, df_hg19_hg38), 
  fisher_test, columns = c('support.x','not_support.x','support.y','not_support.y'))
df_overlap_hg19 = fisherTests[[1]] 
df_overlap_hg38 = fisherTests[[2]]
df_hg19_hg38 = fisherTests[[3]]

df_hg19_hg38$group = 'hg19_vs_hg38'
df_overlap_hg38$group = 'overlap_vs_hg38'
df_overlap_hg19$group = 'overlap_vs_hg19'

# combine different groups together
selected_features = c(
  'segmental duplication', 'assembly problems', 'fix patch sequences',  
  'alternate haplotype' ,'genome assemblies difference', 'gaps in assembly')

df_compare_combined = do.call("rbind", list(df_overlap_hg19, df_overlap_hg38, df_hg19_hg38))
df_compare_combined = df_compare_combined[df_compare_combined$description %in% selected_features,]
df_compare_combined$description = factor(df_compare_combined$description, levels = rev(selected_features))
df_compare_combined$qValue = p.adjust(df_compare_combined$pValue)
df_compare_combined$sig = apply(df_compare_combined, 1, function(x) ifelse(as.numeric(x['qValue']) < 0.01, 'sig', "xno"))
df_compare_combined$group = factor(df_compare_combined$group, 
  labels = c('Overlap vs GRCh37 Unique', 'GRCh37 Unique vs\nGRCh38 Unique', 'Overlap vs GRCh38 Unique'),
  levels = c('overlap_vs_hg19','hg19_vs_hg38','overlap_vs_hg38'))

######################################################################################
########## Plotting of the pairwised comparisons between each DISCREP set ############
######################################################################################

p = ggplot(data = df_compare_combined,
  aes(x = description,  y = log_OR, ymin = log_CI_lower, ymax = log_CI_higher)) +
  geom_pointrange(aes(col=sig), linetype = 0)+
  geom_hline(aes(fill=sig),yintercept =0, linetype=2)+
  xlab('')+ ylab("Odds Ratio and 95% CI (log scale)")+
  geom_errorbar(aes(ymin=log_CI_lower, ymax=log_CI_higher,col=sig, linetype = 'solid'),width=0.2,cex=1)+
  facet_wrap(~group, strip.position="top",nrow=1) +
  theme(plot.title=element_text(size=16,face="bold"),
    axis.text.x=element_text(face="bold", size = 20, vjust = 1),
    axis.text.y=element_text(face="bold", size = 18),
    axis.title.x=element_text(face="bold", size = 20, vjust = 0),
    strip.text.x = element_text(hjust=0.5,vjust = 0.5,angle=180,face="bold"),
    strip.text.x.top = element_text(angle = 0, size = 20),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    legend.position = "none",
  )+
  scale_y_continuous(breaks=c(log2(0.5), log2(1), log2(2), log2(5), log2(10)),
    labels = c('0.5', '1', '2','5','10'))+
  scale_color_manual(values=c("black", "grey70"))+
  guides(linetype=FALSE)+
  coord_flip()
p

################################################
########## Output results in tables ############
################################################

output = df_use_ready[c('group','description','pValue','qValue', 'oddsRatio','CI95_lower','CI95_higher')]
write.table(output, 'results/LOLA_enrichment.tsv', sep = '\t', row.names = F, quote = F)


