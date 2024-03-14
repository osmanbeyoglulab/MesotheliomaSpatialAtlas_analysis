library(tripack)
library(stats)
library(data.table)

library(foreach)
library(doParallel)
library(tidyr)



rm(list=ls())
dir = dirname(rstudioapi::getSourceEditorContext()$path); print(dir)  # [1] "/Users/lees18/OneDrive - UPMC/S86_Covid19_scOmics_Xiaosjun/3_SeuratObject_ClusteringByPatient"
setwd(dir)

# args = commandArgs(trailingOnly=TRUE)
# casetype = args[1]

setwd(system("pwd", intern = T) )

getLogsOddsRatioMtx <- function(x, y, celltypes, ptype1, ptype2){
  
  
  typesetsub = sort(unique(celltypes))
  
  
  tr = tri.mesh(x, y, duplicate = "error")
  
  nbs = neighbours(tr)
  
  counts = data.frame(matrix(0, ncol = length(nbs), nrow = length(typesetsub))) #typeset X nodes
  rownames(counts) = typesetsub
  colnames(counts) = seq(length(nbs))
  
  
  for (i in seq(length(nbs))) {
    for (j in nbs[[i]]) {
      counts[celltypes[j],i] = counts[celltypes[j],i] + 1
    }
  }
  
  correl =data.frame(matrix(0, ncol = length(typesetsub), nrow = length(typesetsub))) # celltype X celltype
  rownames(correl) = typesetsub
  colnames(correl) = typesetsub
  
  for (type1 in c(ptype1)){
    for (type2 in c(ptype2)){
      # for (type1 in typesetsub){
      #   for (type2 in typesetsub){
      
      countCellsOfThisType = 0
      countCellsOtherType = 0
      
      TotalInteractionsOfThisType = 0
      TotalInteractionsOfOtherType = 0
      
      TotalInteractions = 0
      
      InteractionsBtw = 0
      
      for (k in seq(length(nbs))) { #loop of nodes
        
        if ( celltypes[k] == type1) {
          
          countCellsOfThisType = countCellsOfThisType + 1
          
          for (l in typesetsub) { # loop of cell types
            TotalInteractionsOfThisType =  TotalInteractionsOfThisType + counts[l,k];
          }
          
          InteractionsBtw = InteractionsBtw +  counts[type2,k];
        }
        
        if (celltypes[k] == type2) {
          countCellsOtherType = countCellsOtherType + 1
          
          for (l in typesetsub) { # loop of cell types
            TotalInteractionsOfOtherType =  TotalInteractionsOfOtherType + counts[l,k];
          }
        }
        
        for (l in typesetsub) { # loop of cell types
          TotalInteractions =  TotalInteractions + counts[l,k];
        }
      }
      
      
      seq01 = seq(0.01, 0.99, by = 0.01)
      # bd = dbeta(seq01, InteractionsBtw , (TotalInteractionsOfThisType - InteractionsBtw))
      bd = dbeta(seq01, InteractionsBtw , (TotalInteractionsOfThisType + TotalInteractionsOfOtherType - InteractionsBtw))
      oddsRatio = mean(bd) / ((countCellsOtherType / length(nbs)) * (countCellsOfThisType / length(nbs) ))
      # oddsRatio = mean(bd) / ((TotalInteractionsOfThisType / TotalInteractions) * (TotalInteractionsOfOtherType / TotalInteractions ))
      
      correl[type1, type2] = log(oddsRatio+.Machine$double.eps)
    }
  }
  
  return(correl)
}


########################################################################################

# setwd(system("pwd", intern = T) )

# args = commandArgs(trailingOnly=TRUE)

args = commandArgs(trailingOnly=TRUE)
ptype1 = args[1]
ptype2 = args[2]

tissue = "Tumor" #NA #"Tumor" #NA # 'Stroma' #NA # args[3]

casetype = args[3]
subtype = args[4]
# print(casetype)
# print(tissue)


df = fread("../Data_folder/TMA1/panel1_data_cleaned_withCelltypes_merged.csv")
df = as.data.frame(df)

cores = unique(df[df$'Classification'=='Malignant', "Annotation ID"])

pair_scores = c()
cores_sel = c()
n_nullCores = 0
for (core in cores){
  
  if (is.na(tissue) ) {
    df_sub = df[(df$'Annotation ID'==core), c("Cell X Position", "Cell Y Position", "phenotype_combined")]
  }else{
    df_sub = df[(df$'Annotation ID'==core) & (df$'Tissue Category'== tissue), c("Cell X Position", "Cell Y Position", "phenotype_combined")]
  }
  
  if (! is.na(casetype)){
    df_sub = df_sub[df_sub$CaseType == casetype, ]
  }
  
  if (! is.na(subtype)){
    df_sub = df_sub[df_sub$subtype == subtype, ]
  }
  
  
  if (dim(df_sub)[1]<5) {next}
  
  
  
  x = df_sub$"Cell X Position"
  y = df_sub$"Cell Y Position"
  celltypes = df_sub$phenotype_combined
  
  oddsRatios <- getLogsOddsRatioMtx(x,y, celltypes, ptype1, ptype2)
  if (is.null( oddsRatios[ptype1, ptype2])){
    print("oddsRatio is NULL")
    print(unique(celltypes))
    n_nullCores = n_nullCores+1
  }
  pair_scores = c(pair_scores, oddsRatios[ptype1, ptype2])
  cores_sel = c(cores_sel, core)
  
}

df_score = as.data.frame(cbind(cores_sel, pair_scores))
colnames(df_score) = c("core ID", "pair_score")

filename =  paste0("./output/cellContact/",ptype1, "_", ptype2) #"_scores.csv"
if (! is.na(tissue)){
  filename = paste0(filename, "_", tissue)
}
if (! is.na(casetype)){
  filename = paste0(filename, "_", casetype)
}      
if (! is.na(subtype)){
  filename = paste0(filename, "_", subtype)
}      

filename = paste0(file, "_scores.csv")
write.csv(df_score, file=paste0(filename,"_scores.csv"), row.names=FALSE)

          
# ########################################################################################
# if (is.na(tissue) & is.na(casetype)){
#   df_score = read.csv( paste0("./output/cellContact/",ptype1, "_", ptype2, "_scores.csv"), check.names = FALSE)
# }else if (is.na(casetype)){
#   df_score = read.csv( paste0("./output/cellContact/",ptype1, "_", ptype2, "_", tissue,   "_scores.csv"), check.names = FALSE)
# }else if (is.na(tissue)){
#   df_score = read.csv( paste0("./output/cellContact/",ptype1, "_", ptype2, "_", casetype,   "_scores.csv"), check.names = FALSE)
# }else{
#   df_score = read.csv( paste0("./output/cellContact/",ptype1, "_", ptype2, "_", casetype,  "_", tissue,   "_scores.csv"), check.names = FALSE)
# }

df_core_feature = read.csv("./output/cellContact/core_features_allMarkers.csv", check.names = FALSE)

df_out = merge(df_core_feature, df_score, by = 'core ID')
df_out <- na.omit(df_out)

write.csv(df_out, file = paste0(filename,"_scores_wCoreInfo.csv"), row.names=FALSE )

# 
# if (is.na(tissue) ){
#   write.csv(df_out, file = paste0("./output/cellContact/",ptype1, "_", ptype2, "_scores_wCoreInfo.csv"), row.names=FALSE )
# }else {
#   write.csv(df_out, file = paste0("./output/cellContact/",ptype1, "_", ptype2, "_", tissue,  "_scores_wCoreInfo.csv"), row.names=FALSE )
# }
# 


# #######################################################################################
# # group analysis, ridge density plot
# ptype1 = "CD4"
# ptype2 = "CK"
# if (is.na(tissue)){
#   df_score = read.csv( paste0("./output/cellContact/",ptype1, "_", ptype2, "_scores_wCoreInfo.csv"), check.names = FALSE)
# }else{
#   df_score = read.csv( paste0("./output/cellContact/",ptype1, "_", ptype2, "_", tissue ,"_scores_wCoreInfo.csv"), check.names = FALSE)
# }
# # df_score = df_score[(df_score[paste0(ptype2, ' percent')] > 0.05) & (df_score$'tumor density' > 400), ]
# df_score = df_score[(df_score[paste0(ptype2, ' percent')] > 0.05) & (df_score$'CK percent' > 0.2), ]
# markerSep = c(200,800)
# 
# groupname = "BAP1 density"
# highCond = df_score[groupname] >= markerSep[2]
# df_score[highCond, 'BAP1 status'] = "BAP1+"
# lowCond = df_score[groupname] < markerSep[1]
# df_score[lowCond, 'BAP1 status'] = "BAP1-"
# 
# groupname = "NF2 density"
# highCond = df_score[groupname] >= markerSep[2]
# df_score[highCond, 'NF2 status'] = "NF2+"
# lowCond = df_score[groupname] < markerSep[1]
# df_score[lowCond, 'NF2 status'] = "NF2-"
# 
# groupname = "MTAP density"
# highCond = df_score[groupname] >= markerSep[2]
# df_score[highCond, 'MTAP status'] = "MTAP+"
# lowCond = df_score[groupname] < markerSep[1]
# df_score[lowCond, 'MTAP status'] = "MTAP-"
# 
# df_score[is.na(df_score)]=""
# df_score['group'] = paste0(df_score$`BAP1 status`, df_score$`NF2 status`, df_score$`MTAP status`)
# 
# df_score['group'] = df_score$`BAP1 status`
# 
# groups = table(df_score['group'])
# 
# library(ggplot2)
# # ridge density plot
# path_figure = "./output/cellContact/figure"
# if (is.na(tissue) & is.na(casetype) ){
#   title = paste0("CK~", ptype2, " interactions among marker groups in all region")
# }else if (is.na(casetype)){
#   title = paste0("CK~", ptype2, " interactions among marker groups in ",tolower(tissue)," region")
# }else if (is.na(tissue)){
#   title = paste0(casetype, " CK~", ptype2, " interactions among marker groups")
# }else{
#   title = paste0(casetype, " CK~", ptype2, " interactions among marker groups in ",tolower(tissue)," region")
# }
# png(file= file.path(path_figure, paste0("Ridge ",title,".png")),width=600, height=200)
# g = ggplot(df_score, aes(x = pair_score,
#                        y = group,
#                        fill = group,
#                        group= group)) +
#   ggtitle(title)+
#   ggridges::geom_density_ridges()+
#   labs(x=paste0("CK pair interaction score"), y="marker group",fill = "marker group") +
#   theme(panel.background = element_rect(fill = 'white', colour = 'white')) +
#   theme(axis.line = element_line(color="black", size = 0.3))+
#   theme(axis.text=element_text(size=12),
#         axis.title=element_text(size=14,face="bold")) +
#   theme(plot.title = element_text(hjust = 0.5, size=16,face="bold"))+
#   xlim(-5, 10)
# print(g)
# dev.off()
# 
# 
################################################################################################
# pleural vs Peritoneal
#
library(ggplot2)
source("utiles.R")
# group analysis, ridge density plot
ptype1 = "CD8"
ptype2 = "CK"
tissue = "Tumor" #NA # 
map = as.data.frame(fread("../Data_folder/TMA1/mapping.csv",check.names = FALSE))
cores_pleural_ori = map[map$CaseType == "Pleural", "Annotation ID"]
cores_peritoneal_ori = map[map$CaseType == "Peritoneal", "Annotation ID"]

if (is.na(tissue)){
  df_score_ori = read.csv( paste0("./output/cellContact/",ptype1, "_", ptype2, "_scores_wCoreInfo.csv"), check.names = FALSE)
}else{
  df_score_ori = read.csv( paste0("./output/cellContact/",ptype1, "_", ptype2, "_", tissue ,"_scores_wCoreInfo.csv"), check.names = FALSE)
}

# df_score = df_score[(df_score[paste0(ptype2, ' percent')] > 0.05) & (df_score$'tumor density' > 400), ]
# df_score = df_score_ori[(df_score_ori[paste0(ptype1, ' percent')] > 0.01) & (df_score_ori$'tumor percent' > 0.3), ]
df_score = df_score_ori[(df_score_ori[paste0(ptype1, ' percent')] > 0.01) & (df_score_ori$'CK percent' > 0.2 ) ,]


cores_pleural = intersect(cores_pleural_ori, df_score$'core ID')
cores_peritoneal = intersect(cores_peritoneal_ori, df_score$'core ID')

df_score[df_score$'core ID' %in% cores_pleural,'group'] = "Pleural"
df_score[df_score$'core ID' %in% cores_peritoneal,'group'] = "Peritoneal"

df_plot = na.omit(df_score[,c('pair_score','group')])


pvalue = wilcox.test(df_plot[df_plot$group == "Pleural", 'pair_score'], df_plot[df_plot$group == "Peritoneal", 'pair_score'])$p.value
pvalue = formatC(pvalue, format = "e", digits = 2)


title = paste0(ptype1, " to ", ptype2, " contacts")
# g = ggplot(df_plot, aes(x = pair_score,
#                         y = group,
#                         fill = group,
#                         group= group)) +
#   ggtitle(title)+
#   ggridges::geom_density_ridges()+
#   labs(x=paste0("Contact score,  pvalue = ", pvalue), y="group",fill = "group") +
#   theme(panel.background = element_rect(fill = 'white', colour = 'white')) +
#   theme(axis.line = element_line(color="black", size = 0.3))+
#   theme(axis.text=element_text(size=12),
#         axis.title=element_text(size=14,face="bold")) +
#   theme(plot.title = element_text(hjust = 0.5, size=16,face="bold"))+
#   xlim(-6,8)
x = 'pair_score'
y = 'group'
ridge_density_plot(df_plot, x, y, title, pvalue)


#################################################
# subtype density plot
df_score = df_score_ori[(df_score_ori[paste0(ptype1, ' percent')] > 0.01) & (df_score_ori$'CK percent' > 0.2 ) ,]
df_score = merge(df_score_ori, map, by="core ID",by.x = "core ID", 
                  by.y = "Annotation ID", all.x = TRUE, all.y = FALSE)

df_score_sub = df_score[(df_score[paste0(ptype1, ' percent')] > 0.01) & (df_score$'CK percent' > 0.2 ) ,]

subtype = "Peritoneal" #"Pleural" #"Peritoneal"  # #Peritoneal" #
df_score_subtype = df_score_sub[df_score_sub$casetype == subtype, ]

df_score_subtype = df_score_subtype[df_score_subtype$subtype %in% c("epithelial or epithelioid", "biphasic"), ]


pvalue = wilcox.test(df_score_subtype[df_score_subtype$subtype=="epithelial or epithelioid", 'pair_score'], df_score_subtype[ df_score_subtype$subtype=='biphasic', 'pair_score' ])$p.value
pvalue = formatC(pvalue, format = "e", digits = 2)




title = paste0(ptype1, " to ", ptype2, " contacts for " , subtype, " subtypes")
ggplot(df_score_subtype, aes(x = pair_score,
                        y = subtype,
                        fill = subtype,
                        group= subtype)) +
  ggtitle(title)+
  ggridges::geom_density_ridges()+
  labs(x=paste0("contact score, pvalue = ", pvalue), y="group",fill = "group") +
  theme(panel.background = element_rect(fill = 'white', colour = 'white')) +
  theme(axis.line = element_line(color="black", size = 0.3))+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold")) +
  theme(plot.title = element_text(hjust = 0.5, size=16,face="bold"))+
  xlim(-1,8)

##################################################
# Panel 2 markerss
library(ggplot2)
# group analysis, ridge density plot
ptype1 = "CD8"
ptype2 = "CK"

tissue = "Tumor"
if (is.na(tissue)){
  df_score_ori = read.csv( paste0("./output/cellContact/",ptype1, "_", ptype2, "_scores_wCoreInfo.csv"), check.names = FALSE)
}else{
  df_score_ori = read.csv( paste0("./output/cellContact/",ptype1, "_", ptype2, "_", tissue ,"_scores_wCoreInfo.csv"), check.names = FALSE)
}

# df_score = df_score[(df_score[paste0(ptype2, ' percent')] > 0.05) & (df_score$'tumor density' > 400), ]
# df_score = df_score_ori[(df_score_ori[paste0(ptype1, ' percent')] > 0.05) & (df_score_ori$'tumor density' > quantile(df_score_ori$'tumor density', 0.25))& (df_score_ori$'CK density' > quantile(df_score_ori$'CK density', 0.25)), ]
# df_score = df_score_ori[(df_score_ori[paste0(ptype1, ' percent')] > 0.05) & (df_score_ori$'tumor percent' > 0.25) & (df_score_ori$'CK percent' > 0.25), ]

df_score = df_score_ori[(df_score_ori[paste0(ptype1, ' percent')] > 0.05) & (df_score_ori$'CK percent' > 0.2), ]
casetype = "Pleural" #Peritoneal"
df_score = df_score[df_score$casetype == casetype, ]

marker = "BAP1"
groupname = paste0(marker," density")
df_plot = df_score[,c(groupname, 'pair_score')]
# df_plot = df_plot[df_plot$'pair_score' > -6, ]
df_plot = na.omit(df_plot)
# markerSep = c(200,800)
q50 = quantile(df_plot[,groupname], 0.5)
markerSep = c(q50, q50)
highCond = df_plot[groupname] >= markerSep[2]
df_plot[highCond, 'status'] = paste0(marker, " high")
lowCond = df_plot[groupname] < markerSep[1]
df_plot[lowCond, 'status'] = paste0(marker, " low")
df_plot$status = factor(df_plot$status, levels = c('BAP1 high', 'BAP1 low'))

pvalue = wilcox.test(df_plot[highCond, 'pair_score'], df_plot[lowCond, 'pair_score'])$p.value
pvalue = formatC(pvalue, format = "e", digits = 2)



df_plot = na.omit(df_plot)

title = paste0(ptype1, " to ", ptype2, " contacts with ",marker," high and low")
if (! is.na(casetype)){
  title = paste0(casetype, " ", title) 
}
ggplot(df_plot, aes(x = pair_score,
                        y = status,
                        fill = status,
                        group= status)) +
  ggtitle(title)+
  ggridges::geom_density_ridges()+
  labs(x=paste0("CD8~CK contact score"),fill = "group") +
  theme(panel.background = element_rect(fill = 'white', colour = 'white')) +
  theme(axis.line = element_line(color="black", size = 0.3))+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  theme(plot.title = element_text(hjust = 0.5, size=16))+
  guides(fill="none")+
  # theme(legend.text=element_text(size=8),
  #       legend.title=element_blank()) +
  xlim(-3,7) +
  ylab('') + 
  annotate("text", x=-1, y=3, label= paste0("pvalue = ", pvalue), size = 5)
    
ggsave(fig, width = 7, height = 3, dpi = 500, file=paste0("/Paper/Figure/CellContact/CD8~CK contact score with BAP1 high vs low.pdf"))


#
#
#
# groupname = "NF2 density"
# highCond = df_score[groupname] >= markerSep[2]
# df_score[highCond, 'NF2 status'] = "NF2+"
# lowCond = df_score[groupname] < markerSep[1]
# df_score[lowCond, 'NF2 status'] = "NF2-"
#
# groupname = "MTAP density"
# highCond = df_score[groupname] >= markerSep[2]
# df_score[highCond, 'MTAP status'] = "MTAP+"
# lowCond = df_score[groupname] < markerSep[1]
# df_score[lowCond, 'MTAP status'] = "MTAP-"
#
# df_score[is.na(df_score)]=""
# df_score['group'] = paste0(df_score$`BAP1 status`, df_score$`NF2 status`, df_score$`MTAP status`)
#
# df_score['group'] = df_score$`BAP1 status`
#
# groups = table(df_score['group'])
#
# library(ggplot2)
# # ridge density plot
# path_figure = "./output/cellContact/figure"
# if (is.na(tissue) & is.na(casetype) ){
#   title = paste0("CK~", ptype2, " interactions among marker groups in all region")
# }else if (is.na(casetype)){
#   title = paste0("CK~", ptype2, " interactions among marker groups in ",tolower(tissue)," region")
# }else if (is.na(tissue)){
#   title = paste0(casetype, " CK~", ptype2, " interactions among marker groups")
# }else{
#   title = paste0(casetype, " CK~", ptype2, " interactions among marker groups in ",tolower(tissue)," region")
# }
# png(file= file.path(path_figure, paste0("Ridge ",title,".png")),width=600, height=200)
# g = ggplot(df_score, aes(x = pair_score,
#                        y = group,
#                        fill = group,
#                        group= group)) +
#   ggtitle(title)+
#   ggridges::geom_density_ridges()+
#   labs(x=paste0("CK pair interaction score"), y="marker group",fill = "marker group") +
#   theme(panel.background = element_rect(fill = 'white', colour = 'white')) +
#   theme(axis.line = element_line(color="black", size = 0.3))+
#   theme(axis.text=element_text(size=12),
#         axis.title=element_text(size=14,face="bold")) +
#   theme(plot.title = element_text(hjust = 0.5, size=16,face="bold"))+
#   xlim(-5, 10)
# print(g)
# dev.off()
#
#
# ################################################################################################
#
#
