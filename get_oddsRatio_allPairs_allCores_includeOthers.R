# this version is output all contact scores for all cores for all region/ tumor region
library(tripack)
library(stats)
library(data.table)

library(foreach)
library(doParallel)
library(tidyr)

source("utiles.R")
rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# args = commandArgs(trailingOnly=TRUE)
# casetype = args[1]

# setwd(system("pwd", intern = T) )

getLogsOddsRatioMtx <- function(x, y, celltypes){
  
  
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

  for (type1 in typesetsub){
    for (type2 in typesetsub){
      
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


getLogsOddsRatioMtx_parrallel <- function(x, y, celltypes){
  
  
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
  
  #########################################
  
  fun <- function(i){
    type1 = typesetsub[i]
    
    oddsRatio_list = c()
    for (type2 in typesetsub){
      
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
      # bd = dbeta(seq01, InteractionsBtw, (TotalInteractionsOfThisType - InteractionsBtw))
      bd = dbeta(seq01, InteractionsBtw , (TotalInteractionsOfThisType + TotalInteractionsOfOtherType - InteractionsBtw))
      oddsRatio = mean(bd) / ((countCellsOtherType / length(nbs)) * (countCellsOfThisType / length(nbs) ))
      # oddsRatio = mean(bd) / ((TotalInteractionsOfThisType / TotalInteractions) * (TotalInteractionsOfOtherType / TotalInteractions ))
      
      # seq01 = seq(0.01, 0.99, by = 0.01)
      # bd = dbeta(seq01, InteractionsBtw , (TotalInteractionsOfThisType - InteractionsBtw))
      # # bd = dbeta(seq01, InteractionsBtw , (TotalInteractionsOfThisType + TotalInteractionsOfOtherType - InteractionsBtw))
      # # oddsRatio = mean(bd) / ((countCellsOtherType / length(nbs)) * (countCellsOfThisType / length(nbs) ))
      # # oddsRatio = mean(bd) / ((TotalInteractionsOfThisType / TotalInteractions) * (TotalInteractionsOfOtherType / TotalInteractions ))
      # oddsRatio = mean(bd) / (countCellsOtherType / length(nbs))
      
      oddsRatio_list <- c(oddsRatio_list, log(oddsRatio+.Machine$double.eps))
    }
    return(oddsRatio_list)
  }
  
  
  # CPUcores=detectCores()
  # cl <- makeCluster(CPUcores[1]-1) #not to overload your computer
  # registerDoParallel(cl)
  
  correl <- foreach (i = 1: length(typesetsub),  .combine=rbind,  .inorder =TRUE) %dopar% {
    rtn <- fun(i)
    rtn
  }
  correl <- as.data.frame(correl)
  rownames(correl) = typesetsub
  colnames(correl) = typesetsub
  
  
  return(correl)
}
########################################################################################

setwd(system("pwd", intern = T) )

# args = commandArgs(trailingOnly=TRUE)

args = commandArgs(trailingOnly=TRUE)

casetype = NA
# casetype = NA #"Peritoneal" #Pleural" # NA #args[1]
tissue = NA #args[1] #NA #"Tumor" #NA # 'Stroma' #NA # args[3]
print(tissue)
  
# print(casetype)
# print(tissue)


df1 = fread("../Data_folder/TMA1/panel1_data_cleaned_malignant_final.csv")
df1 = as.data.frame(df1)


if ( is.na(tissue)){
  df <- df1
}else{
  df <- df1[df1$'Tissue Category'== tissue , ]
}
cores = unique(df$"Annotation ID")

typeset <-  c('CD20','CD4','CD68','CD8','FOXP3','CK')
# typeset <-  c('CD4','CD68','CD8','CK')

cellContact_list = list()

CPUcores<-detectCores()
cl <- makeCluster(CPUcores[1]-1) #not to overload your computer
registerDoParallel(cl)

for (core in cores[1:10]){

    if (is.na(tissue) ) {
      dfsub = df[(df$'Annotation ID'==core), c("Cell X Position", "Cell Y Position", "phenotype_combined")]
    }else{
      dfsub = df[(df$'Annotation ID'==core) & (df$'Tissue Category'== tissue) ,c("Cell X Position", "Cell Y Position", "phenotype_combined")]
    }
  
    if (dim(dfsub)[1]<5) {next}
  
    x = dfsub$"Cell X Position"
    y = dfsub$"Cell Y Position"
    celltypes = dfsub$phenotype_combined

    oddsRatios <- getLogsOddsRatioMtx_parrallel(x,y, celltypes)
    oddsRatios <- getLogsOddsRatioMtx(x,y, celltypes)
    
    cellContact <- data.frame(matrix(NA, ncol = length(typeset), nrow = length(typeset)))
    rownames(cellContact) <- typeset
    colnames(cellContact) <- typeset
    
    # oddsRatios[oddsRatios < -36]= NA
    
    commSet = intersect(typeset, rownames(oddsRatios) )
    cellContact[commSet, commSet] <- oddsRatios  
    
    cellContact_list[[core]] <- cellContact


}


# results = as.data.frame(apply(simplify2array(lapply(cellContact_list, as.matrix)),1:2,mean, na.rm = TRUE))  

# get the element wise mean ignore NA
if (is.na(tissue)){
  filename = "./output/cellContact/cellContact_allscores_allRegion.RDS"
}else{
  filename = "./output/cellContact/cellContact_allscores_tumorRegion.RDS"
}

saveRDS(cellContact_list, file=filename)


stopCluster(cl)

######################################################################
####################################################################
# read in the outputs and analyses
regions = c('allRegion', 'tumorRegion')
for (region in regions){
  cellContact_list = readRDS(paste0("./output/cellContact/cellContact_allscores_", region,".RDS"))
  # caclulate the median of contact score for each pair
  mtx_median = as.data.frame(apply(simplify2array(lapply(cellContact_list, as.matrix)),1:2,median, na.rm = TRUE))  
  write.csv(mtx_median, file = paste0("./output/cellContact/cellContact_score_median_allCores_", region,".csv"))
  
  # convert list to array [6,6,315]
  cellContact_arr = array(unlist(cellContact_list),dim=c(6,6,length(cellContact_list)))
  # reorder the array
  cellContact_arr_reOrdered = aperm( cellContact_arr, c(3,1,2))  #315,6,6 (row, column, length)
  # for example, plot CK-CD20 with
  # hist(cellContact_arr_reOrdered[,6,1])
  saveRDS(cellContact_arr_reOrdered, file=paste0("./output/cellContact/cellContact_arr_reOrderedt_", region,".RDS"))
  
  # extract pleural,and peritoneal cores
  df_feature = read.csv("./output/cellContact/core_features_allMarkers_withIntensity_core.csv", check.names=FALSE)
  pleural_cores = df_feature$"core ID"[df_feature$casetype == "Pleural"]
  peritoneal_cores = df_feature$"core ID"[df_feature$casetype == "Peritoneal"]
  
  
  
  all_cores = names(cellContact_list)
  pleural_cores = intersect(pleural_cores, all_cores)
  peritoneal_cores = intersect(peritoneal_cores, all_cores)
  
  # cellContact for pleural and peritoneal
  cellContact_list_pleural = cellContact_list[pleural_cores] #232
  cellContact_list_peritoneal = cellContact_list[peritoneal_cores] #67
  
  # cellContact median for pleural and peritoneal
  mtx_median_pleural = as.data.frame(apply(simplify2array(lapply(cellContact_list_pleural, as.matrix)),1:2,median, na.rm = TRUE))  
  write.csv(mtx_median_pleural, file =  paste0("./output/cellContact/cellContact_score_median_pleural", region,".csv"))
  
  
  mtx_median_peritoneal = as.data.frame(apply(simplify2array(lapply(cellContact_list_peritoneal, as.matrix)),1:2,median, na.rm = TRUE))  
  write.csv(mtx_median_peritoneal, file =  paste0("./output/cellContact/cellContact_score_median_peritoneal", region,".csv"))

  #save the reordered cellContact_arr_reOrdered for Pleual and Peritoneal
  pleural_cores_indices = which(names(cellContact_list) %in% pleural_cores)
  peritoneal_cores_indices = which(names(cellContact_list) %in% peritoneal_cores)
  
  cellContact_arr_reOrdered_pleural = cellContact_arr_reOrdered[ pleural_cores_indices,,] #232   6   6
  cellContact_arr_reOrdered_peritoneal = cellContact_arr_reOrdered[ peritoneal_cores_indices,,] # 67  6  6
  
  saveRDS(cellContact_arr_reOrdered_pleural, file=paste0("./output/cellContact/cellContact_array_pleural_", region,".RDS"))
  saveRDS(cellContact_arr_reOrdered_peritoneal, file=paste0("./output/cellContact/cellContact_array_peritoneal_", region,".RDS"))
  
} 

######################################################

# Panel2 marker high vs low
celltypes = c('CD20'=1, 'CD4'=2, 'CD68'=3, 'CD8'=4, 'FOXP3'=5)
region='tumorRegion'

P2marker = "NF2"
measure = paste0(P2marker, " density")
df_feature_ori = read.csv("./output/cellContact/core_features_allMarkers_withIntensity_core.csv", check.names=FALSE)
contact_list = readRDS(paste0("./output/cellContact/cellContact_allscores_", region,".RDS"))



pvalues = c()
num_highs = c()
num_lows = c()
for (study in c("all", "Pleural", "Peritoneal")){
  if (study == "all"){
    df_feature_study = df_feature_ori
  }else{
    df_feature_study = df_feature_ori[df_feature$casetype == study,]
  }
  pvalues = c(pvalues, study)
  num_highs = c(num_highs, NA)
  num_lows = c(num_lows, NA)
  for (i in seq(1:5)){ #CD20        CD4       CD68        CD8       FOXP3  
   
    celltype = names(celltypes[i])
    
    CK_cond =  df_feature_study$`CK percent` > 0.005
    
    BAP1_cond = df_feature_study[measure] > 0
    type_cond = df_feature_study[paste0(celltype, " percent")] > 0
   
  
    df_feature = df_feature_study[CK_cond & type_cond & BAP1_cond, ]
    
    # df_feature = df_feature_study
    
    cond_high = df_feature[,measure] > quantile(df_feature[,measure], 0.6)
    cond_low = df_feature[,measure] <= quantile(df_feature[,measure], 0.4)
    
    high_cores = df_feature$`core ID`[cond_high]
    low_cores = df_feature$`core ID`[cond_low]
    
    
    contact_list_high = contact_list[high_cores ]
    contact_list_low = contact_list[low_cores]
    
    contact_high_arr = array(unlist(contact_list_high),dim=c(6,6,length(contact_list_high)))
    contact_high_arr_reOrdered = aperm(contact_high_arr, c(3,1,2))  #152,6,6 (row, column, length)
    
    contact_low_arr = array(unlist(contact_list_low),dim=c(6,6,length(contact_list_low)))
    contact_low_arr_reOrdered = aperm(contact_low_arr, c(3,1,2))  #315,6,6 (row, column, length)
    
    x = contact_high_arr_reOrdered[,6,i]
    y = contact_low_arr_reOrdered[,6,i]
    
    # x[x<-36] = NA
    # y[y<-36] = NA
    x = x[! is.na(x)]
    y = y[! is.na(y)]
    
    if (length(x)==0 | length(y)==0){
      pvalue=NA
    }else{
      pvalue = wilcox.test(x,y)$p.value
      pvalue = formatC(pvalue, format = "e", digits = 2)
    }
    
    num_high = length(x)
    num_low = length(y)
    
    pvalues = c(pvalues, pvalue)
    num_highs = c(num_highs, num_high)
    num_lows = c(num_lows, num_low)
    # print(paste0("CK ~ ", celltype))
    # print(wilcox.test(x,y)$p.value)
    # print(paste0("high=",length(x), ", low=",length(y)))
    # print("  ")
  }  
    
  }
  
df_pvalues = as.data.frame(cbind(pvalues, num_highs, num_lows)) 
df_pvalues$types = rep(c(NA, "CD20","CD4", "CD68","CD8","FOXP3"), 3)
df_pvalues = df_pvalues[c("types", "pvalues", "num_highs", "num_lows")]

write.csv(df_pvalues, paste0("./output/cellContact/",P2marker, "_high_vs_low_pvalues_", region,".csv"))  
 
############################################################## 
# # save the score for density plot for CD8

celltypes = c('CD20'=1, 'CD4'=2, 'CD68'=3, 'CD8'=4, 'FOXP3'=5)
region='tumorRegion'

P2marker = "BAP1"
measure = paste0(P2marker, " density")
df_feature_ori = read.csv("./output/cellContact/core_features_allMarkers_withIntensity_core.csv", check.names=FALSE)
contact_list = readRDS(paste0("./output/cellContact/cellContact_allscores_", region,".RDS"))





for (study in c("all", "Pleural", "Peritoneal")){
  if (study == "all"){
    df_feature_study = df_feature_ori
  }else{
    df_feature_study = df_feature_ori[df_feature_ori$casetype == study,]
  }

  i = 4  #CD8  (#CD20        CD4       CD68        CD8       FOXP3  )
    
  celltype = names(celltypes[i])
  
  CK_cond =  df_feature_study$`CK percent` > 0.2
  
  BAP1_cond = df_feature_study[measure] > 0
  type_cond = df_feature_study[paste0(celltype, " percent")] > 0.05
  
  
  df_feature = df_feature_study[CK_cond & type_cond & BAP1_cond, ]
  
  # df_feature = df_feature_study
  
  cond_high = df_feature[,measure] > quantile(df_feature[,measure], 0.6)
  cond_low = df_feature[,measure] <= quantile(df_feature[,measure], 0.4)
  
  high_cores = df_feature$`core ID`[cond_high]
  low_cores = df_feature$`core ID`[cond_low]
 
  contact_list_high = contact_list[high_cores ]
  contact_list_low = contact_list[low_cores]
  
  contact_high_arr = array(unlist(contact_list_high),dim=c(6,6,length(contact_list_high)))
  contact_high_arr_reOrdered = aperm(contact_high_arr, c(3,1,2))  #152,6,6 (row, column, length)
  
  contact_low_arr = array(unlist(contact_list_low),dim=c(6,6,length(contact_list_low)))
  contact_low_arr_reOrdered = aperm(contact_low_arr, c(3,1,2))  #315,6,6 (row, column, length)
  
  score_h = contact_high_arr_reOrdered[,6,i]
  score_l = contact_low_arr_reOrdered[,6,i]
    
    # x[x<-36] = NA
    # y[y<-36] = NA
  score_h = score_h[! is.na(score_h)]
  score_l = score_l[! is.na(score_l)]
  
  if (length(score_h)==0 | length(score_l)==0){
    pvalue=NA
  }else{
    pvalue = wilcox.test(score_h,score_l)$p.value
    pvalue = formatC(pvalue, format = "e", digits = 2)
  }
  
  num_high = length(score_h)
  num_low = length(score_l)
  

  
  pair_score = c(score_h, score_l)
  group = c(rep("BAP1+", length(score_h)), rep("BAP1-", length(score_l)))
  
  df_plot = as.data.frame(cbind(pair_score, group))
  df_plot$pair_score = as.numeric(df_plot$pair_score)
  title = "CK~CD8 contacts BAP1 high vs low"
  ridge_density_plot(df_plot, 'pair_score', 'group', title, pvalue)
  
}  



ggplot(df_plot, aes(x = pair_score,
                                              y = group,
                                              fill = group,
                                              group= group)) +
                        ggtitle(title)+
                        ggridges::geom_density_ridges()+
                        labs(x=paste0("Contact score,  pvalue = ", pvalue), y="group",fill = "group") +
                        theme(panel.background = element_rect(fill = 'white', colour = 'white')) +
                        theme(axis.line = element_line(color="black", size = 0.3))+
                        theme(axis.text=element_text(size=12),
                              axis.title=element_text(size=14,face="bold")) +
                        theme(plot.title = element_text(hjust = 0.5, size=16,face="bold"))+
                        xlim(-6,8)



# 
# 
# 
# 
# 
# ridge_plot = function(){
#   pvalue = formatC(pvalue, format = "e", digits = 2)
#   
#   
#   title = paste0(ptype1, " to ", ptype2, " contacts")
#   g = ggplot(df_plot, aes(x = pair_score,
#                           y = group,
#                           fill = group,
#                           group= group)) +
#     ggtitle(title)+
#     ggridges::geom_density_ridges()+
#     labs(x=paste0("Contact score,  pvalue = ", pvalue), y="group",fill = "group") +
#     theme(panel.background = element_rect(fill = 'white', colour = 'white')) +
#     theme(axis.line = element_line(color="black", size = 0.3))+
#     theme(axis.text=element_text(size=12),
#           axis.title=element_text(size=14,face="bold")) +
#     theme(plot.title = element_text(hjust = 0.5, size=16,face="bold"))+
#     xlim(-6,8)
#   print(g)
#   
# }
# 
# 
# 
# 
# 
# 
# 
# # CD20        CD4       CD68        CD8       FOXP3         CK
# # CD20  -36.04365 -36.043653 -36.043653 -36.043653 -11.2160513 -36.043653
# # CD4   -36.04365   2.026591   5.027765   4.417843   2.9307998   2.328015
# # CD68  -36.04365   5.027765   4.580675   4.893693   2.7254391   3.053849
# # CD8   -36.04365   4.417843   4.893693   3.768922   2.6616655   2.517379
# # FOXP3 -11.21605   2.930800   2.725439   2.661666   0.6103592   1.872857
# # CK    -36.04365   2.328015   3.053849   2.517379   1.8728574   0.419957