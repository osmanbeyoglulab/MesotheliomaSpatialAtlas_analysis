
library(stats)
library(data.table)
library(tidyr)
library("dplyr")

rm(list=ls())

# args = commandArgs(trailingOnly=TRUE)
# casetype = args[1]


#########################################################################################


# setwd(system("pwd", intern = T) )
dir = dirname(rstudioapi::getSourceEditorContext()$path); print(dir)  # [1] "/Users/lees18/OneDrive - UPMC/S86_Covid19_scOmics_Xiaosjun/3_SeuratObject_ClusteringByPatient"
setwd(dir)

args = commandArgs(trailingOnly=TRUE)
map = fread("../Data_folder/TMA1/mapping.csv")
map = as.data.frame(map) # 341 cores, 115 MVB
df1 = fread("../Data_folder/TMA1/panel1_data_cleaned_withCelltypes_03-29.csv")
df1 = as.data.frame(df1)
df1 = merge(df1, map,by="Annotation ID") #341 cores, 115 MVB
df1$phenotype_combined[df1$phenotype_combined=="other"] = "Unidentified"


df2 = fread("../Data_folder/TMA2/panel2_data_cleaned_withCelltypes_03-29.csv")
df2 = as.data.frame(df2)
df2 = merge(df2, map,by="Annotation2") #326 cores #113 MVB. #dim 606114    193

# df2_area = fread("../Data_folder/TMA2/densities_tma2_2.csv"). #dim 444  10
# df2_area = as.data.frame(df2_area)
df2_area = as.data.frame(fread("../Data_folder/TMA2/densities_tma2.csv", check.names=FALSE))
df2 =  merge(df2, df2_area, by="Annotation2") #606114    202

# select malignant cores
cores_all = unique(df2[df2$'Classification'=='Malignant', "Annotation ID"])  # 22 non-malignant
# #304 cores
markers1 = c('FOXP3', 'CD4', 'CD8', 'CD68', 'CD20', 'CK', "Unidentified")
markers2 = c('CD11c', 'CD56', 'BAP1', 'NF2', 'MTAP', 'LAG3') 

intensityCol = list()
for (marker in markers1[1:6]){
  if (marker == "FOXP3")
    intensityCol[[marker]] = "Entire Cell FOXP3 (Opal 570) Mean (Normalized Counts, Total Weighting)"
  if (marker == "CD4")
    intensityCol[[marker]] = "Entire Cell CD4 (Opal 690) Mean (Normalized Counts, Total Weighting)"
  if (marker == "CD8")
    intensityCol[[marker]]  = "Entire Cell CD8 (Opal 480) Mean (Normalized Counts, Total Weighting)"
  if (marker == "CD20")
    intensityCol[[marker]]  =  "Entire Cell CD20 (Opal 620) Mean (Normalized Counts, Total Weighting)"
  if (marker == "CK")
    intensityCol[[marker]]  = "Entire Cell PanCK (Opal 780) Mean (Normalized Counts, Total Weighting)"
  if (marker == "CD68")
    intensityCol[[marker]]  = "Entire Cell CD68 (Opal 520) Mean (Normalized Counts, Total Weighting)"
}
for (marker in markers2){
  if (marker == "CD11c")
    intensityCol[[marker]]  = "Entire Cell CD11c (Opal 480) Mean (Normalized Counts, Total Weighting)"
  if (marker == "CD56")
    intensityCol[[marker]]  = "Entire Cell CD56 (Opal 520) Mean (Normalized Counts, Total Weighting)"
  if (marker == "NF2")
    intensityCol[[marker]] = "Entire Cell NF2 (Opal 780) Mean (Normalized Counts, Total Weighting)"
  if (marker == "BAP1")
    intensityCol[[marker]]  = "Entire Cell BAP1 (Opal 690) Mean (Normalized Counts, Total Weighting)"
  if (marker == "MTAP")
    intensityCol[[marker]] = "Entire Cell MTAP (Opal 620) Mean (Normalized Counts, Total Weighting)"
  if (marker == "LAG3")
    intensityCol[[marker]] = "Entire Cell LAG3 (Opal 570) Mean (Normalized Counts, Total Weighting)"
}


cores_sel = c()
percents_CD4 = c()
percents_CD20 = c()
percents_FOXP3 = c()
percents_CD4 = c()
percents_CD68 = c()
percents_CD8 = c()
percents_CK = c()
percents_tumor = c()
percents_BAP1 = c()
percents_NF2 = c()
percents_MTAP = c()
percents_LAG3 = c()
percents_CD11c = c()
percents_CD56 = c()
percents_Unidentified = c()

densities_CD4 = c()
densities_CD20 = c()
densities_FOXP3 = c()
densities_CD68 = c()
densities_CD8 = c()
densities_CK = c()
densities_tumor = c()
densities_BAP1 = c()
densities_NF2 = c()
densities_MTAP = c()
densities_LAG3 = c()
densities_CD11c = c()
densities_CD56 = c()
densities_Unidentified = c()
# 
# densities_CD4_tumor = c()
# densities_CD20_tumor = c()
# densities_FOXP3_tumor = c()
# densities_CD68_tumor = c()
# densities_CD8_tumor = c()
# densities_CK_tumor = c()
# 
# 
# densities_CD4_stroma = c()
# densities_CD20_stroma = c()
# densities_FOXP3_stroma = c()
# densities_CD68_stroma = c()
# densities_CD8_stroma = c()
# densities_CK_stroma = c()

intensities = list()
casetypes = c()
subtypes = c()
MVBs = c()

n1s_cells = c()
n2s_cells = c()
notumor_counts=0
for (c in cores_all){ 
  
  core1 = df1[df1$'Annotation ID' == c, ]
  n1_cells = dim(core1)[1]
  n_tumor = sum(core1$'Tissue Category' == 'Tumor')
 
  area1_tumor = core1$"tumorArea"[1]
  area1_stroma = core1$"stromaArea"[1]
  area1 = area1_tumor +  area1_stroma
  
  casetype = df1[df1$"Annotation ID"== c, 'CaseType'][1]
  MVB = df1[df1$"Annotation ID"== c, 'MVB'][1]
  subtype = df1[df1$"Annotation ID"== c, 'subtype'][1]
  
  core2 = df2[df2$"Annotation ID" == c, ]
  n2_cells = dim(core2)[1]
  area2 = core2[, "Tissue Area (mm2)"][1]
  if( is.na(area2) ) {area2=area1}
  
  percent_tumor = n_tumor/n1_cells
  density_tumor = n_tumor/area1
  
  n_BAP1 = sum(core2[sprintf("Phenotype-%s",'BAP1')]==sprintf("%s+",'BAP1'))
  percent_BAP1 = n_BAP1 / n2_cells
  density_BAP1 = n_BAP1 / area2
    
  n_NF2 = sum(core2[sprintf("Phenotype-%s",'NF2')]==sprintf("%s+",'NF2'))
  percent_NF2 = n_NF2 / n2_cells
  density_NF2 = n_NF2 / area2
  
  n_MTAP = sum(core2[sprintf("Phenotype-%s",'MTAP')]==sprintf("%s+",'MTAP'))
  percent_MTAP = n_MTAP / n2_cells
  density_MTAP = n_MTAP / area2
  
  n_LAG3 = sum(core2[sprintf("Phenotype-%s",'LAG3')]==sprintf("%s+",'LAG3'))
  percent_LAG3 = n_LAG3 / n2_cells
  density_LAG3 = n_LAG3 / area2
  
  n_CD11c = sum(core2[sprintf("Phenotype-%s",'CD11c')]==sprintf("%s+",'CD11C'))
  
  percent_CD11c = n_CD11c / n2_cells
  density_CD11c = n_CD11c / area2
  
  print(n_CD11c)
  print(percent_CD11c)
  
  
  n_CD56 = sum(core2[sprintf("Phenotype-%s",'CD56')]==sprintf("%s+",'CD56'))
  percent_CD56 = n_CD56 / n2_cells
  density_CD56 = n_CD56 / area2
  
  percents_BAP1 = c(percents_BAP1, percent_BAP1)
  percents_NF2 = c(percents_NF2, percent_NF2)
  percents_MTAP = c(percents_MTAP, percent_MTAP)
  percents_LAG3 = c(percents_LAG3, percent_LAG3)
  percents_CD11c = c(percents_CD11c, percent_CD11c)
  percents_CD56 = c(percents_CD56, percent_CD56)
 
  densities_BAP1 = c(densities_BAP1, density_BAP1)
  densities_NF2 = c(densities_NF2, density_NF2)
  densities_MTAP = c(densities_MTAP, density_MTAP)
  densities_LAG3 = c(densities_LAG3, density_LAG3)
  densities_CD11c = c(densities_CD11c, density_CD11c)
  densities_CD56 = c(densities_CD56, density_CD56)
  
  n_CD8 = sum(core1$phenotype_combined == 'CD8')
  n_CK = sum(core1$phenotype_combined == 'CK')
  n_CD4 = sum(core1$phenotype_combined == 'CD4')
  n_CD20 = sum(core1$phenotype_combined == 'CD20')
  n_CD68 = sum(core1$phenotype_combined == 'CD68')
  n_FOXP3 = sum(core1$phenotype_combined == 'FOXP3')
  n_Unidentified = sum(core1$phenotype_combined == 'Unidentified')
  
  
  percents_CD8 = c(percents_CD8, n_CD8/n1_cells)
  percents_CK = c(percents_CK, n_CK/n1_cells)
  percents_CD4 = c(percents_CD4, n_CD4/n1_cells)
  percents_CD68 = c(percents_CD68, n_CD68/n1_cells)
  percents_FOXP3 = c(percents_FOXP3, n_FOXP3/n1_cells)
  percents_CD20 = c(percents_CD20, n_CD20/n1_cells)
  percents_Unidentified = c(percents_Unidentified, n_Unidentified/n1_cells)
  
  
  densities_CD8 = c(densities_CD8, n_CD8/area1)
  densities_CK = c(densities_CK, n_CK/area1)
  densities_CD4 = c(densities_CD4, n_CD4/area1)
  densities_CD68 = c(densities_CD68, n_CD68/area1)
  densities_FOXP3 = c(densities_FOXP3, n_FOXP3/area1)
  densities_CD20 = c(densities_CD20, n_CD20/area1)
  densities_Unidentified = c(densities_Unidentified, n_Unidentified/area1)
  
  
  
  
  cores_sel = c(cores_sel, c)
  casetypes = c(casetypes, casetype)
  MVBs = c(MVBs, MVB)
  subtypes = c(subtypes, subtype)
  
  percents_tumor = c(percents_tumor, percent_tumor)    
  densities_tumor = c(densities_tumor, density_tumor)   
  
  
  n1s_cells = c(n1s_cells, n1_cells)
  n2s_cells = c(n2s_cells, n2_cells)
  
  for (marker in markers1[1:6]){
    intensity = mean(core1[,intensityCol[[marker]]])
    intensities[[marker]] = c(intensities[[marker]],intensity)
  }
  for (marker in markers2){
    intensity = mean(core2[,intensityCol[[marker]]])
    intensities[[marker]] = c(intensities[[marker]],intensity)
  }

}

df_out = as.data.frame(cbind(cores_sel, n1s_cells, n2s_cells, percents_tumor,densities_tumor,
                            percents_BAP1, percents_NF2, percents_MTAP, percents_LAG3,percents_CD11c, percents_CD56,
                            percents_CK, percents_CD8, percents_CD4, percents_CD20, percents_CD68, percents_FOXP3, percents_Unidentified,
                            densities_BAP1, densities_NF2, densities_MTAP, densities_LAG3,densities_CD11c, densities_CD56,
                            densities_CK, densities_CD8, densities_CD4, densities_CD20, densities_CD68, densities_FOXP3,densities_Unidentified,
                            intensities[[markers1[1]]],intensities[[markers1[2]]],intensities[[markers1[3]]],
                            intensities[[markers1[4]]],intensities[[markers1[5]]], intensities[[markers1[6]]],
                            intensities[[markers2[1]]],intensities[[markers2[2]]],intensities[[markers2[3]]],
                            intensities[[markers2[4]]],intensities[[markers2[5]]],intensities[[markers2[6]]],
                            casetypes, subtypes, MVBs))

colnames(df_out) = c('core ID', 'cell count 1', 'cell count 2','tumor percent','tumor_density',
                    'BAP1 percent','NF2 percent','MTAP percent', 'LAG3 percent',  'CD11c percent', 'CD56 percent',
                    'CK percent', 'CD8 percent', 'CD4 percent','CD20 percent', 'CD68 percent','FOXP3 percent', 'Unidentified percent',
                    'BAP1 density','NF2 density','MTAP density', 'LAG3 density', 'CD11c density', 'CD56 density',
                    'CK density', 'CD8 density', 'CD4 density','CD20 density', 'CD68 density','FOXP3 density', 'Unidentified density',
                    'FOXP3 intensity', 'CD4 intensity', 'CD8 intensity', 'CD68 intensity', 'CD20 intensity', 'CK intensity', 
                    'CD11c intensity', 'CD56 intensity','BAP1 intensity','NF2 intensity','MTAP intensity', 'LAG3 intensity',
                    'casetype','subtype', 'MVB')

df_out$subtype[df_out$subtype == "epithelial or epithelioid"] =  "epithelial"


write.csv(df_out, file = "./output/cellContact/core_features_allMarkers_withIntensity_core.csv", row.names=FALSE )


# print out the nubmer of core and MVB for each type and subtype
for (casetype in c("Pleural", "Peritoneal")){
  print(casetype)
  print(length(unique(df_out[(df_out$casetype==casetype) , "core ID"])))
  print(length(unique(df_out[(df_out$casetype==casetype), "MVB"])))
  for (subtype in c("epithelial", "biphasic", "sarcomatoid")){
    
    print(subtype)
    print(length(unique(df_out[(df_out$casetype==casetype) & (df_out$subtype==subtype), "core ID"])))
    print(length(unique(df_out[(df_out$casetype==casetype) & (df_out$subtype==subtype), "MVB"])))
  }
}



# generated 'Unidentified_combined percent' and 'Unidentified_combined density' cols, which are other-CD11c-CD56, and remove the cores 'Unidentified_combined percent'<0
df_feature_core = as.data.frame(fread("./output/cellContact/core_features_allMarkers_withIntensity_core.csv", check.names = FALSE))
df = df_feature_core
df$'Unidentified_combined percent' = df$'Unidentified percent' - df$'CD11c percent' - df$'CD56 percent'
df$'Unidentified_combined density' = df$'Unidentified density' - df$'CD11c density' - df$'CD56 density'

sum(df$'Unidentified_combined percent' <0) #58
sum(df$'Unidentified_combined density' <0) #46

df_out = df[df$'Unidentified_combined percent'>0, ]

# save the cores that CD56 and CD11c no more than unidentified, it is more accurate data
write.csv(df_out, file = "./output/cellContact/core_features_allMarkers_withIntensity_core_accurate.csv", row.names=FALSE )

# print out the nubmer of core and MVB for each type and subtype
for (casetype in c("Pleural", "Peritoneal")){
  print(casetype)
  print(length(unique(df_out[(df_out$casetype==casetype) , "core ID"])))
  print(length(unique(df_out[(df_out$casetype==casetype), "MVB"])))
  for (subtype in c("epithelial", "biphasic", "sarcomatoid")){
    
    print(subtype)
    print(length(unique(df_out[(df_out$casetype==casetype) & (df_out$subtype==subtype), "core ID"])))
    print(length(unique(df_out[(df_out$casetype==casetype) & (df_out$subtype==subtype), "MVB"])))
  }
}


#######################################################################################################################

df_feature = as.data.frame(fread("./output/cellContact/core_features_allMarkers_withIntensity_core.csv", check.names = FALSE))
df_feature_sub = df_feature[c( 'MVB','casetype','subtype','cell count 1', 'cell count 2','tumor percent','tumor_density',
                               'BAP1 percent','NF2 percent','MTAP percent', 'LAG3 percent',  'CD11c percent', 'CD56 percent',
                               'CK percent', 'CD8 percent', 'CD4 percent','CD20 percent', 'CD68 percent','FOXP3 percent', 'Unidentified percent',
                               'BAP1 density','NF2 density','MTAP density', 'LAG3 density', 'CD11c density', 'CD56 density',
                               'CK density', 'CD8 density', 'CD4 density','CD20 density', 'CD68 density','FOXP3 density','Unidentified density',
                               'FOXP3 intensity', 'CD4 intensity', 'CD8 intensity', 'CD68 intensity', 'CD20 intensity', 'CK intensity', 
                               'CD11c intensity', 'CD56 intensity','BAP1 intensity','NF2 intensity','MTAP intensity', 'LAG3 intensity'
                               
)]

df_feature_MVB = df_feature_sub %>% group_by(MVB,casetype,subtype) %>% summarise_all("median",na.rm = TRUE) %>% as.data.frame()

# MERGE DOES NOT WORK WHEN COLUMN VALUE CONTAINS SPACE
# df_feature_MVB2 = merge(x = df_feature_MVB, y = df_feature[c('MVB','casetype','subtype' )], by='MVB', all.x = TRUE)
write.csv(df_feature_MVB, file = "./output/cellContact/core_features_allMarkers_withIntensity_MVB_median.csv", row.names=FALSE )

df_feature_MVB =df_feature_sub %>% group_by(MVB,casetype,subtype) %>% summarise_all("mean", na.rm = TRUE) %>% as.data.frame()
# df_feature_MVB = merge(x = df_feature_MVB, y = df_feature[,c('MVB','casetype','subtype' )], by='MVB', all.x = TRUE)
write.csv(df_feature_MVB, file = "./output/cellContact/core_features_allMarkers_withIntensity_MVB_mean.csv", row.names=FALSE )

#--------------------------
df_feature_sub[df_feature_sub == 0] = NA

df_feature_MVB = df_feature_sub %>% group_by(MVB,casetype,subtype) %>% summarise_all("median",na.rm = TRUE) %>% as.data.frame()
write.csv(df_feature_MVB, file = "./output/cellContact/core_features_allMarkers_withIntensity_MVB_median_no0cores.csv", row.names=FALSE )

df_feature_MVB =df_feature_sub %>% group_by(MVB,casetype,subtype) %>% summarise_all("mean", na.rm = TRUE) %>% as.data.frame()
write.csv(df_feature_MVB, file = "./output/cellContact/core_features_allMarkers_withIntensity_MVB_mean_no0cores.csv", row.names=FALSE )




#######################################################################################################################

df_feature = as.data.frame(fread("./output/cellContact/core_features_allMarkers_withIntensity_core_accurate.csv", check.names = FALSE))
df_feature_sub = df_feature[c( 'MVB','casetype','subtype','cell count 1', 'cell count 2','tumor percent','tumor_density',
                            'BAP1 percent','NF2 percent','MTAP percent', 'LAG3 percent',  'CD11c percent', 'CD56 percent',
                            'CK percent', 'CD8 percent', 'CD4 percent','CD20 percent', 'CD68 percent','FOXP3 percent', 'Unidentified percent',
                            'BAP1 density','NF2 density','MTAP density', 'LAG3 density', 'CD11c density', 'CD56 density',
                            'CK density', 'CD8 density', 'CD4 density','CD20 density', 'CD68 density','FOXP3 density','Unidentified density',
                            'FOXP3 intensity', 'CD4 intensity', 'CD8 intensity', 'CD68 intensity', 'CD20 intensity', 'CK intensity', 
                            'CD11c intensity', 'CD56 intensity','BAP1 intensity','NF2 intensity','MTAP intensity', 'LAG3 intensity',
                            'Unidentified_combined percent', 'Unidentified_combined density'
                           )]

df_feature_MVB = df_feature_sub %>% group_by(MVB,casetype,subtype) %>% summarise_all("median",na.rm = TRUE) %>% as.data.frame()

# MERGE DOES NOT WORK WHEN COLUMN VALUE CONTAINS SPACE
# df_feature_MVB2 = merge(x = df_feature_MVB, y = df_feature[c('MVB','casetype','subtype' )], by='MVB', all.x = TRUE)
write.csv(df_feature_MVB, file = "./output/cellContact/core_features_allMarkers_withIntensity_MVB_median_accurate.csv", row.names=FALSE )

df_feature_MVB =df_feature_sub %>% group_by(MVB,casetype,subtype) %>% summarise_all("mean", na.rm = TRUE) %>% as.data.frame()
# df_feature_MVB = merge(x = df_feature_MVB, y = df_feature[,c('MVB','casetype','subtype' )], by='MVB', all.x = TRUE)
write.csv(df_feature_MVB, file = "./output/cellContact/core_features_allMarkers_withIntensity_MVB_mean_accurate.csv", row.names=FALSE )





##################################################################
# clustering
library(textshape)
library(pheatmap)
library(grid) #textGrob used to rotate row lable

df_features = as.data.frame(fread("./output/cellContact/core_features_allMarkers.csv", check.names = TRUE))
df_features = textshape::column_to_rownames(df_features)


df_map = as.data.frame(fread("../Data_folder/TMA1/mapping.csv"))
df_map = textshape::column_to_rownames((df_map))
df_map = df_map[rownames(df_features),]

df_features['MVB'] = df_map['MVB']

MVB_table = table(df_features$MVB)

MVB_sub = names(MVB_table[MVB_table>2])

df_features_sub = df_features[(df_features$MVB %in% MVB_sub),]


markers1 = c('FOXP3', 'CD4', 'CD8', 'CD68', 'CD20', 'CK', 'CD11c', 'CD56')
# markers2 = c('CD11c', 'CD56', 'BAP1', 'NF2', 'MTAP') #, 'LAG3'

# cols_sel = paste0(markers1, " density")
# df_plt = log(as.matrix(df_features_sub[,cols_sel]) + .Machine$double.eps)
# 
# MVB = factor(df_features_sub[,'MVB'])
# 
# annotation_row = data.frame(
#   "MVB" = MVB, check.names = FALSE)
# rownames(annotation_row) = rownames(df_plt)
# # group_colors <- as.character(palette_mapper[MVB])
# # names(group_colors) <- unique(MVB)
# 
# # ann_colors = list(
# #   "MVB" = group_colors
# # )
# 
# df_plt[is.na(df_plt)] <- 0
# 
# save_pheatmap_png <- function(x, filename, width=1000, height=2000) {
#   stopifnot(!missing(x))
#   stopifnot(!missing(filename))
#   png(filename, width=width, height=height,res=1200)
#   grid::grid.newpage()
#   grid::grid.draw(x$gtable)
#   dev.off()
# }
# 
# draw_colnames_45 <- function (coln, gaps, ...) {
#   coord = pheatmap:::find_coordinates(length(coln), gaps)
#   x = coord$coord - 0.5 + coord$size
#   res = textGrob(coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"), vjust = 0.5, hjust = 1, rot = 45, gp = gpar(...))
#   return(res)}
# 
# ## 'Overwrite' default draw_colnames with your own version 
# assignInNamespace(x="draw_colnames", value="draw_colnames_45",
#                   ns=asNamespace("pheatmap"))
# 
# xx = pheatmap(df_plt, annotation_row  = annotation_row, show_colnames=TRUE, show_rownames=FALSE) #, annotation_colors = ann_colors, main=substr(filename,1,nchar(filename)-4 ))
# save_pheatmap_png(xx, file.path("../figure/test", "heatmap.png"), width=10000, height=ncol(df)+170+1500)
# 

cols_sel = paste0(markers1, ".density")
df_plt = log(df_features_sub[,cols_sel]) + .Machine$double.eps
df_plt[,"MVB"] = df_features_sub[,'MVB']
library(dplyr)
df_plt2 <- df_plt %>% 
  filter_all(all_vars(!is.infinite(.)))

library(nlme)
df_plt2[,'res'] = 1

fm = lme(fixed = rep(1, nrow(df_plt2)) ~ FOXP3.density+CD4.density+CD8.density+CD68.density+CD20.density+CD11c.density+CD56.density  , data = df_plt2, random = ~ 1 | MVB)
fm = lme(fixed = CK.density+FOXP3.density+CD4.density+CD8.density+CD68.density ~ 1, data = df_plt2, random = ~ 1 |MVB )

##############################################################################
# calculate pvalues for casetype(/subtype) difference

inte_type = "mean"
df_feature_ori = as.data.frame(fread(paste0("./output/cellContact/core_features_allMarkers_withIntensity_MVB_", inte_type,".csv"), check.names = FALSE))

df_feature_ori[(df_feature_ori$'subtype'=='epithelial') | (df_feature_ori$'subtype'=='epithelial or epithelioid'), 'subtype' ] = "epithelioid"
# violin_plot('Density of phenotype_intensity',df_plot, 'phenotype', 'density', 'casetype')

markers = c("CD20", "CD4", "CD8", "FOXP3", "CD68", "CD11c", "CD56", "CK", "Unidentified")

grouptype = "density" #"percent" # #"percent" #"density" # # # #

cols = paste0(markers," ",grouptype)
df_feature = df_feature_ori[,c('MVB','casetype','subtype', cols)]
colnames(df_feature) = c('MVB','casetype','subtype',  markers)
df_feature_long = pivot_longer(df_feature, 4:12, names_to="phenotype", values_to=grouptype)

df_feature_long[[grouptype]][df_feature_long[[grouptype]]==0] = NA #do not plot the density of 0 cores
if (grouptype == "density"){
  df_feature_long['density'] = log(df_feature_long$'density')
}

#---- calculte pvalues for all tissue---------------------------------------------   


library(stringr)
title = paste0(str_to_title(grouptype)," of Phenotype")



get_pvalue <- function(x,y, z=NA, one_sided=TRUE , non_parametric=TRUE ){
  
  if ( ! non_parametric) {
    if (is.na(z)){
      if (one_sided){
        if (median(x, na.rm=TRUE) >= median(y, na.rm=TRUE)){
          alt="greater"
        }else{
          alt="less"
        }
        pvalue = t.test(x, y, alternative=alt)$p.value
      }else{
        pvalue = t.test(x, y)$p.value
      }
    }else{
      a=1
    }
  }else{
    if (is.na(z)){
      if (one_sided){
        if (median(x, na.rm=TRUE) >= median(y, na.rm=TRUE)){
          alt="greater"
        }else{
          alt="less"
        }
        pvalue = wilcox.test(x, y,mu=-0.02, alternative=alt, exact=FALSE,correct=TRUE)$p.value
      }else{
        pvalue = wilcox.test(x, y)$p.value
      }
    }else{
        a=1
    }
  }
 
  return(formatC(pvalue, format = "e", digits = 2))
}


generate_pvalues <- function(df_plot, markers, casetype, groups, grouptype,title, non_parametric=TRUE){
  pvalues = list()   
  for (marker in markers){
    x = df_plot[[grouptype]][(df_plot[[casetype]]==groups[1]) & (df_plot[['phenotype']]==marker)]
    y = df_plot[[grouptype]][(df_plot[[casetype]]==groups[2]) & (df_plot[['phenotype']]==marker)]
    if (length(groups) == 3){
      z = df_plot[[grouptype]][(df_plot[[casetype]]==groups[3]) & (df_plot[['phenotype']]==marker)]
      pvalues[[marker]] = get_pvalue(x,y,z, non_parametric=non_parametric)
    }
      
    else if (length(groups) == 2){
      pvalues[[marker]] = get_pvalue(x,y, non_parametric=non_parametric)
    }
  }   



  if (non_parametric){
    sink( paste0("../../Paper/Figure/BarPlot/", title, "pvalue non_parametric2.csv"))
  }else{
    sink( paste0("../../Paper/Figure/BarPlot/", title, "pvalue2.csv"))
  }
  print(pvalues)
  sink()
}
    


#Pleural and Peritoneal
casetype = NA
hue_order=c("Pleural", "Peritoneal")
generate_pvalues(df_feature_long, markers, "casetype", hue_order, grouptype,title, non_parametric=TRUE)


df_plot=df_feature_long
casetype = "casetype"
groups=hue_order
_