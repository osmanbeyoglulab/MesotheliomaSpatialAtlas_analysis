library(stats)
library(data.table)
library(tidyr)
library("survival")
library("survminer")
library(km.ci)
library("coin") #logrank_testt
library(ggplot2)
library(dplyr)

rm(list=ls())
# setwd(system("pwd", intern = T) )
dir = dirname(rstudioapi::getSourceEditorContext()$path); print(dir)  # [1] "/Users/lees18/OneDrive - UPMC/S86_Covid19_scOmics_Xiaosjun/3_SeuratObject_ClusteringByPatient"
setwd(dir)

args = commandArgs(trailingOnly=TRUE)


map = fread("../Data_folder/TMA1/mapping.csv")
map = as.data.frame(map)
map = map[map$Classification == "Malignant", ]
colnames(map)[colnames(map) == 'Annotation ID'] <- 'core.ID'

df_core_feature = read.csv("./output/cellContact/core_features_allMarkers_withIntensity_core.csv", check.names = TRUE)

# log density
markers = c("CD4", "CD8", "CK", "CD20", "FOXP3", "CD68", "CD11c", "CD56")
for (marker in markers){
  df_core_feature[,paste0(marker,".density")] = log(df_core_feature[,paste0(marker,".density")]  + .Machine$double.eps)
}

cores = intersect(df_core_feature$core.ID, map$core.ID)
df_core_feature = df_core_feature[df_core_feature$core.ID %in% cores,]
df_core_feature =  df_core_feature[match(cores, df_core_feature$core.ID),]

map_sel = map[map$core.ID %in% cores, c('core.ID','SurvivalPeriod', 'Grade')]
map_sel = map_sel[match(cores, map_sel$core.ID), ]

df_merged = df_core_feature
df_merged$'SurvivalPeriod' = map_sel$'SurvivalPeriod'
df_merged$'Grade' = map_sel$'Grade'
# df_merged = merge(x = map, y = df_core_feature, by = "core.ID", all = TRUE)

df_merged = df_merged[df_merged$SurvivalPeriod != "Unknown",]
# df_merged = df_merged %>% drop_na()
df_merged['status'] = 1

addGrade = FALSE
tumorcore = FALSE
if (tumorcore){
  df_merged_sub = df_merged[df_merged$CK.percent > 0.2, ]
}else{
  df_merged_sub = df_merged
}
  
# survival analysis based on patient instead of core, so get the mean value of the cores for each patient
df_survival <- df_merged_sub %>% group_by(MVB, Grade, SurvivalPeriod, casetype, status) %>% 
    summarise(mean_CD4.density=mean(CD4.density,na.rm=TRUE),
              mean_CD8.density=mean(CD8.density,na.rm=TRUE),
              mean_CD20.density=mean(CD20.density,na.rm=TRUE),
              mean_CD68.density=mean(CD68.density,na.rm=TRUE),
              mean_FOXP3.density=mean(FOXP3.density,na.rm=TRUE),
              mean_CD11c.density=mean(CD11c.density,na.rm=TRUE),
              mean_CD56.density=mean(CD56.density,na.rm=TRUE),
              mean_CK.density=mean(CK.density,na.rm=TRUE),
              mean_BAP1.density = mean(BAP1.density, na.rm=TRUE),
              mean_NF2.density = mean(NF2.density, na.rm=TRUE),
              mean_MTAP.density = mean(MTAP.density, na.rm=TRUE),
              mean_LAG3.density = mean(LAG3.density, na.rm=TRUE),
              mean_CD4.intensity=mean(CD4.intensity,na.rm=TRUE),
              mean_CD8.intensity=mean(CD8.intensity,na.rm=TRUE),
              mean_CD20.intensity=mean(CD20.intensity,na.rm=TRUE),
              mean_CD68.intensity=mean(CD68.intensity,na.rm=TRUE),
              mean_FOXP3.intensity=mean(FOXP3.intensity,na.rm=TRUE),
              mean_CD11c.intensity=mean(CD11c.intensity,na.rm=TRUE),
              mean_CD56.intensity=mean(CD56.intensity,na.rm=TRUE),
              mean_CK.intensity=mean(CK.intensity,na.rm=TRUE),
              mean_BAP1.intensity = mean(BAP1.intensity, na.rm=TRUE),
              mean_NF2.intensity = mean(NF2.intensity, na.rm=TRUE),
              mean_MTAP.intensity = mean(MTAP.intensity, na.rm=TRUE),
              mean_LAG3.intensity = mean(LAG3.intensity, na.rm=TRUE),
              mean_CK.percent=mean(CK.percent,na.rm=TRUE),
              mean_CD4.percent=mean(CD4.percent,na.rm=TRUE),
              mean_CD8.percent=mean(CD8.percent,na.rm=TRUE),
              mean_CD20.percent=mean(CD20.percent,na.rm=TRUE),
              mean_CD68.percent=mean(CD68.percent,na.rm=TRUE),
              mean_FOXP3.percent=mean(FOXP3.percent,na.rm=TRUE),
              mean_CD11c.percent=mean(CD11c.percent,na.rm=TRUE),
              mean_CD56.percent=mean(CD56.percent,na.rm=TRUE),
              mean_BAP1.percent = mean(BAP1.percent, na.rm=TRUE),
              mean_NF2.percent = mean(NF2.percent, na.rm=TRUE),
              mean_MTAP.percent = mean(MTAP.percent, na.rm=TRUE),
              mean_LAG3.percent = mean(LAG3.percent, na.rm=TRUE),
              .groups = 'drop') %>%
  as.data.frame()


# df_type <- df_type %>% group_by(MVB, Grade, SurvivalPeriod, casetype, status) %>% 
#   summarise(mean_CD4.density=mean(CD4.density,na.rm=TRUE),
#             mean_CD8.density=mean(CD8.density,na.rm=TRUE),
#             mean_CD20.density=mean(CD20.density,na.rm=TRUE),
#             mean_CD68.density=mean(CD68.density,na.rm=TRUE),
#             mean_FOXP3.density=mean(FOXP3.density,na.rm=TRUE),
#             mean_CD11c.density=mean(CD11c.density,na.rm=TRUE),
#             mean_CD56.density=mean(CD56.density,na.rm=TRUE),
#             mean_CK.density=mean(CK.density,na.rm=TRUE),
#             mean_BAP1.density = mean(BAP1.density, na.rm=TRUE),
#             mean_NF2.density = mean(NF2.density, na.rm=TRUE),
#             mean_MTAP.density = mean(MTAP.density, na.rm=TRUE),
#             mean_LAG3.density = mean(LAG3.density, na.rm=TRUE),
#             mean_CD4.intensity=mean(CD4.intensity,na.rm=TRUE),
#             mean_CD8.intensity=mean(CD8.intensity,na.rm=TRUE),
#             mean_CD20.intensity=mean(CD20.intensity,na.rm=TRUE),
#             mean_CD68.intensity=mean(CD68.intensity,na.rm=TRUE),
#             mean_FOXP3.intensity=mean(FOXP3.intensity,na.rm=TRUE),
#             mean_CD11c.intensity=mean(CD11c.intensity,na.rm=TRUE),
#             mean_CD56.intensity=mean(CD56.intensity,na.rm=TRUE),
#             mean_CK.intensity=mean(CK.intensity,na.rm=TRUE),
#             mean_BAP1.intensity = mean(BAP1.intensity, na.rm=TRUE),
#             mean_NF2.intensity = mean(NF2.intensity, na.rm=TRUE),
#             mean_MTAP.intensity = mean(MTAP.intensity, na.rm=TRUE),
#             mean_LAG3.intensity = mean(LAG3.intensity, na.rm=TRUE),
#             mean_CK.percent=mean(CK.percent,na.rm=TRUE),
#             mean_CD4.percent=mean(CD4.percent,na.rm=TRUE),
#             mean_CD8.percent=mean(CD8.percent,na.rm=TRUE),
#             mean_CD20.percent=mean(CD20.percent,na.rm=TRUE),
#             mean_CD68.percent=mean(CD68.percent,na.rm=TRUE),
#             mean_FOXP3.percent=mean(FOXP3.percent,na.rm=TRUE),
#             mean_CD11c.percent=mean(CD11c.percent,na.rm=TRUE),
#             mean_CD56.percent=mean(CD56.percent,na.rm=TRUE),
#             mean_BAP1.percent = mean(BAP1.percent, na.rm=TRUE),
#             mean_NF2.percent = mean(NF2.percent, na.rm=TRUE),
#             mean_MTAP.percent = mean(MTAP.percent, na.rm=TRUE),
#             mean_LAG3.percent = mean(LAG3.percent, na.rm=TRUE),
#             .groups = 'drop') %>%
#   as.data.frame()
df_survival$SurvivalPeriod = as.integer(df_survival$SurvivalPeriod)

#  markrs
# markers = c("CD4", "CD8", "CK", "FOXP3", "CD20", "CD68", "CD11c", "CD56") #, ) #, 
# markers = c("CD4", "CD8", "CK",  "FOXP3", "CD68") #"CD20",

markers = c("CD4", "CD8", "CK", "FOXP3", "CD20", "CD68", "CD11c", "CD56")
#########################################################################
# pleural vs peritoneal
pathout = "../../Paper/Figure/Survival/"

df = df_survival[df_survival$casetype != "Other", ]
my.fit <- survfit(Surv(as.numeric(SurvivalPeriod), status) ~ casetype   , data=df)

g = ggsurvplot(my.fit,
              pval = TRUE, conf.int = FALSE,
              risk.table = TRUE, # Add risk table
              risk.table.col = "strata", # Change risk table color by groups
              linetype = "strata", # Change line type by groups
              # surv.median.line = "hv", # Specify median survival
              ggtheme = theme_bw()) # Change ggplot2 theme
              # palette = c("#E7B800", "#2E9FDF")))

pdf(paste0(pathout, "suvival_curves_Pleural_vs_Pritoneal.pdf"))
print(g, newpage = FALSE)
dev.off()


pdf(paste0(pathout, "suvival_curves_Pleural_vs_Pritoneal_noTable.pdf"))
print(g$plot, newpage = FALSE)
dev.off()
#########################################################################
pathout  = "../../Paper/Figure/NEW_ORGANIZATION/Survival/Phenotype/"
vars = c()
hazards = c()
conf_lows = c()
conf_highs = c()
pvalues = c()
casetypes = c()

addGrade = FALSE
tumorcore = FALSE

counts = list()

for (casetype in c("All", "Pleural", "Peritoneal")){
  if (casetype == "All"){
    df_sub = df_survival
  }else{
    df_sub = df_survival[ (df_survival$casetype ==casetype),]
  }
  counts[[casetype]] = dim(df_sub)[1]
  # for (grouptype in c( ".density")){ #, ".intensity",".percent"
  grouptype = ".percent" #".density" ".intensity" # #
  
  for (marker in markers){
    groupname = paste0("mean_",marker, grouptype)
    # replace 0 with NA ( don't count 0 patients)
    
    df_sub[[groupname]][df_sub[[groupname]]==0]=NA
    
    if(addGrade){
      res <- coxph( as.formula(paste0('Surv(SurvivalPeriod, status) ~ ',groupname, "+ Grade")), data = df_sub)
    }else{
      res <- coxph( as.formula(paste0('Surv(SurvivalPeriod, status) ~ ',groupname)), data = df_sub)
    }
    print(res)
    if (grouptype == ".density"){
      vars = c(vars, paste0("log(", marker, ")(", counts[[casetype]], ")") ) #    c(vars, marker)
    }else{
      vars = c(vars, paste0(marker, "(", counts[[casetype]], ")") ) #    c(vars, marker)
    }
    hazards = c(hazards, summary(res)$coefficient[2])
    pvalues = c(pvalues, summary(res)$coefficient[5])
    conf_lows = c(conf_lows, summary(res)$conf.int[3]) 
    conf_highs = c(conf_highs, summary(res)$conf.int[4]) 
    casetypes = c(casetypes, casetype)

  }
}

pvalues = formatC(pvalues, format = "e", digits = 2)
hazards = formatC(hazards, format = "e", digits = 2) #format(round(hazards, 2), nsmall = 2)
conf_lows = formatC(conf_lows, format = "e", digits = 2) #format(round(conf_lows, 2), nsmall = 2)
conf_highs = formatC(conf_highs, format = "e", digits = 2) #format(round(conf_highs, 2), nsmall = 2)

hazards_conf = paste0(hazards, " (", conf_lows, ", ", conf_highs, ")")
df = as.data.frame(cbind(vars, hazards_conf, pvalues, casetypes))
df <- apply(df,2,as.character)


if (tumorcore){
  filename = paste0("markers_hazard_ratio_tumorcores",grouptype)
}else{
  filename = paste0("markers_hazard_ratio_allcores", grouptype)
}
if (addGrade){
  filename = paste0(filename,  "_addGrade")
}
write.csv(df, paste0(pathout,filename,".csv") )
 

    
    
########################################################################
# survfit  discritize phenotype
tumorcore = FALSE

# for (grouptype in c(".intensity", ".percent", ".density")){
  
  grouptype = ".percent"
  if (tumorcore){
    outpath = paste0("./output/survival2/tumorcores/")
  }else{
    outpath = paste0("./output/survival2/allcores/")
  }
  dir.create(outpath, showWarnings = FALSE, recursive = TRUE)

  for (casetype in c("All", "Pleural", "Peritoneal")){
    if (casetype == "All"){
      df_sub = df_survival
    }else{
      df_sub = df_survival[ (df_survival$casetype ==casetype),]
    }
    # for (grouptype in c( ".density")){ #, ".intensity",".percent"

    pvalues = c()
    vars = c()
    for (marker in markers){
      groupname = paste0("mean_",marker, grouptype)

      df_plot = df_sub
      df_plot[[groupname]][df_plot[[groupname]]==0] = NA
      markerSep = quantile( df_plot[[groupname]], c(0.4, 0.6), na.rm=TRUE)
      # markerSep = quantile( df_plot[groupname], c(0.3, 0.7), na.rm=TRUE)
      # markerSep = c(mean(df_plot[,groupname]), mean(df_plot[,groupname]))
      # markerSep = quantile( df_plot[groupname], c(0.4, 0.6), na.rm=TRUE)
      highCond = which(df_plot[groupname] >= markerSep[2])
      df_plot[highCond, 'group'] = paste0(marker ,"+")
      lowCond = which(df_plot[groupname] < markerSep[1])
      df_plot[lowCond, 'group'] = paste0(marker ,"-")
      
      df_plot = df_plot[! is.na(df_plot$'group'),]

      my.fit <- survfit(Surv(as.numeric(SurvivalPeriod), status) ~ group   , data=df_plot)
      ptext = surv_pvalue(my.fit)$pval
      ptext = formatC(ptext, format = "e", digits = 2)

      
      pvalues = c(pvalues, ptext)
      vars = c(vars, marker)
      # if (casetype == "All"){
      #   title = paste0("Survival ", marker, " high vs low in all patients ")
      # }else{
      #   title = paste0("Survival ", marker, " high vs low in ", casetype)
      # }
      # g = ggsurvplot(my.fit,
      #            pval = TRUE, conf.int = FALSE,
      #            title = title,
      #            risk.table = TRUE, # Add risk table
      #            risk.table.col = "strata", # Change risk table color by groups
      #            linetype = "strata", # Change line type by groups
      #            # surv.median.line = "hv", # Specify median survival
      #            ggtheme = theme_bw()) # Change ggplot2 theme
      #            # palette = c("#E7B800", "#2E9FDF")))
      #
      # print(g)
    }
    
    
    
    
    df= cbind.data.frame(vars, pvalues)
    write.csv(df, paste0(outpath, "survfit_", grouptype,"_",casetype, "_0.4_0.6.csv") )

  }
# }

plots.dir.path <- list.files(tempdir(), pattern="rs-graphics", full.names = TRUE);
plots.png.paths <- list.files(plots.dir.path, pattern=".png", full.names = TRUE)

file.copy(from=plots.png.paths, to="./output/survival/plots/")


#   
#     if (casetype == "All"){
#       title = paste0(marker,grouptype, " for all patients")
#     }else{
#       title = paste0(marker,grouptype, " for ", casetype)
#     }
#     if(tumorcore){
#       title = paste0(title, " all cores")
#     }else{
#       title = paste0(title, " tumore cores")
#     }
#     
#     
#     # cols = c("Magenta","dodger blue")
#     # ltys = c(1,4)
#     # plot(my.fit, lwd = 2.5, lty = ltys,  cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5, xlab = "Months", ylab="Survival Probability", main=title, col=cols, ylim=c(0,1), xlim=c(0,150)) #xmax=100
#     # 
#     # legend(x=55,y=1.03, ptext, box.lwd = 0, box.col = "white",bg = "white" , cex = 1.4) # no lty then no bars on lengend
#     # 
#     # legend(x=60,y=0.85, c(paste0(marker, " high (n=",sum(highCond),")"), paste0(marker, " low (n=",sum(lowCond),")")), cex = 1.3, lty =ltys, col=cols,  box.lwd = 0, box.col = "white")
#     
#     # dev.copy(png,paste0(outpath,title,".png"))
#     # dev.off()
#     
#     #############################
#     df_plot$SurvivalPeriod = as.integer(df_plot$SurvivalPeriod)
#     
#     res <- coxph( Surv(SurvivalPeriod, status) ~ group + Grade, data = df_plot)
#     
#     group_df <- with(df_plot,
#                    data.frame(group = c(1, 2), 
#                               Grade = rep('Low', 2)
#                              
#                    )
#     )
#     fit <- survfit(res , newdata = group_df)
#     ggsurvplot(fit, conf.int = TRUE, legend.labs=paste0(marker, c(" low", " high")),
#                ggtheme = theme_minimal())
#     
#     
#     vars = c(vars, groupname)
#     hazards = c(hazards, summary(res)$coefficient[2])
#     pvalues = c(pvalues, summary(res)$coefficient[5])
#     # pvalues = c(pvalues, summary(res)$sctest[3])
#     conf_lows = c(conf_lows, summary(res)$conf.int[3]) 
#     conf_highs = c(conf_highs, summary(res)$conf.int[4]) 
#     seps = c(seps, sep)
#     tumorcores = c(tumorcores, tumorcore)
#     casetypes = c(casetypes, casetype)
#     grouptypes = c(grouptypes, grouptype)
#     }
# }
# 
# 
# 
# pvalues = formatC(pvalues, format = "e", digits = 2)
# hazards = format(round(hazards, 2), nsmall = 2)
# conf_lows = format(round(conf_lows, 2), nsmall = 2)
# conf_highs = format(round(conf_highs, 2), nsmall = 2)
# hazards_conf = paste0(hazards, " (", conf_lows, ", ", conf_highs, ")")
# df = as.data.frame(cbind(vars, hazards_conf, pvalues, seps, tumorcores, casetypes, grouptypes))
# write.csv(df, "./output/survival2/markers_hazard_ratio_discritisized_CKpercent_density.csv")
# 




# # for Panel2 markers
# markerSep = c(200,200)
# 
# df_plot = df_sub[,c('BAP1.density','NF2.density', 'MTAP.density', 'BAP1.percent','NF2.percent', 'MTAP.percent','SurvivalPeriod')]
# 
# grouptype = 'density'
# # grouptype = 'percent'
# 
# groupname = paste0("BAP1.", grouptype)
# # markerSep = quantile( df_plot[groupname], c(0.1, 0.1), na.rm=TRUE)
# highCond = df_plot[groupname] >= markerSep[2]
# df_plot[highCond, 'BAP1 status'] = "BAP1+"
# lowCond = df_plot[groupname] < markerSep[1]
# df_plot[lowCond, 'BAP1 status'] = "BAP1-"
# 
# groupname = paste0("NF2.", grouptype)
# # markerSep = quantile( df_plot[groupname], c(0.1, 0.1), na.rm=TRUE)
# highCond = df_plot[groupname] >= markerSep[2]
# df_plot[highCond, 'NF2 status'] = "NF2+"
# lowCond = df_plot[groupname] < markerSep[1]
# df_plot[lowCond, 'NF2 status'] = "NF2-"
# 
# groupname = paste0("MTAP.", grouptype)
# # markerSep = quantile( df_plot[groupname], c(0.1, 0.1), na.rm=TRUE)
# highCond = df_plot[groupname] >= markerSep[2]
# df_plot[highCond, 'MTAP status'] = "MTAP+"
# lowCond = df_plot[groupname] < markerSep[1]
# df_plot[lowCond, 'MTAP status'] = "MTAP-"
# 
# 
# df_plot[is.na(df_plot)]=""
# 
# df_plot['group'] = paste0(df_plot$`BAP1 status`, df_plot$`NF2 status`, df_plot$`MTAP status`)
# df_plot['status'] = 1
# 
# my.fit <- survfit(Surv(as.numeric(SurvivalPeriod), status) ~ group, data=df_plot)
# 
# ggsurvplot(my.fit,
#            pval = TRUE, conf.int = FALSE,
#            risk.table = TRUE, # Add risk table
#            risk.table.col = "strata", # Change risk table color by groups
#            linetype = "strata", # Change line type by groups
#            # surv.median.line = "hv", # Specify median survival
#            ggtheme = theme_bw(), # Change ggplot2 theme
#            # palette = c("#E7B800", "#2E9FDF"))
# )
for (i in seq(24)){
  print(format(round(conf_highs[i], 2), nsmall = 2))
}

for (i in seq(24)){
  print(format(round(conf_highs[i], 2), nsmall = 2))
}
