
# 
# library(hrbrthemes)
# library(viridis)
library(tidyverse) #include ggplot2
library(ComplexHeatmap)

rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("utiles.R")


# Waterfall plot

pathout = "./output/waterfall_plot/"

markers =  c('CD8','CK','CD20', 'CD4', 'FOXP3',  'CD68', 'CD11c', 'CD56',  'Unidentified' )

# markers = c("CD20", "CD4", "CD8", "FOXP3", "CD68", "CD11c", "CD56", "CK", "Unidentified") # this is for generating legned for pieplot


data_wide = read.csv("./output/waterfall_plot/waterfall_compostion_by_MVB_wide.csv")

data_wide_sorted = data_wide %>%
                  mutate(CaseType = factor(CaseType, levels=c("Pleural", "Peritoneal")))%>%
                  mutate(subtype = factor(subtype, levels = c("epithelioid", "biphasic", "sarcomatoid"))) %>%
                  group_by(CaseType,subtype)  %>% # groupby will sort group by its factor level
                  arrange(desc(CD8), .by_group = TRUE)

# factor used inner order, like for groupby and plot
# arrange is for visible order

# write.csv(data_wide_sorted, "./output/waterfall_plot/data_wide_sorted.csv")  
  
data_long = data_wide_sorted  %>%
            replace(is.na(.), 0)  %>% #replce null to 0
             pivot_longer(2:10, names_to="Phenotype", values_to="Percent") %>% # select pivot columns to pivot
             mutate( MVB= factor(MVB, levels=data_wide_sorted$MVB))%>% # when plot, the order will change, so need to make it factor to keep the original order when plot
             mutate( Phenotype = factor(Phenotype, 
                                        levels=markers))
# write.csv(data_long, "./output/waterfall_plot/data_long.csv")



# plot

colors = setNames(  #paletteer_d("tidyquant::tq_light")
                    # c("#6A3D9A",
                    #  "#FF7F00",
                    #  "#FDBF6F",
                    #  "#B2DF8A",
                    #  "#1F78B4",
                    #  "#A6CEE3",
                    #  "#18BC9C",
                    #  "#E31A1C",
                    #  "#4B4C4E"),
  
                   # c(#paletteer_d("RColorBrewer::Set1")
                   #   "#E41A1C",
                   #   "#377EB8",
                   #   "#4DAF4A",
                   #   "#984EA3",
                   #   "#FF7F00",
                   #   "#FFFF33",
                   #   "#A65628",
                   #   "#F781BF",
                   #   "#999999"
                   # ),
                   
                  sapply(markers, get_marker_color),
                  markers)

# my_labs <- c( expression(CD8^{"+"}~'T cells'),'Pan CK cells','B cells (CD20)',expression(CD4^{"+"}~'T cells'),
#               'Treg cells (FOXP3)','Macrophages (CD68)', 'Dendritic cells (CD11c)', 'NK cells (CD56)','Unidentified cells' )

my_labs = as.vector(sapply(markers, get_marker_label))
g1 = ggplot(data_long) +
  geom_bar(mapping =  aes( x=MVB,y=Percent, fill=Phenotype), 
    position="stack", stat="identity") + 
  # scale_fill_viridis(discrete=TRUE,option = "D", begin = 0, end = 1,direction = -1, name="") +
  # theme_ipsum() +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(), 
        panel.grid.major.x = element_blank(),
        legend.text.align=0,
        legend.title = element_text(face = "bold") ,
        panel.background = element_blank(), # remove backgroud
        text=element_text(family="sans")
  )+
  
  #left aligh legend
  scale_fill_manual(values = colors, #change the legend color and labels
                    labels=my_labs,
                    name="Cell type") +
  labs(x = "Patient", y= "Cell type composition")
  
  
 
print(g1)

pdf(paste0(pathout, "waterfall_main2.pdf"), width = 8, height = 2.5) #height = 3
# png(paste0(pathout, "waterfall_main.png"), width = 8, height = 3,units="in", res=600, family="Arial")

# pdf(paste0(pathout, "legend_forPiePlot.pdf"), width = 8, height = 3)
# png(paste0(pathout, "legend_forPiePlot.png"), width = 8, height = 3,units="in", res=600, family="Arial")

print(g1)
dev.off()

# png(paste0(pathout, "waterfall_main.png"), width = 8, height = 3,units="in", res=600, family="Arial")
# print(g1)
# dev.off()

# ggsave(g1, width = 8, height = 3, dpi = 500,  file=paste0(pathout, "waterfall_main.png"))

# genereate the annotation part

data_heatmap  = data_wide_sorted %>%
  column_to_rownames(var="MVB")  %>%
  select(-c(CaseType,subtype))%>%
  t()

pushViewport(viewport(gp = gpar(fontfamily = "sans")))
ha = HeatmapAnnotation("Case type" =data_wide_sorted$CaseType, 
                       Subtype = data_wide_sorted$subtype,
                       col = list(
                         "Case type" = c("Pleural" = "#8f6798", "Peritoneal" = "#fef17c"),
                         Subtype = c("epithelioid" = "#8da0cb", "biphasic" = "#7db988", "sarcomatoid" = "#e59d76")
                         ),
                                  
                       annotation_name_side = "left",
                       annotation_name_gp= gpar(fontsize = 9)
                       
    )



ht = Heatmap (data_heatmap,
                 top_annotation = ha,
                 cluster_rows = FALSE,
                 cluster_columns = FALSE,
                  show_heatmap_legend = FALSE
                 )

draw(ht, newpage = FALSE)


png(file=paste0(pathout, "waterfall_main_anno.png"), width=8.39,height=3.25,units="in",res=600, family="Arial")
print(ht)
dev.off()

pdf(file=paste0(pathout, "waterfall_main_anno.pdf"), width=8.39,height=3.25) #7.55
print(ht)
dev.off()

ggsave(ht, width = 12, height = 6, dpi = 500, file=paste0(pathout, "waterfall_main_anno.png"))


###############################
# using pheatmap
library(pheatmap)
colors2 = list(
  "Case type" =  c("Pleural" = "#8f6798", "Peritoneal" = "#fef17c"),
  "Subtype" = c("epithelioid" = "#8da0cb", "biphasic" = "#7db988", "sarcomatoid" = "#e59d76"))
  
annotation_col = data.frame(
         "Case type" = data_wide_sorted$CaseType,
        "Subtype" = data_wide_sorted$subtype,
        check.names = FALSE)
rownames(annotation_col) = colnames(data_heatmap)

png(paste0(pathout, "waterfall_main_anno1.png"), width = 12, height = 6, units="in", res=600, family="Arial")
g = pheatmap(data_heatmap,border_color=NA, annotation_col  = annotation_col, annotation_colors = colors2,cluster_rows=FALSE, cluster_cols=FALSE, annotation_names_row=F,show_colnames=FALSE) 
dev.off()

pdf(paste0(pathout, "waterfall_main_anno1.pdf"), width = 12, height = 6)
g = pheatmap(data_heatmap,border_color=NA, annotation_col  = annotation_col, annotation_colors = colors2,cluster_rows=FALSE, cluster_cols=FALSE, annotation_names_row=F,show_colnames=FALSE) 
dev.off()

                  