library(ggplot2)
library(ggridges) #geom_density_ridges

#########################################################################################################
min_max_normalize <- function(x,  ...) {
  return((x- min(x,  ...)) /(max(x,  ...)-min(x,  ...)))
}

zscore <- function(x, ...){
  return( (x-mean(x, ...))/sd(x, ...))
}
###########################################################################################################
#scatter plot
get_marker_label_2 = function(marker){
  markers =  c('CD8','CK','CD20', 'CD4', 'FOXP3',  'CD68', 'CD11c', 'CD56',  'Unidentified' )
  labs <- c(expression(CD8^{"+"}~'T cells intensity'),'Tumor cells (Pan-CK) intensity','B cells (CD20) intensity', expression(CD4^{"+"}~'T cells intensity'), 
             'Treg cells (FOXP3) intensity','Macrophages (CD68) intensity', 'Dendritic cells (CD11c) intensity', 'NK cells (CD56) intensity','Unidentified cells intensity' )
  
  return(labs[which(markers==marker)])
}

scatter_plot <- function(df, x, y, color, corr,  title, outpath, save=TRUE){
  
  df[[color]] = factor(df[[color]], levels=c("epithelioid","biphasic", "sarcomatoid", "Multicystic",  "Papillary", "Desmoplastic", "Not specified"))

  labels = c( "epithelioid"="Epithelioid","biphasic"= "Biphasic", "sarcomatoid"="Sarcomatoid", "Multicystic"="Multicystic",  "Papillary"="Papillary", "Desmoplastic"="Desmoplastic", "Not specified"="Not specified")
  colors = c("epithelioid"="#7db988", "biphasic"= "#8da0cb", "sarcomatoid"="#e59d76", "Multicystic"="#DAA520", "Papillary"="#D65B5A","Desmoplastic"="#B47880","Not specified"="#4D4D4D")
  
  scale_color = scale_color_manual(labels = labels, values = colors) 
  
  maxy = max(df[[y]], na.rm=TRUE) 
  ypos = maxy + maxy/20
  
  maxx = max(df[[x]], na.rm=TRUE) 
  xpos = maxx - maxx/5
  
  fig = ggplot(df, aes(x =.data[[x]] , y = .data[[y]], color = .data[[color]]))+
    geom_point(size = 3, #size of point
               alpha = 1, #transparency of point
    )+
    cowplot::theme_cowplot() +
    theme(text=element_text(family="sans")) +

    # scale_color_manual(values=c("#7db988", "#8da0cb", "#e59d76"))+
    scale_color +
    xlab(get_marker_label_2(marker)) +  #(paste(get_marker_label_2(marker),  sub(".", "", grouptype))) + #x-axis title
    ylab(paste(gene, sub(".", "", grouptype))) + #y-axis title
    labs(color='Subtype') + #legend title
    
    annotate("text", x=xpos, y=ypos, label= paste0("corr = ", format(corr, nsmall=2)), size = 13/.pt )
  
  
  print(fig)
  if (save){
    ggsave(fig, width = 5.7, height = 4, dpi = 500, file=paste0(outpath, title, ".pdf"))
  }
}

###########################################################
# ridge density plot

#
ridge_plot_BAP1 <- function(df_plot,title, x, y, xlim,pvalue,num_high, num_low, pathout,annox=4.5, save=TRUE){
  
  g = ggplot(df_plot, aes(x=.data[[x]], y=.data[[y]], fill=.data[[y]]))+
    theme_classic() +
    geom_density_ridges()+
    xlim(xlim) +
    
    scale_fill_manual(values =c("BAP1-"="#008FD5", "BAP1+"="#FF2700" ), # # c("BAP1-"="#5773CC", "BAP1+"="#FFB900" ),
                      labels = c("BAP1+", "BAP1-"),
                      name="",
                      breaks=c("BAP1+","BAP1-") #change the legend order
    )+
    scale_y_discrete(labels = c(paste0("BAP1-\n(n=",num_low,")"), paste0("BAP1+\n(n=",num_high,")") )  ) + #change y ticks
    labs(x = "Cell contact score",
         y= expression("Numver of cells per"~mm^{"2"}) )+
    ggtitle(title)+
    theme( text = element_text(size = 15, family="sans"),
           
           legend.position = c(0.95, 0.9),
           legend.text = element_text(colour="black", size = 13),
           
           plot.title = element_text(hjust = 0.5, #put the title in the middle
                                     margin = margin(0,0,20,0), #increase the space between tilte and plot
                                     size=16,face="bold"),
           
           axis.text.y = element_text( size=13, color="black"), 
           plot.margin = unit(c(1,2,1,1), "lines")) +   #top, right, bottom, left
    
    
    annotate("text", x=annox, y=3.5, label= paste0("pvalue = ", pvalue), size = 15/.pt )
  
  print(g)
  if (save){
    pdf(paste0(pathout, title,".pdf"), width = 8.5, height = 3)
    print(g)
    dev.off()
    print(paste0("file save to", pathout, title,".pdf"))
  }
}  


ridge_plot_casetype <- function(df_plot,title, x, y, xlim,pvalue,num_high, num_low, pathout,annox=4.5, save=TRUE){
  
  g = ggplot(df_plot, aes(x=.data[[x]], y=.data[[y]], fill=.data[[y]]))+
    theme_classic() +
    geom_density_ridges()+
    xlim(xlim) +
    
    # scale_fill_manual(values = c("BAP1-"="#5773CC", "BAP1+"="#FFB900" ),
    #                   labels = c("BAP1+", "BAP1-"),
    #                   name="",
    #                   breaks=c("BAP1+","BAP1-") #change the legend order
    # )+
    # scale_y_discrete(labels = c(paste0("BAP1-\n(n=",num_low,")"), paste0("BAP1+\n(n=",num_high,")") )  ) + #change y ticks
    labs(x = "Cell contact score",
         y= expression("Numver of cells per"~mm^{"2"}) )+
    ggtitle(title)+
    theme( text = element_text(size = 15, family="sans"),
           
           legend.position = c(0.95, 0.9),
           legend.text = element_text(colour="black", size = 13),
           
           plot.title = element_text(hjust = 0.5, #put the title in the middle
                                     margin = margin(0,0,20,0), #increase the space between tilte and plot
                                     size=16,face="bold"),
           
           axis.text.y = element_text( size=13, color="black"), 
           plot.margin = unit(c(1,2,1,1), "lines")) +   #top, right, bottom, left
    
    
    annotate("text", x=annox, y=3.5, label= paste0("pvalue = ", pvalue), size = 15/.pt )
  
  print(g)
  if (save){
    pdf(paste0(pathout, title,".pdf"), width = 8.5, height = 3)
    print(g)
    dev.off()
  }
} 
##########################################################################

ridge_plot = function(df_plot,title, x, y, xlim,pvalue,num_pl, num_pe, posx, pathout, fill_cor, save=TRUE){
  
  g = ggplot(df_plot, aes(x=.data[[x]], y=.data[[y]], fill=.data[[y]]))+
    theme_classic() +
    geom_density_ridges()+
    xlim(xlim) +
    
    scale_fill_manual(values = fill_cor)+
    #                   labels = c("BAP1+", "BAP1-"),
    #                   name="",
    #                   breaks=c("BAP1+","BAP1-") #change the legend order
    # )+
    # scale_y_discrete(labels = c(paste0("BAP1-\n(n=",num_low,")"), paste0("BAP1+\n(n=",num_high,")") )  ) + #change y ticks
    labs(x = "Cell contact score",
         y= expression("Numver of cells per"~mm^{"2"}) )+
    ggtitle(title)+
    theme( text = element_text(size = 15),
           
           legend.position = c(0.95, 0.9),
           legend.text = element_text(colour="black", size = 13),
           
           plot.title = element_text(hjust = 0.5, #put the title in the middle
                                     margin = margin(0,0,20,0), #increase the space between tilte and plot
                                     size=16,face="bold"),
           
           axis.text.y = element_text( size=13, color="black"), 
           plot.margin = unit(c(1,2,1,1), "lines"),
           text=element_text(family="sans")) +   #top, right, bottom, left
    
    
    annotate("text", x=posx, y=2.9, label= paste0("pvalue = ", pvalue), size = 15/.pt ) +
     annotate("text", x=posx, y=2.5, label= paste0("N count (", num_pl, ",", num_pe, ")"), size = 15/.pt )
  
  print(g)
  if (save){
    png(paste0(pathout, title,".png"), width = 8.5, height = 3, )
    print(g)
    dev.off()
  }
}  

######################################################
# dot plot

two_group_dotplot = function(df_counts, clinic){
my_labs <- c( expression(CD4^{"+"}~'T cells'), expression(CD8^{"+"}~'T cells'),'Pan CK cells','B cells (CD20)',
              'Macrophages (CD68)', 'Treg cells (FOXP3)','Dendritic cells (CD11c)', 'NK cells (CD56)')
if (clinic == "c"){
  scale_colors = scale_color_manual(labels = c(
    
    paste0("Exposure (n=",df_counts[df_counts$catgys=='AsbestosExposure','Ng1'], ')'),
    paste0("Smoker (n=",df_counts[df_counts$catgys=='smoking','Ng1'], ')'),
    paste0("Male (n=",df_counts[df_counts$catgys=='Gender','Ng1'], ')'),
    paste0("White (n=",df_counts[df_counts$catgys=='Race','Ng1'], ')'),
    paste0("<=60 (n=",df_counts[df_counts$catgys=='agesplit','Ng1'], ')'),
    paste0("High (n=",df_counts[df_counts$catgys=='Grade','Ng1'], ')'),
    paste0("<=T1 (n=",df_counts[df_counts$catgys=='cTStage','Ng1'], ')'),
    paste0("Lymph- (n=",df_counts[df_counts$catgys=='cNStage','Ng1'], ')'),
    paste0("no Spread (n=",df_counts[df_counts$catgys=='cMStage','Ng1'], ')'),
    
    paste0("not Exposure (n=",df_counts[df_counts$catgys=='AsbestosExposure','Ng2'], ')'),
    paste0("not Smoker (n=",df_counts[df_counts$catgys=='smoking','Ng2'], ')'),
    paste0("Female (n=",df_counts[df_counts$catgys=='Gender','Ng2'], ')'),
    paste0("not White (n=",df_counts[df_counts$catgys=='Race','Ng2'], ')'),
    paste0(">60 (n=",df_counts[df_counts$catgys=='agesplit','Ng2'], ')'),
    paste0("not High (n=",df_counts[df_counts$catgys=='Grade','Ng2'], ')'),
    paste0(">T1 (n=",df_counts[df_counts$catgys=='cTStage','Ng2'], ')'),
    paste0("Lymph+ (n=",df_counts[df_counts$catgys=='cNStage','Ng2'], ')'),
    paste0("Spread (n=",df_counts[df_counts$catgys=='cMStage','Ng2'], ')')
  ), 
  values = c('Asbestos exposure1' = rgb(254, 123, 125, maxColorValue = 255),
             Smoking1= rgb(255, 181, 87, maxColorValue = 255),
             Gender1 = rgb(255, 222, 31, maxColorValue = 255),
             Race1 = rgb(191, 181, 54, maxColorValue = 255),
             'Age split1' = rgb(147, 198, 139, maxColorValue = 255),
             Grade1 = rgb(34, 187, 127, maxColorValue = 255),
             cTStage1 = rgb(76, 182, 244, maxColorValue = 255),
             cNStage1 = rgb(173, 126, 247, maxColorValue = 255),
             cMStage1 = rgb(229, 109, 220, maxColorValue = 255),
             
             'Asbestos exposure0' =rgb(204, 73, 75, maxColorValue = 255),
             Smoking0 =rgb(205, 131, 57, maxColorValue = 255),
             Gender0 = rgb(205, 172, 20, maxColorValue = 255),
             Race0 = rgb(131, 121, 4, maxColorValue = 255),
             'Age split0' =rgb(97, 148, 89, maxColorValue = 255),
             Grade0 = rgb(0, 137, 77, maxColorValue = 255),
             cTStage0 = rgb(26, 132, 194, maxColorValue = 255),
             cNStage0 = rgb(123, 76, 197, maxColorValue = 255),
             cMStage0 = rgb(179, 59, 170, maxColorValue = 255)))
}else if(clinic=="p"){
  scale_colors = scale_color_manual(labels = c(
    
    paste0("Exposure (n=",df_counts[df_counts$catgys=='AsbestosExposure','Ng1'], ')'),
    paste0("Smoker (n=",df_counts[df_counts$catgys=='smoking','Ng1'], ')'),
    paste0("Male (n=",df_counts[df_counts$catgys=='Gender','Ng1'], ')'),
    paste0("White (n=",df_counts[df_counts$catgys=='Race','Ng1'], ')'),
    paste0("<=60 (n=",df_counts[df_counts$catgys=='agesplit','Ng1'], ')'),
    paste0("High (n=",df_counts[df_counts$catgys=='Grade','Ng1'], ')'),
    paste0("<=T1 (n=",df_counts[df_counts$catgys=='pTStage','Ng1'], ')'),
    paste0("Lymph- (n=",df_counts[df_counts$catgys=='pNStage','Ng1'], ')'),
    paste0("no Spread (n=",df_counts[df_counts$catgys=='pMStage','Ng1'], ')'),
    
    paste0("not Exposure (n=",df_counts[df_counts$catgys=='AsbestosExposure','Ng2'], ')'),
    paste0("not Smoker (n=",df_counts[df_counts$catgys=='smoking','Ng2'], ')'),
    paste0("Female (n=",df_counts[df_counts$catgys=='Gender','Ng2'], ')'),
    paste0("not White (n=",df_counts[df_counts$catgys=='Race','Ng2'], ')'),
    paste0(">60 (n=",df_counts[df_counts$catgys=='agesplit','Ng2'], ')'),
    paste0("not High (n=",df_counts[df_counts$catgys=='Grade','Ng2'], ')'),
    paste0(">T1 (n=",df_counts[df_counts$catgys=='pTStage','Ng2'], ')'),
    paste0("Lymph+ (n=",df_counts[df_counts$catgys=='pNStage','Ng2'], ')'),
    paste0("Spread (n=",df_counts[df_counts$catgys=='pMStage','Ng2'], ')')
  ), 
  values = c('Asbestos exposure1' = rgb(254, 123, 125, maxColorValue = 255),
             Smoking1= rgb(255, 181, 87, maxColorValue = 255),
             Gender1 = rgb(255, 222, 31, maxColorValue = 255),
             Race1 = rgb(191, 181, 54, maxColorValue = 255),
             'Age split1' = rgb(147, 198, 139, maxColorValue = 255),
             Grade1 = rgb(34, 187, 127, maxColorValue = 255),
             pTStage1 = rgb(76, 182, 244, maxColorValue = 255),
             pNStage1 = rgb(173, 126, 247, maxColorValue = 255),
             pMStage1 = rgb(229, 109, 220, maxColorValue = 255),
             
             'Asbestos exposure0' =rgb(204, 73, 75, maxColorValue = 255),
             Smoking0 =rgb(205, 131, 57, maxColorValue = 255),
             Gender0 = rgb(205, 172, 20, maxColorValue = 255),
             Race0 = rgb(131, 121, 4, maxColorValue = 255),
             'Age split0' =rgb(97, 148, 89, maxColorValue = 255),
             Grade0 = rgb(0, 137, 77, maxColorValue = 255),
             pTStage0 = rgb(26, 132, 194, maxColorValue = 255),
             pNStage0 = rgb(123, 76, 197, maxColorValue = 255),
             pMStage0 = rgb(179, 59, 170, maxColorValue = 255)))
}


fig = df_dotplot1 %>% 
  ggplot(aes(x=phenotypes, y = items, size=sig , color=combine) )+ 
  geom_point() + 
  cowplot::theme_cowplot() + 
  # theme(axis.line  = element_blank()) +
  theme(axis.text.x = element_text(angle = 45, hjust=1, size=8)) +
  theme(axis.text.y = element_text(size=8)) +
  xlab('') +
  ylab('') +
  
  scale_x_discrete(labels = my_labs) +
  scale_size(range = c(0, 4+plus), labels=c("not significant", "P < 0.05", "P < 0.005", "P < 0.0005"), name="Pvalue") +
  
  scale_colors +
  guides( colour = guide_legend(ncol=2)) + # change the lengend order
  theme(legend.text=element_text(size=8),
        legend.title=element_blank()) +
  theme(legend.key.height=unit(0.85,'cm')) +
  theme(text=element_text(family="sans")) +
  # guides(color = FALSE)
  guides(size = FALSE)


print(fig)
ggsave(fig, width = 8, height = 4.5, dpi = 500, file=paste0(pathout, "two_group_sep_dotplot_marker", grouptype,"_",clinic,"_stage_sep1.pdf"))
ggsave(fig, width = 8, height = 4.5, dpi = 500, file=paste0(pathout, "two_group_sep_dotplot_marker", grouptype,"_",clinic,"_stage_sep1.png"))
}

get_marker_label = function(marker){
  markers =  c('CD8','CK','CD20', 'CD4', 'FOXP3',  'CD68', 'CD11c', 'CD56',  'Unidentified' )
  labs <- c( expression(CD8^{"+"}~'T cells'),'Tumor cells (Pan-CK)','B cells (CD20)',expression(CD4^{"+"}~'T cells'), 
                'Treg cells (FOXP3)','Macrophages (CD68)', 'Dendritic cells (CD11c)', 'NK cells (CD56)','Unidentified cells' )

  return(labs[which(markers==marker)])
}

get_marker_color = function(marker){
  markers =  c('CD8','CK','CD20', 'CD4', 'FOXP3',  'CD68', 'CD11c', 'CD56',  'Unidentified' )
  colors =  c(#paletteer_d("RColorBrewer::Set1")
    "#E41A1C",
    "#377EB8",
    "#4DAF4A",
    "#984EA3",
    "#FF7F00",
    "#FFFF33",
    "#A65628",
    "#F781BF",
    "#999999"
  )
  return(colors[which(markers==marker)])
}


#################################################
# used in the similarity_withn_among_patient.R
box_plot_similarity <- function(df, x, y, pvalue, title, pathout){
  
  g = ggplot(df, aes_string(x=x, y= y, fill=x)) +
    geom_boxplot( width = 0.7,   #box width
                  lwd=0.4,       #   line width of the box
                  outlier.size = 0.07)  +   #outlier dot size  #outlier.shape = 2, can chage the shape of outlier
    geom_jitter(color="black", size=0.2, alpha=0.2, width = 0.2) +
    
    theme_bw()+
    
    theme(axis.text.x = element_text(angle = 0, vjust = 0, hjust=0.5, size=7)) +
    theme(axis.text.y = element_text(size=7)) +
    theme(plot.title = element_text(hjust = 0.5)) # ,size=8, face = "bold")) +
    theme(text=element_text(family="sans"))+
    My_Theme +
      
    ggtitle(paste("Pvalue = ",pvalue)) +
    ylab("Corr of density distr.of cores ") +
    xlab("")

   
  
  print(g)
  ggsave(g, width = 5, height = 4, file=paste0(pathout, title,".pdf"))
  ggsave(g, width = 5, height = 4, dpi = 500, file=paste0(pathout, title,".png"))
  dev.off()
  
}


singleboxplot <- function(df, marker, x, y, grouptype,type, pvalue, opt, pathout){
  
  df$casetype =  factor(df$casetype, levels=c("Pleural", "Peritoneal"))
  title = paste0(type, grouptype,".",marker,"." , opt)
  num_pl = length(x)
  num_pe = length(y)
  
  maxy = max(df[ paste0(marker, grouptype)], na.rm=TRUE) 
  ylim = maxy + maxy/5
  ypos = maxy + maxy/6
  
  colorder = c("#8f6798", "#fef17c")
  labels = c(paste0("Pleural\n(n=",num_pl,")"), paste0("Peritoneal\n(n=",num_pe,")"))
  
  g = ggplot(df, aes_string(x="casetype", y= paste0(marker, grouptype), fill="casetype")) +

    geom_boxplot( width = 0.7,   #box width
                  lwd=0.4,       #   line width of the box
                  outlier.size = 0.07)  +   #outlier dot size  #outlier.shape = 2, can chage the shape of outlier
    geom_jitter(color="black", size=0.2, alpha=0.2, width = 0.2) +
    
    theme_bw()+
    
    theme(axis.text.x = element_text(angle = 0, vjust = 0, hjust=0.5, size=7)) +
    theme(axis.text.y = element_text(size=7)) +
    theme(plot.title = element_text(hjust = 0.5,size=8, face = "bold")) +
    My_Theme +
    
    ggtitle(paste(type,marker, opt)) +
    ylab(tools::toTitleCase(substring(grouptype, 2, nchar(grouptype)))) +
    xlab("")+
    ylim (0, ylim) +
    
    scale_fill_manual( values = colorder, 
                       labels = labels) +
    
    # scale_fill_manual( values = colorder, 
    #                    labels = c("Pleural" = paste0("Pleural\n(n=",num_pl,")"), "Peritoneal" =  paste0("Peritoneal\n(n=",num_pe,")"))) + 
    scale_x_discrete(labels=labels) +
    
    theme(legend.position="none") +
    annotate("text", x=1.5, y=ypos, label= paste0("pvalue = ", pvalue), size = 8/.pt )
  
  
  print(g)
  ggsave(g, width = 2, height = 3, dpi = 500, file=paste0(pathout, title,".pdf"))
  dev.off()
  
}

#######################################################
# this function used to make x 2 decimal points
scaleDecimalPoints2 <- function(x)  sprintf("%.2f", x)
scaleDecimalPoints1 <- function(x)  sprintf("%.1f", x)

#################################################
# this is used in PLvsPE_phtenotype.R 
# phenotype PL vs PE without clinic data

singleboxplot_phenotype_noClinic <- function(df, marker, x, y, grouptype,pvalue, pathout, title="", save=TRUE){
  
  df$casetype =  factor(df$casetype, levels=c("Pleural", "Peritoneal"))
  if (title=="") title = paste(marker)
  if (title=="CK") title = "Pan-CK"
  
  num_pl = length(x)
  num_pe = length(y)
  
  maxy = max(df[ paste0(marker, grouptype)], na.rm=TRUE) 
  ylim = maxy + maxy/5
  ypos = maxy - maxy/50
  
  colorder = c("#8f6798", "#fef17c")
  labels = c(paste0("MPM\n(n=",num_pl,")"), paste0("MPeM\n(n=",num_pe,")"))
  
  g = ggplot(df, aes_string(x="casetype", y= paste0(marker, grouptype), fill="casetype")) +
    geom_boxplot( width = 0.7,   #box width
                  lwd=0.4,       #   line width of the box
                  outlier.size = 0.07)  +   #outlier dot size  #outlier.shape = 2, can chage the shape of outlier
    
    # scale_fill_viridis(discrete = TRUE, alpha=0.6) +  set box color

    geom_jitter(color="black", size=0.2, alpha=0.2, width = 0.2) +
    theme_bw()+  # white + grid with box
    # theme_classic()+ # white with left and bottom axis
    
    # panel.border = element_blank())
    theme(axis.text.x = element_text(angle = 0, vjust = 0, hjust=0.5, size=7))+
    theme(axis.text.y = element_text(size=7)) +
    theme(legend.position="none") +
    My_Theme +
    
    ggtitle("") + #paste0(title, "\n","pvalue = ", pvalue)
    theme(plot.title = element_text(hjust = 0.5,size=8, face = "bold")) +
    
    ylab(tools::toTitleCase(substring(grouptype, 2, nchar(grouptype)))) +
    xlab("")+
    ylim (0, ylim) +
    
    scale_fill_manual( values = colorder,  # this is for legned
                       labels = labels) +
   
    scale_x_discrete(labels=labels) +
    
    # annotate("text", x=1.5, y=ypos, label= paste0(), size = 8/.pt ) +
    scale_y_continuous(labels=scaleDecimalPoints2)   #force the y ticker 2 decimal points, so all the plots will be the same width
  
  
  print(g)
  if (save){
    ggsave(g, width = 1.5, height = 1.97, dpi = 500, file=paste0(pathout, title,".pdf"))
    ggsave(g, width = 1.5, height = 1.97, dpi = 500, file=paste0(pathout, title,".png"))
    dev.off()
  }
 
  
}


###################################################
# used to generate the boxplot for phenotype PL vs PE with clinic data
singleboxplot_phenotype <- function(df, marker, x, y, grouptype,type, pvalue, opt, pathout){
  
  df$casetype =  factor(df$casetype, levels=c("Pleural", "Peritoneal"))
  title = paste0(type, grouptype,".",marker,"." , opt)
  num_pl = length(x)
  num_pe = length(y)
  
  maxy = max(df[ paste0(marker, grouptype)], na.rm=TRUE) 
  ylim = maxy + maxy/5
  ypos = maxy + maxy/20
  
  colorder = c("#8f6798", "#fef17c")
  labels = c(paste0("MPM\n(n=",num_pl,")"), paste0("MPeM\n(n=",num_pe,")"))
  
  g = ggplot(df, aes_string(x="casetype", y= paste0(marker, grouptype), fill="casetype")) +
    geom_boxplot( width = 0.7,   #box width
                  lwd=0.4,       #   line width of the box
                  outlier.size = 0.07)  +   #outlier dot size  #outlier.shape = 2, can chage the shape of outlier
    geom_jitter(color="black", size=0.2, alpha=0.2, width = 0.2) +
    
    theme_bw()+
    
    theme(axis.text.x = element_text(angle = 0, vjust = 0, hjust=0.5, size=7))+
    theme(axis.text.y = element_text(size=7)) +
    theme(plot.title = element_text(hjust = 0.5,size=8, face = "bold")) +
    My_Theme +
    
    ggtitle("") + #paste(type,marker, opt)
    ylab(tools::toTitleCase(substring(grouptype, 2, nchar(grouptype)))) +
    xlab("")+
    ylim (0, ylim) +
    
    scale_fill_manual( values = colorder, 
                       labels = labels) +
    
    # scale_fill_manual( values = colorder, 
    #                    labels = c("Pleural" = paste0("Pleural\n(n=",num_pl,")"), "Peritoneal" =  paste0("Peritoneal\n(n=",num_pe,")"))) + 
    scale_x_discrete(labels=labels) +
    
    theme(legend.position="none") +
    # annotate("text", x=1.5, y=ypos, label= paste0("pvalue = ", pvalue), size = 8/.pt )+
    scale_y_continuous(labels=scaleDecimalPoints2)
  
  
  print(g)
  ggsave(g, width = 1.5, height = 1.97, dpi = 500, file=paste0(pathout, title,".pdf"))
  ggsave(g, width = 1.5, height = 1.97, dpi = 500, file=paste0(pathout, title,".png"))
  
  dev.off()
  
}


############################################################################################

singleboxplot_contactScore <- function(df, score, x, y, type, pvalue,  pathout, pair){
  
  title = paste0(pair," ", type)
  num_pl = length(x)
  num_pe = length(y)
  
  maxy = max(df[score], na.rm=TRUE) 
  ylim = maxy + maxy/5
  ypos = maxy + maxy/6
  
  colorder = c("#A6D854", "#8DA0CB")
  
  g = ggplot(df, aes_string(x=type, y= score, fill=type)) +
    # geom_boxplot() +
    # geom_jitter(color="black", size=0.2, alpha=0.5) +
    
    geom_boxplot( width = 0.7,   #box width
                  lwd=0.4,       #   line width of the box
                  outlier.size = 0.07)  +   #outlier dot size  #outlier.shape = 2, can chage the shape of outlier
    
    geom_jitter(color="black", size=0.2, alpha=0.2, width = 0.2) +
    
    theme_bw()+
    
    theme(axis.text.x = element_text(angle = 0, vjust = 0, hjust=0.5, size=7)) +
    theme(axis.text.y = element_text(size=7)) +
    theme(plot.title = element_text(hjust = 0.5,size=6, face = "bold")) +
    My_Theme +
    
    ggtitle(title) +
    ylab("Contact score") +
    xlab("")+
    ylim (0, ylim) +
    
    scale_fill_manual( values = colorder) +   #change legend color, actually is the filled color
  
  # scale_fill_manual( values = colorder, 
  #                   labels = c("Pleural" = paste0("Pleural\n(n=",num_pl,")"), "Peritoneal" =  paste0("Peritoneal\n(n=",num_pe,")"))) + 
  # scale_x_discrete(labels=labels) +   #change ticker label
  
  theme(legend.position="none") +
    # annotate("text", x=1.5, y=ypos, label= paste0("pvalue = ", pvalue), size = 8/.pt )+
    scale_y_continuous(labels=scaleDecimalPoints1)
  
  
  print(g)
  ggsave(g, width = 1.5, height = 2, dpi = 500, file=paste0(pathout, title,".pdf"))
  ggsave(g, width = 1.5, height = 2, dpi = 500, file=paste0(pathout, title,".png"))
  dev.off()
  
}

#########################################################################################################
My_Theme = theme(
  axis.title = element_text(size = 8),
  axis.text = element_text(size = 7))+
  theme(axis.line = element_line(color='black'),  #remove grid
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

####################################################################################################
singleboxplot_contactScore_plVSpe_noClinic <- function(df, score, x, y, pvalue, pathout="", pair="", ypos=0, save=F){

  df$casetype =  factor(df$casetype, levels=c("Pleural", "Peritoneal")) # change the x lable order
  title = paste0(pair)
  num_pl = length(x)
  num_pe = length(y)
  
  maxy = max(df[score], na.rm=TRUE) 
  ylim = maxy + maxy/5
  
  if (ypos == 0) {
    ypos = maxy + maxy/50
  }
  
  
  colorder = c("#8f6798", "#fef17c")
  labels = c(paste0("MPM\n(n=",num_pl,")"), paste0("MPeM\n(n=",num_pe,")"))
  
  g = ggplot(df, aes_string(x="casetype", y= score, fill="casetype")) +
    geom_boxplot( width = 0.5,   #box width
                  lwd=0.4,       #   line width of the box
                  outlier.size = 0.07)  +   #outlier dot size  #outlier.shape = 2, can change the shape of outlier
    
    
    geom_jitter(color="black", size=0.2, alpha=0.1, width = 0.2) +
    
    theme_bw()+
    
    theme(axis.text.x = element_text(angle = 0, vjust = 0, hjust=0.5, size=7))+
    theme(axis.text.y = element_text(size=7))+
    theme(plot.title = element_text(hjust = 0.5,size=8, face = "bold")) +
    My_Theme +
    
    
    ggtitle("") +  #paste0(title , "\n", "pvalue = ", pvalue)
    ylab("Contact score") +
    xlab("")+
    ylim (0, ylim) +
    
    scale_fill_manual( values = colorder, 
                       labels = labels) +    #change legend color, actually is the filled color
    
    # scale_fill_manual( values = colorder, 
    #                   labels = c("Pleural" = paste0("Pleural\n(n=",num_pl,")"), "Peritoneal" =  paste0("Peritoneal\n(n=",num_pe,")"))) + 
    scale_x_discrete(labels=labels) +   #change ticker label
    
    theme(legend.position="none") +
    # annotate("text", x=1.5, y=ypos, label= paste0("pvalue = ", pvalue), size = 8/.pt ) +
    scale_y_continuous(labels=scaleDecimalPoints1)
  
  
  print(g)
  if (save){
    ggsave(g, width = 1.5, height = 2, dpi = 500, file=paste0(pathout, title,".pdf"))
    ggsave(g, width = 1.5, height = 2, dpi = 500, file=paste0(pathout, title,".png"))
    dev.off()
  }
}


####################################################################################
singleboxplot_contactScore_plVSpe <- function(df, score, x, y, type, pvalue, opt, pathout, pair){
  
 
  # opt means one clinic category sep. eg. high, or low
  
  df$casetype =  factor(df$casetype, levels=c("Pleural", "Peritoneal")) # change the x lable order
  title = opt #paste0(pair," ", type," " , opt)
  num_pl = length(x)
  num_pe = length(y)
  
  maxy = max(df[score], na.rm=TRUE) 
  ylim = maxy + maxy/5
  ypos = maxy + maxy/20
  
  colorder = c("#8f6798", "#fef17c")
  labels = c(paste0("MPM\n(n=",num_pl,")"), paste0("MPeM\n(n=",num_pe,")"))
  
  g = ggplot(df, aes_string(x="casetype", y= score, fill="casetype")) +
    
    geom_boxplot( width = 0.7,   #box width
                  lwd=0.4,       #   line width of the box
                  outlier.size = 0.07)  +   #outlier dot size  #outlier.shape = 2, can chage the shape of outlier
    geom_jitter(color="black", size=0.2, alpha=0.2, width = 0.2) +
    
    theme_bw()+
    
    theme(axis.text.x = element_text(angle = 0, vjust = 0, hjust=0.5, size=7))+
    theme(axis.text.y = element_text(size=7))+
    theme(plot.title = element_text(hjust = 0.5,size=8, face = "bold")) +
    theme(axis.text=element_text(size=6)) +
    My_Theme +
    
    ggtitle("") +  #title
    ylab("Contact score") +
    xlab("")+
    ylim (0, ylim) +

    scale_fill_manual( values = colorder, 
                       labels = labels) +    #change legend color, actually is the filled color
    
    # scale_fill_manual( values = colorder, 
    #                   labels = c("Pleural" = paste0("Pleural\n(n=",num_pl,")"), "Peritoneal" =  paste0("Peritoneal\n(n=",num_pe,")"))) + 
    scale_x_discrete(labels=labels) +   #change ticker label
    
    theme(legend.position="none") +
    # annotate("text", x=1.5, y=ypos, label= paste0("pvalue = ", pvalue), size = 8/.pt ) +
    scale_y_continuous(labels=scaleDecimalPoints1)
  
  
  print(g)
  ggsave(g, width = 1.5, height = 1.97, dpi = 500, file=paste0(pathout, pair, "_", title,".pdf"))
  ggsave(g, width = 1.5, height = 1.97, dpi = 500, file=paste0(pathout, pair, "_", title,".png"))
  dev.off()
  
}
####################################################################

generate_pvalues_phenotype_plVSpe <- function(df, type, grouptype, opt, pathout, plot=TRUE, wdiff=FALSE){
  
  markers = c("CD4", "CD8", "CK", "CD20", "FOXP3", "CD68", "CD11c", "CD56")
  ps = list()
  diffs = list()
  for (marker in markers){
    # process high
    x = df[ df["casetype"] == "Pleural", paste0(marker, grouptype)]
    y = df[ df["casetype"] == "Peritoneal", paste0(marker, grouptype)]
    if ( (length(x) <= 3) | (length(y) <= 3)){
      
      print(paste(marker, opt, "less than 3"))
      diffs[[marker]] = NA
      ps[[marker]] = NA

    }else{
      
      if (median(x) > median(y)){
        alter = "greater"
      }else{
        alter = "less"
      }
      
      pvalue = t.test(x, y, alternative = alter )$p.value
      pvalue_f =  formatC(pvalue, format = "e", digits = 2)
      
      diffs[[marker]] = median(x) - median(y)
      ps[[marker]] = pvalue_f
      
      if (plot & ( (pvalue <0.05)  | (marker == "CD8") | (marker == "CK") ) ) {
        singleboxplot_phenotype(df, marker, x, y, grouptype, type, pvalue_f, opt, pathout)
      }
    }
  }
  
  
  if (wdiff){
    ret = list()
    ret[[1]] = ps
    ret[[2]] = diffs
    ret[[3]] = c(sum(df["casetype"] == "Pleural"), sum(df["casetype"] == "Peritoneal"))
    return(ret)
  }
  else{
    return(ps)    
  }
      
}
  

########################################################################################
generate_pvalues_contactScore_plVSpe <- function(df, type, score, opt, pathout, pair, plot=T,  parametric = F){
  # this is PL vs PE
  
    x = df[ df["casetype"] == "Pleural", score] %>% na.omit()
    y = df[ df["casetype"] == "Peritoneal", score]  %>% na.omit()
    if ( (length(x) <= 3) | (length(y) <= 3)){
      pvalue_f = NA
      print(paste(type, opt, "less than 3"))
      diff = NA
    }else{
      if (median(x) > median(y)){
        alter = "greater"
      }else{
        alter = "less"
      }
      diff = median(x) - median(y)
      if (parametric){
        pvalue = t.test(x, y, alternative = alter )$p.value
      }else{
        pvalue = wilcox.test(x, y, alternative = alter )$p.value
      }
      pvalue_f =  formatC(pvalue, format = "e", digits = 2)
      
      # if (plot & pvalue <0.05){
      if (plot){  
      singleboxplot_contactScore_plVSpe(df, score, x, y, type, pvalue_f, opt, pathout, pair)
        
      }
    }
    
    
    ret = list()
    ret[[1]] = pvalue_f
    ret[[2]] = diff
    ret[[3]] = c(sum(df["casetype"] == "Pleural"), sum(df["casetype"] == "Peritoneal"))
    return(ret)
      
    
}




generate_pvalues_contactScore_clinic <- function(df, type, type_values, score,  pathout, pair){
  # this is for clinic category high vs low
  x = df_type[df_type[type] == type_values[1], "pair_score"]
  y = df_type[df_type[type] == type_values[2], "pair_score"]
  
  if ( (length(x) <= 3) | (length(y) <= 3)){
    pvalue_f = NA
    print(paste(type, "less than 3"))
  }else{
    if (median(x) > median(y)){
      alter = "greater"
    }else{
      alter = "less"
    }
    pvalue = t.test(x, y, alternative = alter )$p.value
    pvalue_f =  formatC(pvalue, format = "e", digits = 2)
    
    # if (pvalue <0.05){
      
      singleboxplot_contactScore(df, score, x, y, type, pvalue_f, pathout, pair)
      
    # }
  }
  
  return(pvalue_f)  
  
}
  