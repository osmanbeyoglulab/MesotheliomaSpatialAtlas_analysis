import matplotlib.pyplot as plt
import collections
import pandas as pd
import numpy as np
from utils import PiePlot, GetTypeCompForCol, BarPlot, ridge_plot, violin_plot, bar_plot, PiePlot2, \
    PiePlot_combined, heatmap_plot, bar_plot_series, scatter_plot,bar_plot_plus, bar_plot_plus_withAnno, \
    waterfall_plot, plot_with_anno, PiePlot_combined_byPatient, get_marker_label, get_pvalue,\
        generate_pvalues, generate_barplot_and_pvalues, PiePlot_df1_byPatient,\
            violin_plot_plus, box_plot_plus, box_plot_plus_horizontal
import seaborn as sns
from scipy import stats
import sys
# from statannotations.Annotator import Annotator
import itertools

# df_core = pd.read_csv("../Data_folder/TMA1/panel1_data_cleaned_withCelltypes.csv")
df_core = pd.read_csv("../Data_folder/TMA1/panel1_data_cleaned_withCelltypes_03-29.csv")
df_map = pd.read_csv("../Data_folder/TMA1/mapping.csv")

df_map = df_map.loc[df_map['Classification']=='Malignant', :]
df_map.loc[df_map['subtype']== "epithelial or epithelioid", 'subtype'] = "epithelioid"
# df_core2 = pd.read_csv("../Data_folder/TMA2/panel2_data_cleaned_withMarkerCelltypes.csv")

df_core2 = pd.read_csv("../Data_folder/TMA2/panel2_data_cleaned_withCelltypes_03-29.csv")
df_density2 = pd.read_csv("../Data_folder/TMA2/densities_tma2.csv")
df_map_density2 = pd.merge(df_map, df_density2, how='inner', on="Annotation2")

df = pd.merge(df_core, df_map, how="left", on="Annotation ID") #(781689, 197)
df2 = pd.merge(df_core2, df_map_density2, how="left", on="Annotation2")

df.dropna(subset = ['Institute'], inplace=True) #(576032, 197)  #105 MVB, 315 cores
df2.dropna(subset = ['Institute'], inplace=True) #(606114, 191) #104 MVB 304 cores

# only select the cores of df that exists id df2, df has the same core/MVB with df2
df = df.loc[df['Annotation ID'].isin(df2['Annotation ID']), :] ##104 MVB 304 cores

# completely remove other cells
df.loc[df['phenotype_combined'] == 'other', 'phenotype_combined'] = 'Unidentified'


# df.to_csv("../Data_folder/TMA1/panel1_data_cleaned_withCelltypes_merged.csv")
# df2.to_csv("../Data_folder/TMA1/panel2_data_cleaned_withCelltypes_merged.csv")

df.to_parquet("../Data_folder/TMA1/panel1_data_cleaned_final.parquet.gzip")
df2.to_parquet("../Data_folder/TMA1/panel2_data_cleaned_final.parquet.gzip"     )
#---------------------------------------------------------------------------------------
#############################################################################################
#read in data

df = pd.read_parquet("../Data_folder/TMA1/panel1_data_cleaned_final.parquet.gzip")
df2 = pd.read_parquet("../Data_folder/TMA1/panel2_data_cleaned_final.parquet.gzip")
#----------------------------------------------------------
types = df['phenotype_CD8']
title = "Phenotype Composition of all data"
PiePlot(types, title)

#------------------------
#panel2
types = df2['phenotype_marker']
title = "Panel2 Phenotype Composition of markers"
PiePlot2(types, title)
#--------------------------------
# combine panel1 and panel2
markers = ["CD20", "CD4", "CD8", "FOXP3", "CD68", "CD11c", "CD56", "CK", "Unidentified"]
types1 = df['phenotype_combined']
types2 = df2['phenotype_combined']
title = "Phenotype Composition of all markers"
PiePlot_combined(types1, types2, title, order=markers)



# casetype combine panel1 and panel2
casetype = "Pleural" #"Peritoneal" #" # # #"Peritoneal"# "epithelioid"   #"sarcomatoid" "biphasic"
types1 = df.loc[df['CaseType']== casetype, 'phenotype_combined']
types2 = df2.loc[df2['CaseType']== casetype, 'phenotype_combined']
title = f"{casetype} - Phenotype Composition of all markers"
PiePlot_combined(types1, types2, title, order=markers, angle=0)   

#---------------------------------------------
# subtype combine panel1 and panel2
subtype = "epithelioid" #"biphasic" # #"sarcomatoid"#    #"sarcomatoid" 
types1 = df.loc[df['subtype']== subtype, 'phenotype_combined']
types2 = df2.loc[df2['subtype']== subtype, 'phenotype_combined']
                
title = f"{subtype}: Phenotype Composition of all markers"
PiePlot_combined(types1, types2, title, order=markers, angle=0)        

## subtype+casetype combine panel1 and panel2
subtype ="epithelioid" # "biphasic" #sarcomatoid"  #epithelioid" #"biphasic"  #   , ,]:
casetype = "Pleural" #"Peritoneal" # # "Peritoneal" # #
types1 = df.loc[(df['subtype']== subtype) & (df['CaseType']==casetype), 'phenotype_combined']
types2 = df2.loc[(df2['subtype']== subtype) & (df2['CaseType']==casetype), 'phenotype_combined']
                
title = f"{casetype}-{subtype}: Phenotype Composition of all markers"
PiePlot_combined(types1, types2, title, order=markers, angle=0)    
  
###############################################
# composition combined by patient 
markers = ["CD20", "CD4", "CD8", "FOXP3", "CD68", "CD11c", "CD56", "CK", "Unidentified"]
title = "Phenotype Composition of all markers"
PiePlot_combined_byPatient(df,df2, title, markers, angle=20, save=True)
   
for casetype in ["Pleural", "Peritoneal"]:  
    title = f"{casetype} - Phenotype Composition of all markers"   
    PiePlot_combined_byPatient(df,df2, title, markers, angle=25, casetype = casetype, nolegend=True, save=True)  

for casetype in ["Pleural", "Peritoneal"]:  
    for subtype in ["epithelioid" ,"biphasic" ,"sarcomatoid"]:
        title = f"{casetype}-{subtype} Phenotype Composition of all markers"
        if casetype == "Pleural" and subtype=="sarcomatoid":
            angle = 80
        else:
            angle = 40
        PiePlot_combined_byPatient(df,df2, title, markers, angle=angle, casetype = casetype, subtype=subtype, nolegend=True, save=True)  
    
###############################################
# composition tumor vs stroma by patient 
markers = ["CD20", "CD8", "CD4", "FOXP3", "CD68", "CK", "Unidentified"]

for comp in ["Tumor", "Stroma"]:
    df1 =  df.loc[df['Tissue Category']== comp, :] 
 
    for casetype in ["All", "Pleural", "Peritoneal"]:  
        if casetype == "All":
            title = f"Phenotype Composition for {comp}"
            PiePlot_df1_byPatient(df1 ,title, markers, angle=20)
        else:   
            
            title = f"{casetype} - Phenotype Composition for {comp}"  
            PiePlot_df1_byPatient(df1, title, markers, angle=25, casetype = casetype, nolegend=True)  

# for casetype in ["Pleural", "Peritoneal"]:  
#     for subtype in ["epithelioid" ,"biphasic" ,"sarcomatoid"]:
#         title = f"{casetype}-{subtype}: Phenotype Composition of all markers"
#         if casetype == "Pleural" and subtype=="sarcomatoid":
#             angle = 80
#         else:
            angle = 40
        PiePlot_combined_byPatient(df,df2, title, markers, angle=angle, casetype = casetype, subtype=subtype, nolegend=True)  

#----------------------------
# title = "Phenotype Composition without other"
# df_sub = df.loc[df['phenotype_CD8']!= 'other',:]
# types = df_sub['phenotype_CD8']
# PiePlot(types, title)

title = "Tissue Category composition"
types = df['Tissue Category']
PiePlot(types, title)

title = "CaseType composition"
types = df['CaseType']
PiePlot(types, title)

title = "Institute composition"
types = df['Institute']
PiePlot(types, title)

title = "Subtype composition"
types = df['subtype']
PiePlot(types, title, angle=210)

casetype = 'Peritoneal' #"Pleural"
title = f"{casetype} - subtype composition"
df_sub = df.loc[df['CaseType']==casetype,:]
types = df_sub['subtype']
PiePlot(types, title, angle=210)


#----------------------------------------------------------------------
#panel2
title = "Panel2 Phenotype Composition on Pleural"
df_sub = df2.loc[df2['CaseType']== 'Pleural',:]
types = df_sub['phenotype_marker']
PiePlot(types, title)

title = "Panel2 Phenotype Composition on Peritoneal"
df_sub = df2.loc[df2['CaseType']== 'Peritoneal',:]
types = df_sub['phenotype_marker']
PiePlot(types, title)

#-------------------------------------------
title = "Phenotype Composition on Stroma"
df_sub = df.loc[df['Tissue Category']== 'Stroma',:]
types = df_sub['phenotype_CD8']
PiePlot(types, title)

title = "Phenotype Composition on Tumor"
df_sub = df.loc[df['Tissue Category']== 'Tumor',:]
types = df_sub['phenotype_CD8']
PiePlot(types, title)

#--------------------------------------
#--------------------------------------

title = "Subtype Composition of Pleural"
df_sub = df.loc[df['CaseType']== 'Pleural',:]
types = df_sub['subtype']
PiePlot(types, title)

title = "Subtype Composition of Peritoneal"
df_sub = df.loc[df['CaseType']== 'Peritoneal',:]
types = df_sub['subtype']
PiePlot(types, title)


#--------------------------------------

title = "Phenotype Composition on Pleural Stroma"
df_sub = df.loc[(df['Tissue Category']== 'Stroma') & (df['CaseType']=="Pleural"),:]
types = df_sub['phenotype_CD8']
PiePlot(types, title)

title = "Phenotype Composition on Pleural Tumor"
df_sub = df.loc[(df['Tissue Category']== 'Tumor') & (df['CaseType']=="Pleural"),:]
types = df_sub['phenotype_CD8']
PiePlot(types, title)


#--------------------------------------

title = "Phenotype Composition on Peritoneal Stroma"
df_sub = df.loc[(df['Tissue Category']== 'Stroma') & (df['CaseType']=="Peritoneal"),:]
types = df_sub['phenotype_CD8']
PiePlot(types, title)

title = "Phenotype Composition on Peritoneal Tumor"
df_sub = df.loc[(df['Tissue Category']== 'Tumor') & (df['CaseType']=="Peritoneal"),:]
types = df_sub['phenotype_CD8']
PiePlot(types, title)
###############################################################

col='CaseType'
_, df_count, df_percent = GetTypeCompForCol(df,col,'phenotype_CD8')

title = "Phenotype counts composition by CaseType "
BarPlot(df_count, title, angle=0)

title = "Phenotype percentage composition by CaseType"
BarPlot(df_percent, title, angle=0)

title = "Phenotype percentage composition by CaseType    Other"
BarPlot(df_percent.loc[df_percent.index !='Other',:],title, angle=0)


# exlcude CK and other
col="CaseType"
df_sub = df.loc[~ df['phenotype_CD8'].isin(["CK", "other"]),:]
_, df_count, df_percent = GetTypeCompForCol(df_sub,col,'phenotype_CD8')

title = "phenotype percentage composition(exclude CK and other) by CaseType"
BarPlot(df_percent.loc[df_percent.index !='Other',:],title, angle=0)

#---------------------------------------------
# Panel2
col='CaseType'
df2_sub = df2.loc[df2['phenotype_marker']!="other",:]
_, df_count, df_percent = GetTypeCompForCol(df2_sub,col,'phenotype_marker')

title = "Panel2 Phenotype percentage composition by CaseType without Other"
BarPlot(df_percent.loc[df_percent.index !='Other',:],title, angle=0)

#--------------------------------------------
# Tumor
col='CaseType'
df_sub = df.loc[df['Tissue Category']=="Tumor",:]
_, df_count, df_percent = GetTypeCompForCol(df_sub,col,'phenotype_CD8')

title = "Tumor phenotype counts composition by CaseType "
BarPlot(df_count, title, angle=0)

title = "Tumor phenotype percentage compositionby CaseType"
BarPlot(df_percent, title, angle=0)

title = "Tumor phenotype percentage composition by CaseType without Other"
BarPlot(df_percent.loc[df_percent.index !='Other',:],title, angle=0)


# Tumore exlcude CK and other
col="CaseType"
df_sub = df.loc[(~ df['phenotype_CD8'].isin(["CK", "other"])) & (df['Tissue Category']=="Tumor"),:]
_, df_count, df_percent = GetTypeCompForCol(df_sub,col,'phenotype_CD8')

title = "Tumor phenotype percentage composition(exclude CK and other) by CaseType"
BarPlot(df_percent.loc[df_percent.index !='Other',:],title, angle=0)


# Stroma
col='CaseType'
df_sub = df.loc[df['Tissue Category']=="Stroma",:]
_, df_count, df_percent = GetTypeCompForCol(df_sub,col,'phenotype_CD8')

title = "Stroma phenotype counts composition by CaseType "
BarPlot(df_count, title, angle=0)

title = "Stroma phenotype percentage compositionby CaseType"
BarPlot(df_percent, title, angle=0)

title = "Stroma phenotype percentage composition by CaseType without Other"
BarPlot(df_percent.loc[df_percent.index !='Other',:],title, angle=0)

# Stroma exlcude CK and other
col="CaseType"
df_sub = df.loc[(~ df['phenotype_CD8'].isin(["CK", "other"])) & (df['Tissue Category']=="Stroma"),:]
_, df_count, df_percent = GetTypeCompForCol(df_sub,col,'phenotype_CD8')

title = "Stroma phenotype percentage composition(exclude CK and other) by CaseType"
BarPlot(df_percent.loc[df_percent.index !='Other',:],title, angle=0)

#----------------------------------------------------------
col='CaseType'
_, df_count, df_percent = GetTypeCompForCol(df,col,'subtype')

title = "Subtype counts composition by CaseType "
BarPlot(df_count, title, angle=0, changeColor=True)

title = "Subtype percentage composition by CaseType"
BarPlot(df_percent, title, angle=0, changeColor=True)

title = "Subtype percentage composition by CaseType without Other"
BarPlot(df_percent.loc[df_percent.index !='Other',:],title, angle=0, changeColor=True)

# Tumor
df_sub = df.loc[df['Tissue Category']=="Tumor",:]
_, df_count, df_percent = GetTypeCompForCol(df_sub,col,'subtype')

title = "Tumor subtype counts composition by CaseType "
BarPlot(df_count, title, angle=0, changeColor=True)

title = "Tumor subtype percentage composition by CaseType without Other"
BarPlot(df_percent.loc[df_percent.index !='Other',:],title, angle=0, changeColor=True)

# Stroma
df_sub = df.loc[df['Tissue Category']=="Stroma",:]
_, df_count, df_percent = GetTypeCompForCol(df_sub,col,'subtype')

title = "Stroma subtype counts composition by CaseType "
BarPlot(df_count, title, angle=0, changeColor=True)

title = "Stroma subtype percentage composition by CaseType without Other"
BarPlot(df_percent.loc[df_percent.index !='Other',:],title, angle=0, changeColor=True)


#------------------------------------------------------------------------
#phenotype composition in common subtype between pleural and and peritoneal
col='CaseType'
df_sub = df.loc[df['subtype']=="Papillary", :]
_, df_count, df_percent = GetTypeCompForCol(df_sub, col,'phenotype_CD8')

title = "Papillary phenotype counts composition by CaseType "
BarPlot(df_count, title, angle=0, changeColor=False)

title = "Papillary phenotype percentage composition by CaseType without other"
BarPlot(df_percent.loc[df_percent.index !='Other',:],title, angle=0, changeColor=False)

#------------------------------------------
df_sub = df.loc[df['subtype']=="biphasic", :]
_, df_count, df_percent = GetTypeCompForCol(df_sub, col,'phenotype_CD8')

title = "biphasic phenotype counts composition by CaseType "
BarPlot(df_count, title, angle=0, changeColor=False)

title = "biphasic phenotype percentage composition by CaseType without other"
BarPlot(df_percent.loc[df_percent.index !='Other',:],title, angle=0, changeColor=False)

#------------------------------------------
df_sub = df.loc[df['subtype']=="epithelioid", :]
_, df_count, df_percent = GetTypeCompForCol(df_sub, col,'phenotype_CD8')

title = "epithelioid phenotype counts composition by CaseType "
BarPlot(df_count, title, angle=0, changeColor=False)

title = "epithelioid phenotype percentage composition by CaseType without other"
BarPlot(df_percent.loc[df_percent.index !='Other',:],title, angle=0, changeColor=False)



################################################################################
# Pleural subtype phenotype composition
col = "subtype"
df_sub = df.loc[df['CaseType']=="Pleural", :]
_, df_count, df_percent = GetTypeCompForCol(df_sub, col,'phenotype_CD8')

title = "Pleural phenotype counts composition by subtype "
BarPlot(df_count, title, angle=45, changeColor=False)

title = "Pleural phenotype percentage composition by subtype "
BarPlot(df_percent, title, angle=45, changeColor=False)

# Peritoneal subtype phenotype composition
col = "subtype"
df_sub = df.loc[df['CaseType']=="Peritoneal", :]
_, df_count, df_percent = GetTypeCompForCol(df_sub, col,'phenotype_CD8')

title = "Peritoneal phenotype counts composition by subtype "
BarPlot(df_count, title, angle=45, changeColor=False)

title = "Peritoneal phenotype percentage composition by subtype "
BarPlot(df_percent, title, angle=45, changeColor=False)

#######################################################
# tumore vs stroma
col='Tissue Category'
_, df_count, df_percent = GetTypeCompForCol(df,col,'phenotype_CD8')

title = "Phenotype counts composition by Tissue Category"
BarPlot(df_count, title, angle=0)

title = "Phenotype percentage composition by Tissue Category"
BarPlot(df_percent, title, angle=0)

#-----------------------------------------------------------------
# Peritoneal
col='Tissue Category'
df_sub = df_sub = df.loc[df['CaseType']=="Peritoneal", :]
_, df_count, df_percent = GetTypeCompForCol(df_sub,col,'phenotype_CD8')

title = "Peritoneal - Phenotype counts composition by Tissue Category"
BarPlot(df_count, title, angle=0)

title = "Peritoneal - Phenotype percentage composition by Tissue Category"
BarPlot(df_percent, title, angle=0)


#-----------------------------------------------------------------
# Pleural
col='Tissue Category'
df_sub = df_sub = df.loc[df['CaseType']=="Pleural", :]
_, df_count, df_percent = GetTypeCompForCol(df_sub,col,'phenotype_CD8')

title = "Pleural - Phenotype counts composition by Tissue Category"
BarPlot(df_count, title, angle=0)

title = "Pleural - Phenotype percentage composition by Tissue Category"
BarPlot(df_percent, title, angle=0)
#-----------------------------------------------------------------
'''
set(df_map.loc[df_map['CaseType']=="Pleural", 'subtype'])
Out[122]: 
{'Benign Fibrous ',
 'Fibrocystic',
 'Multicystic',
 'Not specified',
 'Papillary',
 'biphasic',
 'epithelioid',
 'sarcomatoid'}

set(df_map.loc[df_map['CaseType']=="Peritoneal", 'subtype'])
Out[123]: 
{'Desmoplastic',
 'Not Specified',
 'Not specified',
 'Papillary',
 'biphasic',
 'epithelioid'}
'''
################################################################
# density


#-------------------------------------------------



# # ridge_plot(df_plot,  "density_stroma", "phenotype")
# violin_plot('Stoma Density of phenotype_CD8',df_plot, phenotype, density, casetype)

#------------------------------------------------------
# cores = list(set(df['Annotation ID'].values))

# density_stroma = []
# density_tumor = []
# density = []
# phenotype = []
# casetype = []
# for k in cores:
#     area_stroma = df_map.loc[(df_map['Annotation ID']==k) , 'stromaArea'].values[0]
#     area_tumor = df_map.loc[df_map['Annotation ID']==k, 'tumorArea'].values[0]
#     area = area_stroma + area_tumor
#     core_casetype = df_map.loc[(df_map['Annotation ID']==k) , 'CaseType'].values[0]
    
#     type_counts = df.loc[(df['Annotation ID']==k) & (df['CaseType']==core_casetype), 'phenotype_intensity'].value_counts()
#     for tk in type_counts.keys():
#         phenotype.append(tk)    
#         casetype.append(core_casetype)
#         density.append(np.log(type_counts.loc[tk]/area))

        
# # panel2 count from density file
# for k in list(set(df2["Annotation2"])):

#     for marker in ["CD11C", "CD56"]:
        
#         area = df_density2.loc[df_density2["Annotation2"]==k, "Tissue Area (mm2)"].values[0]
        
#         count = 
#         count = df_density2.loc[df_density2["Annotation2"]==k, marker+"+"].values[0]
#         if count< 1e-5:
#             print (k, marker, "count=0")
#             continue
            
#         d = area/count *1000*1000
#         core_casetype = df2.loc[df2["Annotation2"]==k, "CaseType"].values[0]
        
#         phenotype.append(marker)
#         casetype.append(core_casetype)
#         density.append( np.log(d))
        
 #-----------------------------------------------------------------------------------
 # panel1
cores = list(set(df['Annotation ID'].values))
density_stroma = []
density_tumor = []
density = []
phenotype = []
casetype = []
for k in cores:
    area_stroma = df_map.loc[(df_map['Annotation ID']==k) , 'stromaArea'].values[0]
    area_tumor = df_map.loc[df_map['Annotation ID']==k, 'tumorArea'].values[0]
    area = area_stroma + area_tumor
    core_casetype = df_map.loc[(df_map['Annotation ID']==k) , 'CaseType'].values[0]
    
    type_counts = df.loc[(df['Annotation ID']==k), 'phenotype_combined'].value_counts()
    for tk in type_counts.keys():
        phenotype.append(tk)    
        casetype.append(core_casetype)
        density.append( np.log(type_counts.loc[tk]/area + sys.float_info.epsilon))
        # density_tumor.append(np.log(type_counts.loc[tk]/area_tumor + sys.float_info.epsilon))
        

       
# panel2 counts from statistics of df2
for k in list(set(df2["Annotation2"])):

    
    area = df_density2.loc[df_density2["Annotation2"]==k, "Tissue Area (mm2)"].values[0]
    
    core_casetype = df_map.loc[(df_map['Annotation2']==k) , 'CaseType'].values[0]
    
    type_counts = df2.loc[(df2['Annotation2']==k) , 'phenotype_combined'].value_counts()
    
    for tk in  ["CD11c", "CD56", "other"]:
        if tk in type_counts.keys():
            phenotype.append(tk)    
            casetype.append(core_casetype)
            density.append(np.log(type_counts.loc[tk]/area + sys.float_info.epsilon))
    
        

                
        
df_plot = pd.DataFrame({'density':density,
                        'casetype':casetype,
                        # 'density_tumor': density_tumor,
                        'phenotype': phenotype
    })  


df_plot = df_plot.loc[df_plot['casetype']!='Other',:]

df_plot.to_csv("./output/density_of_phenotype.csv")


# df_plot = pd.read_csv("./output/density_of_phenotype.csv", index_col=0)
#patient based plot
df_feature_ori = pd.read_csv("./output/cellContact/core_features_allMarkers_withIntensity_MVB_"+inte_type+".csv")

# violin_plot('Density of phenotype_intensity',df_plot, 'phenotype', 'density', 'casetype')

markers = ["CD20", "CD4", "CD8", "FOXP3", "CD68", "CD11c", "CD56", "CK", "Unidentified"]
grouptype = " density"

cols = [marker+grouptype for marker in markers]
df_feature = df_feature_ori[['MVB']+ cols]
df_feature.volumns = ['MVB' + marker]
df_plot = df_feature.melt(id_vars=['MVB'], var_name='phenotype', value_name='density')

bar_plot_plus_1('Density of phenotype',df_plot, 'phenotype', 'density', 'casetype', markers, hue_order=["Pleural", "Peritoneal"])


import statannot
statannot.add_stat_annotation(
    ax,
    data=df_plot,
    x=phenotype,
    y=density,
    hue=casetype,
    box_pairs=[
        (("Biscoe", "Male"), ("Torgersen", "Female")),
        (("Dream", "Male"), ("Dream", "Female")),
    ],
    test="t-test_ind",
    text_format="star",
    loc="outside",
)

#########################################################################
# common subtypes: epithelioid  and biphasic
subtype =  'epithelioid' #"biphasic" #"sarcomatoid"

# cores = list(set(df.loc[df['subtype']==subtype, 'Annotation ID'].values))

density_stroma = []
density_tumor = []
density = []
phenotype = []
casetype = []

for k in list(set(df.loc[df['subtype']==subtype, 'Annotation ID'].values)): #

    core_subtype = df_map.loc[(df_map['Annotation ID']==k) , 'subtype'].values[0]
    # if core_subtype != subtype:
    #      continue

    area_stroma = df_map.loc[(df_map['Annotation ID']==k) , 'stromaArea'].values[0]
    area_tumor = df_map.loc[df_map['Annotation ID']==k, 'tumorArea'].values[0]
    area = area_stroma + area_tumor
    
    core_casetype = df_map.loc[(df_map['Annotation ID']==k) , 'CaseType'].values[0]
    

    type_counts = df.loc[(df['Annotation ID']==k) , 'phenotype_combined'].value_counts()
    for tk in type_counts.keys():
        phenotype.append(tk)    
        casetype.append(core_casetype)
        density.append( np.log(type_counts.loc[tk]/area + sys.float_info.epsilon))
            


# panel2 counts from statistics of df2
for k in list(set(df2.loc[df2['subtype']==subtype, 'Annotation2'].values)):
   
    core_subtype = df_map.loc[(df_map['Annotation2']==k) , 'subtype'].values[0]
     # if core_subtype != subtype:
     #     continue
    area = df_density2.loc[df_density2["Annotation2"]==k, "Tissue Area (mm2)"].values[0]
    
    core_casetype = df_map.loc[(df_map['Annotation2']==k) , 'CaseType'].values[0]
    
    type_counts = df2.loc[(df2['Annotation2']==k), 'phenotype_combined'].value_counts()
    
    for tk in  ["CD11c", "CD56", "other"]:
        if tk in type_counts.keys():
            phenotype.append(tk)    
            casetype.append(core_casetype)
            density.append(np.log(type_counts.loc[tk]/area + sys.float_info.epsilon))
    
        

                
        
df_plot = pd.DataFrame({'density':density,
                        'casetype':casetype,
                        # 'density_tumor': density_tumor,
                        'phenotype': phenotype
    })  

# df_plot.to_csv("../Data_folder/df_log_density_plot.csv")

# ridge_plot(df_plot,  "density_stroma", "phenotype")

df_plot = df_plot.loc[df_plot['casetype']!='Other',:]
# violin_plot('Density of phenotype_intensity',df_plot, 'phenotype', 'density', 'casetype')

markers = ["CD20", "CD4", "CD8", "FOXP3", "CD68", "CD11c", "CD56", "CK", "other"]

bar_plot_plus(f'{subtype} - Density of phenotype',df_plot, 'phenotype', 'density', 'casetype', markers, hue_order=["Pleural", "Peritoneal"]) #


######################################################################################################
# ----------------------------------------------------------------------------------------------------
## subtypes

mycasetype =  "Pleural"  # "Peritoneal" # "Pleural"

density_stroma = []
density_tumor = []
density = []
phenotype = []
subtype = []
casetype = []

subtypes_set = ["epithelioid", "biphasic", "sarcomatoid"]
cores = list(set(df.loc[(df['subtype'].isin(subtypes_set) ) & (df['CaseType']==mycasetype), 'Annotation ID'].values))
# cores = list(set(df.loc[ df['subtype'].isin(subtypes_set), 'Annotation ID'].values))

for k in cores: #

    core_casetype = df_map.loc[(df_map['Annotation ID']==k),'CaseType'].values[0]
    core_subtype = df_map.loc[(df_map['Annotation ID']==k) , 'subtype'].values[0]


    area_stroma = df_map.loc[(df_map['Annotation ID']==k) , 'stromaArea'].values[0]
    area_tumor = df_map.loc[df_map['Annotation ID']==k, 'tumorArea'].values[0]
    area = area_stroma + area_tumor
    
    # core_casetype = df_map.loc[(df_map['Annotation ID']==k) , 'CaseType'].values[0]
    

    type_counts = df.loc[(df['Annotation ID']==k) , 'phenotype_intensity'].value_counts()
    for tk in type_counts.keys():
        casetype.append(core_casetype)
        phenotype.append(tk)    
        subtype.append(core_subtype)
        density.append( np.log(type_counts.loc[tk]/area + sys.float_info.epsilon))
            


# panel2 counts from statistics of df2
cores = list(set( df2.loc[ (df2['subtype'].isin(subtypes_set) ) & (df2['CaseType']==mycasetype),  'Annotation2'].values))
# cores = list(set(df2.loc[ df2['subtype'].isin(subtypes_set),  'Annotation2'].values ))

for k in cores:
   
    core_subtype = df_map.loc[(df_map['Annotation2']==k) , 'subtype'].values[0]
    
    core_casetype = df_map.loc[(df_map['Annotation2']==k),'CaseType'].values[0]
   
    area = df_density2.loc[df_density2["Annotation2"]==k, "Tissue Area (mm2)"].values[0]
    
    # core_casetype = df_map.loc[(df_map['Annotation2']==k) , 'CaseType'].values[0]
    
    type_counts = df2.loc[(df2['Annotation2']==k), 'phenotype_intensity'].value_counts()
    
    for tk in  ["CD11c", "CD56", "other"]:
        if tk in type_counts.keys():
            casetype.append(core_casetype)
            phenotype.append(tk)    
            subtype.append(core_subtype)
            density.append(np.log(type_counts.loc[tk]/area + sys.float_info.epsilon))
    
        

                
        
df_plot = pd.DataFrame({'density':density,
                        'subtype':subtype,
                        'casetype':casetype,
                        # 'density_tumor': density_tumor,
                        'phenotype': phenotype
    })  

df_plot.to_csv("../Data_folder/df_log_density_composition_by_subtype.csv")

# ridge_plot(df_plot,  "density_stroma", "phenotype")

# violin_plot('Density of phenotype_intensity',df_plot, 'phenotype', 'density', 'casetype')

markers = ["CD20", "CD4", "CD8", "FOXP3", "CD68", "CD11c", "CD56", "CK", "other"]

#---------------------------------------------------------------
mycasetype =  "Pleural" 

df_plot = pd.read_csv("../Data_folder/df_log_density_composition_by_subtype.csv", index_col=0)
df_plot = df_plot.loc[df_plot['casetype'] == mycasetype, :]

# pvalues
####################################################################################################
# pvalue  3subtypes Kruskal-Wallis test
pvalues_dict= dict()
for marker in markers:
    x = df_plot.loc[(df_plot['subtype']=="epithelioid") & (df_plot['phenotype']==marker),'density'].values
    y = df_plot.loc[ (df_plot['subtype']=="biphasic") & (df_plot['phenotype']==marker),'density'].values
    z = df_plot.loc[ (df_plot['subtype']=="sarcomatoid") & (df_plot['phenotype']==marker),'density'].values
    pvalues_dict[marker] = round(stats.kruskal(x,y,z).pvalue,5)
 
#plot

title = f'{mycasetype} - Density of phenotype'
data = df_plot
x= 'phenotype'
y='density'
hue = 'subtype'
hue_order=["epithelioid", "biphasic", "sarcomatoid"]
colors = [ "#00E5EE","#FF3E96","#F0F8FF"]
bbox_to_anchor=(1.5, .8)
pvalues = list(pvalues_dict.values())

if mycasetype == "Pleural":
    bar_plot_plus_withAnno(title,df_plot, 'phenotype', 'density', 'subtype', markers, \
                  hue_order=["epithelioid", "biphasic", "sarcomatoid"], \
                      colors = [ "#00E5EE","#FF3E96","#F0F8FF"], \
                          pvalues = list(pvalues_dict.values()), \
                              bbox_to_anchor=(1.5, .8))  
else:  
    bar_plot_plus(title,df_plot, 'phenotype', 'density', 'subtype', markers, hue_order=["epithelioid", "biphasic" ], colors = [ "#00E5EE","#FF3E96"], bbox_to_anchor=(1.5, .8)) #



#####################################################################################################
#pvalue 2 subtypes wilconxon rank sum test

# import seaborn as sns
# sns.histplot(df_plot.loc[df_plot['casetype']=="Peritoneal",'density'].values, bins=100)




#####################################################################################################
#pvalue 2 casetype wilconxon rank sum test
alternatives =pd.Series( ["less", "less","less","greater"])
alternatives =pd.Series( ["greater","less", "greater","greater","greater","greater","less","greater","less"])
alternatives =pd.Series( ["less","greater", "less","greater","greater","less","less","greater","less"])
alternatives.index=markers

                
pvalues = dict()
for marker in markers:
    x = df_plot.loc[(df_plot['casetype']=="Pleural") & (df_plot['phenotype']==marker),'density'].values
    y = df_plot.loc[ (df_plot['casetype']=="Peritoneal") & (df_plot['phenotype']==marker),'density'].values
    pvalues[marker] = round(stats.ranksums(x, y, alternative=alternatives[marker]).pvalue,5)
    
    
 ####################################################################################################
 # pvalue  2subtypes Kruskal-Wallis test
pvalues = dict()
for marker in markers:
    x = df_plot.loc[(df_plot['subtype']=="epithelioid") & (df_plot['phenotype']==marker),'density'].values
    y = df_plot.loc[ (df_plot['subtype']=="not epithelioid and epithelioid") & (df_plot['phenotype']==marker),'density'].values
    
    pvalues[marker] = round(stats.ranksums(x, y, alternative=alternatives[marker]).pvalue,5)
     

  
#---------------------------------------------------------------------------
# combine biphasic and sarcomatoid for Pleural-Density of phenotype plots
casetype =  "Pleural"  # "Peritoneal" # "Pleural"
# cores = list(set(df.loc[df['subtype']==subtype, 'Annotation ID'].values))

density_stroma = []
density_tumor = []
density = []
phenotype = []
subtype = []

subtypes_set = ["epithelioid", "biphasic", "sarcomatoid"]
cores = list(set(df.loc[(df['subtype'].isin(subtypes_set) ) & (df['CaseType']==casetype), 'Annotation ID'].values))

for k in cores: #

    core_subtype = df_map.loc[(df_map['Annotation ID']==k) , 'subtype'].values[0]
    if (core_subtype == "biphasic") | (core_subtype == "sarcomatoid") :
        core_subtype = "not epithelioid and epithelioid"

    area_stroma = df_map.loc[(df_map['Annotation ID']==k) , 'stromaArea'].values[0]
    area_tumor = df_map.loc[df_map['Annotation ID']==k, 'tumorArea'].values[0]
    area = area_stroma + area_tumor
    
    # core_casetype = df_map.loc[(df_map['Annotation ID']==k) , 'CaseType'].values[0]
    

    type_counts = df.loc[(df['Annotation ID']==k) , 'phenotype_combined'].value_counts()
    for tk in type_counts.keys():
        phenotype.append(tk)    
        subtype.append(core_subtype)
        density.append( np.log(type_counts.loc[tk]/area + sys.float_info.epsilon))
            


# panel2 counts from statistics of df2
for k in list(set(df2.loc[(df2['subtype'].isin(subtypes_set) ) & (df2['CaseType']==casetype),  'Annotation2'].values)):
   
    core_subtype = df_map.loc[(df_map['Annotation2']==k) , 'subtype'].values[0]
    if (core_subtype == "biphasic") | (core_subtype == "sarcomatoid") :
        core_subtype = "not epithelioid and epithelioid"

    
    area = df_density2.loc[df_density2["Annotation2"]==k, "Tissue Area (mm2)"].values[0]
    
    # core_casetype = df_map.loc[(df_map['Annotation2']==k) , 'CaseType'].values[0]
    
    type_counts = df2.loc[(df2['Annotation2']==k), 'phenotype_combined'].value_counts()
    
    for tk in  ["CD11c", "CD56", "other"]:
        if tk in type_counts.keys():
            phenotype.append(tk)    
            subtype.append(core_subtype)
            density.append(np.log(type_counts.loc[tk]/area + sys.float_info.epsilon))
    
        

                
        
df_plot = pd.DataFrame({'density':density,
                        'subtype':subtype,
                        # 'density_tumor': density_tumor,
                        'phenotype': phenotype
    })  

# df_plot.to_csv("../Data_folder/df_log_density_plot.csv")

# ridge_plot(df_plot,  "density_stroma", "phenotype")

# violin_plot('Density of phenotype_intensity',df_plot, 'phenotype', 'density', 'casetype')

markers = ["CD20", "CD4", "CD8", "FOXP3", "CD68", "CD11c", "CD56", "CK", "other"]

bar_plot_plus(f'{casetype} - Density of phenotype',df_plot, 'phenotype', 'density', 'subtype', markers, hue_order=["epithelioid", "not epithelioid and epithelioid"]) #"sarcomatoid" 





    
######################################################################################################   
######################################################################################################
# [0502]
'''
df2['Classification'].value_counts()/len(df2['Classification'])
Out[1100]: 
Malignant        567736
Not specified     19631
Benign            18747
----------------------------------------
Malignant        0.936682
Not specified    0.032388
Benign           0.030930


df2['Phenotype-LAG3'].value_counts()
Out[1110]: 
LAG3-    589863
LAG3+     15411                 LAG3 has not much assignments
Name: Phenotype-LAG3, dtype: int64

df2['Phenotype-BAP1'].value_counts()
Out[1111]: 
BAP1-    487314
BAP1+    117960
Name: Phenotype-BAP1, dtype: int64

df2['Phenotype-NF2'].value_counts()
Out[1112]: 
NF2-    411417
NF2+    193857
Name: Phenotype-NF2, dtype: int64

df2['Phenotype-MTAP'].value_counts()
Out[1113]: 
MTAP-    414847
MTAP+    190427
Name: Phenotype-MTAP, dtype: int64

'''
#-----------------------------------------------------------------------
# prepare data
# filter out the cores with less tumor cells

gb = df.loc[df['Classification']=='Malignant', :].groupby('Annotation ID')

tumor_percent = []
stroma_percent = []
for k, gp in gb:
    valuecounts = gp['Tissue Category'].value_counts()/len(gp['Tissue Category'])
    if 'Tumor' in valuecounts.index:
        tumor_percent.append(valuecounts['Tumor'])
        
    if 'Stroma' in valuecounts.index:
        stroma_percent.append(valuecounts['Stroma'])
        
sns.histplot(tumor_percent, bins=100).set(title='Tumor cell percentage in the core')
sns.histplot(stroma_percent, bins=100).set(title='Stroma percentage in the corev')

qt30 = pd.Series(tumor_percent).quantile(0.3)

# select cores based on qt30(tumor percentage at 0.3 quantile) 
# decide to use a absolute value to cut

# check each subtype core tumor cell percentage
# subtype = 'sarcomatoid' # 'epithelioid'#'biphasic' # # 'sarcomatoid'
# gb = df.loc[(df['Classification']=='Malignant') & (df['subtype'] == subtype), :].groupby('Annotation ID')
# tumor_percent = []
# stroma_percent = []
# for k, gp in gb:
#     valuecounts = gp['Tissue Category'].value_counts()/len(gp['Tissue Category'])
#     if 'Tumor' in valuecounts.index:
#         tumor_percent.append(valuecounts['Tumor'])
        
#     if 'Stroma' in valuecounts.index:
#         stroma_percent.append(valuecounts['Stroma'])
        
# sns.histplot(tumor_percent, bins=50).set(title=f'{subtype} cores tumor percentage')
# plt.show()
# sns.histplot(stroma_percent, bins=50).set(title=f'{subtype} cores stroma percentage')



     
cores_sel = []  #149
for k, gp in gb:
    valuecounts = gp['Tissue Category'].value_counts()/len(gp['Tissue Category'])
    if ('Tumor' in valuecounts.index)  :
        if (valuecounts['Tumor'] > 0.25) :
            cores_sel.append((k))
     
cores_sel = []  #149
for k, gp in gb:
    valuecounts = gp['Tissue Category'].value_counts()/len(gp['Tissue Category'])
    valuecounts2 = gp['phenotype_combined'].value_counts()/len(gp['phenotype_combined'])
    if ('Tumor' in valuecounts.index) and ('CK' in valuecounts2.index) :
        if (valuecounts['Tumor'] > 0.25) and (valuecounts2['CK'] > 0.25):
            cores_sel.append((k))


# only consider selected cores (142) in which tumor percentage above 0.3 quantile                          # only consider malignant cores 315
#calculate the density of "BAP1", "NF2", "MTAP", "LAG3"
 # panel1

#-----------------------------------------------------------------------------------------
# start correalation analysis


# markers1 = ["CD20", "CD4", "CD8", "FOXP3", "CD68",  "CK"]
# # genes= ["BAP1", "NF2", "MTAP", "LAG3"] #, "LAG3"
# markers2 = ["CD11c", "CD56", "BAP1", "NF2", "MTAP"]
# intensityCol = {}
# for marker in markers1:
#     if marker == "FOXP3":
#         intensityCol[marker] = "Entire Cell FOXP3 (Opal 570) Mean (Normalized Counts, Total Weighting)"
#     if marker == "CD4":
#         intensityCol[marker] = "Entire Cell CD4 (Opal 690) Mean (Normalized Counts, Total Weighting)"
#     if marker == "CD8":
#         intensityCol[marker] = "Entire Cell CD8 (Opal 480) Mean (Normalized Counts, Total Weighting)"
#     if marker == "CD20":
#         intensityCol[marker] =  "Entire Cell CD20 (Opal 620) Mean (Normalized Counts, Total Weighting)"
#     if marker == "CK":
#         intensityCol[marker] = "Entire Cell PanCK (Opal 780) Mean (Normalized Counts, Total Weighting)"
#     if marker == "CD68":
#         intensityCol[marker] = "Entire Cell CD68 (Opal 520) Mean (Normalized Counts, Total Weighting)"
# for marker in markers2:
#     if marker == "CD11c":
#         intensityCol[marker] = "Entire Cell CD11c (Opal 480) Mean (Normalized Counts, Total Weighting)"
#     if marker == "CD56":
#         intensityCol[marker] = "Entire Cell CD56 (Opal 520) Mean (Normalized Counts, Total Weighting)"
#     if marker == "NF2":
#         intensityCol[marker] = "Entire Cell NF2 (Opal 780) Mean (Normalized Counts, Total Weighting)"
#     if marker == "BAP1":
#         intensityCol[marker] = "Entire Cell BAP1 (Opal 690) Mean (Normalized Counts, Total Weighting)"
#     if marker == "MTAP":
#         intensityCol[marker] = "Entire Cell MTAP (Opal 620) Mean (Normalized Counts, Total Weighting)"
#     if marker == "LAG3":
#         intensityCol[marker] = "Entire Cell LAG3 (Opal 570) Mean (Normalized Counts, Total Weighting)"
        
# density_stroma = []
# density_tumor = []
# density = []
# phenotype = []
# casetype = []
# markercore = []

# percent = []
# intensity = []
# mixtype = []

# # for k in list(set(df.loc[df['Classification']=='Malignant', 'Annotation ID'].values)):
# # for k in list(set(df.loc[df['Annotation ID'].isin(cores_sel), 'Annotation ID'].values)):
# #just select malignant cores, and do the further core selection at the later stage
# for k in list(set(df.loc[df['Classification']=='Malignant', 'Annotation ID'].values)):
#     area_stroma = df_map.loc[(df_map['Annotation ID']==k) , 'stromaArea'].values[0]
#     area_tumor = df_map.loc[df_map['Annotation ID']==k, 'tumorArea'].values[0]
#     area = area_stroma + area_tumor
#     core_casetype = df_map.loc[(df_map['Annotation ID']==k) , 'CaseType'].values[0]
#     # marker counts
#     type_counts = df.loc[(df['Annotation ID']==k), 'phenotype_combined'].value_counts()
#     for tk in type_counts.keys(): # other cause pivot issue with duplicated cunnts in two panels
#         if tk == "Unidentified":
#             continue
#         markercore.append(k)
#         phenotype.append(tk)    
#         casetype.append(core_casetype)
#         density.append( np.log(type_counts.loc[tk]/area + sys.float_info.epsilon))
        
#         intensity.append(df.loc[(df['Annotation ID']==k), intensityCol[tk]].mean())
#         percent.append(type_counts.loc[tk]/type_counts.sum())
#         mixtype.append(np.log(type_counts.loc[tk]/area + sys.float_info.epsilon))
    

 
# # panel2 counts from statistics of df2
# gene_name = []
# gene_density=[]
# gene_casetype = []
# gene_core = []

# gene_intensity = []
# gene_percent = []
# gene_mixtype = []

# genes = ["BAP1", "NF2", "MTAP"] #, "LAG3"]
# for k2 in list(set(df2.loc[df2['Classification']=='Malignant',"Annotation2"].values)):
# # for k2 in list(set(df2.loc[df2['Annotation ID'].isin(cores_sel),"Annotation2"].values)):    
#     area = df_density2.loc[df_density2["Annotation2"]==k2, "Tissue Area (mm2)"].values[0]
    
#     core_casetype = df_map.loc[(df_map['Annotation2']==k2) , 'CaseType'].values[0]
    
#     type_counts = df2.loc[(df2['Annotation2']==k2) , 'phenotype_combined'].value_counts()
    
#     #convert k to k1
#     k1 = df_map.loc[df_map["Annotation2"]==k2, "Annotation ID"].values[0]
#     for tk in  ["CD11c", "CD56"]:
#         if tk in type_counts.keys():
#             phenotype.append(tk)    
#             casetype.append(core_casetype)
#             density.append(np.log(type_counts.loc[tk]/area + sys.float_info.epsilon))
#             markercore.append(k1)
    
#             intensity.append(df2.loc[(df2['Annotation2']==k2),intensityCol[tk]].mean())
#             percent.append(type_counts.loc[tk]/type_counts.sum())
#             mixtype.append(np.log(type_counts.loc[tk]/area + sys.float_info.epsilon))
        
#     # gene counts
    
#     for gene in genes:
        
#         gene_count = sum( (df2['Annotation2']==k2) & (df2[f"Phenotype-{gene}"]==f"{gene}+") )
#         if gene_count == 0:
#             print (f"{gene} in {k2} of {core_casetype} has 0 count")
#             continue
#         gene_density.append(np.log(gene_count/area + sys.float_info.epsilon))
#         gene_casetype.append(core_casetype)
#         gene_name.append(gene)  
#         gene_core.append(k1)
        
#         gene_intensity.append(df2.loc[(df2['Annotation2']==k2), intensityCol[gene]].mean())
#         gene_percent.append(gene_count/type_counts.sum())
#         gene_mixtype.append(df2.loc[(df2['Annotation2']==k2), intensityCol[gene]].mean())
        
# # density plot of other genes        
# df_plot = pd.DataFrame({'density':gene_density,
#                         'casetype':gene_casetype,
#                         'phenotype': gene_name,
#                         'intensity': gene_intensity,
#                         'percent': gene_percent,
#                         'mixtype': gene_mixtype
#     })  

# df_plot = df_plot.loc[df_plot['casetype']!='Other',:]
# # violin_plot('Density of phenotype_intensity',df_plot, 'phenotype', 'density', 'casetype')

# genes= ["BAP1", "NF2", "MTAP"]#, "LAG3"]

# bar_plot_plus('Density of other genes',df_plot, 'phenotype', 'density', 'casetype', genes, hue_order=["Pleural", "Peritoneal"])



# # df_plot.to_csv("../figure/Density of other genes.csv")
# #-----------------------------------------------------------------------------------------------------
# #df markers
# df_marker = pd.DataFrame({'density':density,
#                         'casetype':casetype,
#                         'phenotype': phenotype,
#                         "core":markercore,
#                         'intensity':intensity,
#                         'percent':percent,
#                         'mixtype':mixtype})

                          
# df_gene = pd.DataFrame({'density':gene_density,
#                         'casetype':gene_casetype,
#                         'phenotype': gene_name,
#                         'core':gene_core,
#                         'intensity':gene_intensity,
#                         'percent': gene_percent,
#                         'mixtype': gene_mixtype})

# studyTypes = ["Density",  "Intensity", "Percent",'mixtype']
# for idx in [1]:
#     studyType = studyTypes[idx]
#     #density
#     if studyType == studyTypes[0]:
#         df_gene_pleural = df_gene.loc[df_gene['casetype']=="Pleural",:].pivot(index='core', columns='phenotype')['density'].reset_index()
#         df_gene_peritoneal = df_gene.loc[df_gene['casetype']=="Peritoneal",:].pivot(index='core', columns='phenotype')['density'].reset_index()
        
#         df_marker_pleural = df_marker.loc[df_marker['casetype']=="Pleural",:].pivot(index='core', columns='phenotype')['density'].reset_index()
#         df_marker_peritoneal = df_marker.loc[df_marker['casetype']=="Peritoneal",:].pivot(index='core', columns='phenotype')['density'].reset_index()
    
#         #---------------
        
#         df_gene_pivoted = df_gene.pivot(index=['core','casetype'], columns='phenotype', values='density').reset_index()
#         df_marker_pivoted = df_marker.pivot(index=['core','casetype'], columns='phenotype', values='density').reset_index()
             
#     # intensity
#     if studyType == studyTypes[1]:
#         df_gene_pleural = df_gene.loc[df_gene['casetype']=="Pleural",:].pivot(index='core', columns='phenotype')['intensity'].reset_index()
#         df_gene_peritoneal = df_gene.loc[df_gene['casetype']=="Peritoneal",:].pivot(index='core', columns='phenotype')['intensity'].reset_index()
        
#         df_marker_pleural = df_marker.loc[df_marker['casetype']=="Pleural",:].pivot(index='core', columns='phenotype')['intensity'].reset_index()
#         df_marker_peritoneal = df_marker.loc[df_marker['casetype']=="Peritoneal",:].pivot(index='core', columns='phenotype')['intensity'].reset_index()
    
#         df_gene_pivoted = df_gene.pivot(index=['core','casetype'], columns='phenotype', values='intensity').reset_index()
#         df_marker_pivoted = df_marker.pivot(index=['core','casetype'], columns='phenotype', values='intensity').reset_index()
          
#     # percent
#     if studyType == studyTypes[2]:
#         df_gene_pleural = df_gene.loc[df_gene['casetype']=="Pleural",:].pivot(index='core', columns='phenotype')['percent'].reset_index()
#         df_gene_peritoneal = df_gene.loc[df_gene['casetype']=="Peritoneal",:].pivot(index='core', columns='phenotype')['percent'].reset_index()
        
#         df_marker_pleural = df_marker.loc[df_marker['casetype']=="Pleural",:].pivot(index='core', columns='phenotype')['percent'].reset_index()
#         df_marker_peritoneal = df_marker.loc[df_marker['casetype']=="Peritoneal",:].pivot(index='core', columns='phenotype')['percent'].reset_index()
    
#         df_gene_pivoted = df_gene.pivot(index=['core','casetype'], columns='phenotype', values='percent').reset_index()
#         df_marker_pivoted = df_marker.pivot(index=['core','casetype'], columns='phenotype', values='percent').reset_index()
          
    
#     # mixtype
#     if studyType == studyTypes[3]:
#         df_gene_pleural = df_gene.loc[df_gene['casetype']=="Pleural",:].pivot(index='core', columns='phenotype')['mixtype'].reset_index()
#         df_gene_peritoneal = df_gene.loc[df_gene['casetype']=="Peritoneal",:].pivot(index='core', columns='phenotype')['mixtype'].reset_index()
        
#         df_marker_pleural = df_marker.loc[df_marker['casetype']=="Pleural",:].pivot(index='core', columns='phenotype')['mixtype'].reset_index()
#         df_marker_peritoneal = df_marker.loc[df_marker['casetype']=="Peritoneal",:].pivot(index='core', columns='phenotype')['mixtype'].reset_index()
    
#     # start-- this part used in comparition between BAP1-NF2-MTAP- and BAP1-NF2+MTAP+
#     df_merged = pd.merge(df_gene_pivoted,df_marker_pivoted, how="inner", on=["core", "casetype"])
#     df_merged.set_index('core', inplace=True) 
#     # end --

    
#     markers = ["CD20", "CD4", "CD8", "FOXP3", "CD68", "CD11c", "CD56", "CK"]
#     # genes= ["BAP1", "NF2", "MTAP", "LAG3"] #, "LAG3"
#     genes= ["BAP1", "NF2", "MTAP"]
     
#     df_pleural = pd.merge(df_gene_pleural,df_marker_pleural, how="inner", on="core")
#     df_pleural.set_index('core', inplace=True)
#     corr_pleural = df_pleural.corr(method="spearman")
#     # corr_pleural = df_pleural.dropna(axis=0).corr()
#     corr_pleural = corr_pleural.loc[genes, markers]
#     corr_pleural.to_csv("../figure/corr pleural between marker and genes.csv")
    
#     heatmap_plot(corr_pleural, studyType + " Spearman Corr in Pleural")
#     plt.show()
    
    
#     df_peritoneal = pd.merge(df_gene_peritoneal,df_marker_peritoneal, how="inner", on="core")
#     df_peritoneal.set_index('core', inplace=True)
#     corr_peritoneal = df_peritoneal.corr(method="spearman")
#     # corr_peritoneal = df_peritoneal.dropna(axis=0).corr()
#     corr_peritoneal = corr_peritoneal.loc[genes, markers]
#     corr_peritoneal.to_csv("../figure/corr peritoneal between marker and genes.csv")
    
#     heatmap_plot(corr_peritoneal, studyType + " Spearman Corr in Peritoneal")
#     plt.show()
    
    
# # --------------------------------------------------------------------------- 
# # compare between combined groups
# ## -------------------------------------------
# # compare in total

# BAP1_q50 = df_merged['BAP1'].quantile(0.5)
# NF2_q50 = df_merged['NF2'].quantile(0.5)
# MTAP_q50 = df_merged['MTAP'].quantile(0.5)

# df_merged['BAP1status'] = np.where(df_merged['BAP1'] > BAP1_q50, 'BAP1+', 'BAP1-')
# df_merged['NF2status'] = np.where(df_merged['NF2'] > NF2_q50, 'NF2+', 'NF2-')
# df_merged['MTAPstatus'] = np.where(df_merged['MTAP'] > MTAP_q50, 'MTAP+', 'MTAP-')
# df_merged['status'] = df_merged['BAP1status'] + df_merged['NF2status'] + df_merged['MTAPstatus']
# '''
# df_merged['status'].value_counts()
# BAP1-NF2-MTAP-    28
# BAP1+NF2+MTAP+    24
# BAP1-NF2+MTAP+    17
# BAP1+NF2+MTAP-    17
# BAP1+NF2-MTAP+    16
# BAP1-NF2-MTAP+    15
# BAP1-NF2+MTAP-    14
# BAP1+NF2-MTAP-    14
# '''
# # total
# df_sub = df_merged.loc[(df_merged['status'] == 'BAP1-NF2-MTAP-') | (df_merged['status'] == 'BAP1-NF2+MTAP+'), ['CD11c', 'CD20', 'CD4', 'CD56', 'CD68', 'CD8',
#        'CK', 'FOXP3', 'status']]

# df_sub = df_merged[['CD11c', 'CD20', 'CD4', 'CD56', 'CD68', 'CD8',
#        'CK', 'FOXP3', 'status']]

# df_sub_plot = df_sub.melt(id_vars = 'status', var_name='phenotype', value_name=studyType)

# bar_plot(studyType+' comparison ',df_sub_plot, 'phenotype', studyType, 'status')
# plt.show()
# # Pleural
# df_sub_pleural = df_merged.loc[(df_merged['casetype']=="Pleural") & ((df_merged['status'] == 'BAP1-NF2-MTAP-') | (df_merged['status'] == 'BAP1-NF2+MTAP+')) , ['CD11c', 'CD20', 'CD4', 'CD56', 'CD68', 'CD8',
#        'CK', 'FOXP3', 'status']]

# df_sub_pleural = df_merged.loc[(df_merged['casetype']=="Pleural")  , ['CD11c', 'CD20', 'CD4', 'CD56', 'CD68', 'CD8',
#        'CK', 'FOXP3', 'status']]

# df_sub_plot = df_sub_pleural.melt(id_vars = 'status', var_name='phenotype', value_name=studyType)

# bar_plot(studyType+' comparison in Pleural ',df_sub_plot, 'phenotype', studyType, 'status')
# plt.show()
# # Perttoneal
# df_sub_peritoneal = df_merged.loc[(df_merged['casetype']=="Peritoneal") & ((df_merged['status'] == 'BAP1-NF2-MTAP-') | (df_merged['status'] == 'BAP1-NF2+MTAP+')) , ['CD11c', 'CD20', 'CD4', 'CD56', 'CD68', 'CD8',
#        'CK', 'FOXP3', 'status']]

# df_sub_peritoneal = df_merged.loc[(df_merged['casetype']=="Peritoneal") , ['CD11c', 'CD20', 'CD4', 'CD56', 'CD68', 'CD8',
#        'CK', 'FOXP3', 'status']]

# df_sub_plot = df_sub_peritoneal.melt(id_vars = 'status', var_name='phenotype', value_name=studyType)

# bar_plot(studyType+' comparison in Peritoneal ',df_sub_plot, 'phenotype', studyType, 'status')
# plt.show()
# '''
# df_sub_peritoneal['status'].value_counts()
# Out[373]: 
# BAP1-NF2-MTAP-    7
# BAP1-NF2+MTAP+    7
# Name: status, dtype: int64

# df_sub_pleural['status'].value_counts()
# Out[374]: 
# BAP1-NF2-MTAP-    21
# BAP1-NF2+MTAP+    10
# Name: status, dtype: int64
# '''
# #-----------------------------------------
# # anther way of doing correlation
# markers = ["CD20", "CD4", "CD8", "FOXP3", "CD68", "CD11c", "CD56", "CK"] 
# markers_density = [item + " density" for item in markers]
# genes= ["BAP1", "NF2", "MTAP"]
# genes_density = [item + " density" for item in genes]

# df_core_feature = pd.read_csv("./output/cellContact/core_features_allMarkers.csv")
# cond1 = df_core_feature['tumor density'] > df_core_feature['tumor density'].quantile(0.25)
# cond2 = df_core_feature['CK density'] > df_core_feature['CK density'].quantile(0.3)
# cond1 = df_core_feature['tumor percent'] > 0.3
# cond2 = df_core_feature['CK percent'] > 0.2
# df_core_feature = df_core_feature.loc[cond2  ,]

# df_pleural = df_core_feature.loc[df_core_feature['casetype']=="Pleural", markers_density + genes_density]
# df_peritoneal= df_core_feature.loc[df_core_feature['casetype']=="Peritoneal",markers_density + genes_density]

# corr_pleural = df_pleural.corr()
# corr_pleural = corr_pleural.loc[genes_density, markers_density]
# corr_peritoneal = df_peritoneal.corr()
# corr_peritoneal = corr_peritoneal.loc[genes_density, markers_density]

# heatmap_plot(corr_pleural.round(2), "Pleural Intensity correlation in CK intensive cores",-0.6,0.6)

# #---------------------------------------
# # The third way of doing correlation




# from collections import defaultdict

# markers1 = ['FOXP3', 'CD4', 'CD8', 'CD68', 'CD20', 'CK']
# markers2 = ['CD11c', 'CD56', 'BAP1', 'NF2', 'MTAP'] #, 'LAG3']

# intensityCol = {}
# for marker in markers1:
#     if marker == "FOXP3":
#         intensityCol[marker] = "Entire Cell FOXP3 (Opal 570) Mean (Normalized Counts, Total Weighting)"
#     if marker == "CD4":
#         intensityCol[marker] = "Entire Cell CD4 (Opal 690) Mean (Normalized Counts, Total Weighting)"
#     if marker == "CD8":
#         intensityCol[marker] = "Entire Cell CD8 (Opal 480) Mean (Normalized Counts, Total Weighting)"
#     if marker == "CD20":
#         intensityCol[marker] =  "Entire Cell CD20 (Opal 620) Mean (Normalized Counts, Total Weighting)"
#     if marker == "CK":
#         intensityCol[marker] = "Entire Cell PanCK (Opal 780) Mean (Normalized Counts, Total Weighting)"
#     if marker == "CD68":
#         intensityCol[marker] = "Entire Cell CD68 (Opal 520) Mean (Normalized Counts, Total Weighting)"
# for marker in markers2:
#     if marker == "CD11c":
#         intensityCol[marker] = "Entire Cell CD11c (Opal 480) Mean (Normalized Counts, Total Weighting)"
#     if marker == "CD56":
#         intensityCol[marker] = "Entire Cell CD56 (Opal 520) Mean (Normalized Counts, Total Weighting)"
#     if marker == "NF2":
#         intensityCol[marker] = "Entire Cell NF2 (Opal 780) Mean (Normalized Counts, Total Weighting)"
#     if marker == "BAP1":
#         intensityCol[marker] = "Entire Cell BAP1 (Opal 690) Mean (Normalized Counts, Total Weighting)"
#     if marker == "MTAP":
#         intensityCol[marker] = "Entire Cell MTAP (Opal 620) Mean (Normalized Counts, Total Weighting)"
#     if marker == "LAG3":
#         intensityCol[marker] = "Entire Cell LAG3 (Opal 570) Mean (Normalized Counts, Total Weighting)"
        


# densities = defaultdict(lambda: [])
# percents = defaultdict(lambda: [])
# intensities = defaultdict(lambda: [])
# others = defaultdict((lambda: []))

# all_cores = set(df.loc[df['Classification']=='Malignant', 'Annotation ID'])

# # for k in list(set(df.loc[df['Classification']=='Malignant', 'Annotation ID'].values)):
# for k in all_cores:
#     k2 = df_map.loc[(df_map['Annotation ID']==k) ,'Annotation2'].values[0]

#     core1 = df_core.loc[df_core['Annotation ID'] == k, ]
#     core2 = df_core2.loc[df_core2['Annotation2'] == k2, ]
#     if (core2.shape[0] == 0): continue
    
#     area_stroma = df_map.loc[(df_map['Annotation ID']==k) , 'stromaArea'].values[0]
#     area_tumor = df_map.loc[df_map['Annotation ID']==k, 'tumorArea'].values[0]
#     area = area_stroma + area_tumor
#     core_casetype = df_map.loc[(df_map['Annotation ID']==k) , 'CaseType'].values[0]
#     core_subtype = df_map.loc[(df_map['Annotation ID']==k) , 'subtype'].values[0]
    
#     area2 = df_density2.loc[df_density2["Annotation2"]==k2, "Tissue Area (mm2)"]
#     if area2.empty:
#         area2 = area
#     else:
#         area2 = area2.values[0]
#         if area2 == 0:
#             area2 = area


    
#     n_tumor = sum(core1['Tissue Category'] == 'Tumor')
#     others['Tumor density'].append(np.log(n_tumor/area + sys.float_info.epsilon))
#     others['Tumor percent'].append(n_tumor/core1.shape[1])
#     others['Core ID'].append(k)
#     others['casetype'].append(core_casetype)
#     others['subtype'].append(core_subtype)
    
#     for marker in markers1:
        
#         n_marker = sum(core1['phenotype_combined'] == marker)
        
#         densities[marker].append(np.log(n_marker/area + sys.float_info.epsilon))
#         percents[marker].append(n_marker/core1.shape[0])
#         intensities[marker]. append(core1[intensityCol[marker]].mean())
    
#     for marker in['CD11c', 'CD56']:
 
#         n_marker = sum(core2['phenotype_combined'] == marker)
        
#         densities[marker].append(np.log(n_marker/area2 + sys.float_info.epsilon))
#         percents[marker].append(n_marker/core2.shape[0])
#         intensities[marker].append(core2[intensityCol[marker]].mean())

#     for marker in ['BAP1', 'NF2', 'MTAP']:
        
#         gene_count = sum( core2[f"Phenotype-{marker}"]==f"{marker}+") 
#         densities[marker].append(np.log(gene_count/area + sys.float_info.epsilon))
#         percents[marker].append(gene_count/core2.shape[0])
#         intensities[marker].append(core2[intensityCol[marker]].mean())
        
        
# df_density = pd.DataFrame(densities)
# df_percent = pd.DataFrame(percents)
# df_intensity = pd.DataFrame(intensities)

# df_others = pd.DataFrame(others)

# df_density = pd.concat([df_density,df_others],axis=1) 
# df_percent = pd.concat([df_percent,df_others],axis=1) 
# df_intensity = pd.concat([df_intensity,df_others],axis=1) 

# df_intensity.set_index('Core ID', inplace=True)
# df_percent.set_index('Core ID', inplace=True)
# df_density.set_index('Core ID', inplace=True)

# df_intensity.to_csv("./output/marker_intensity_malignant.csv")
# df_density.to_csv("./output/marker_density_malignant.csv")
# df_percent.to_csv("./output/marker_percent_malignant.csv")

# #------------------------------------------------
# # read back 
# df_intensity = pd.read_csv("./output/marker_intensity_malignant.csv", index_col=0)
# df_density = pd.read_csv("./output/marker_density_malignant.csv", index_col=0)
# df_percent = pd.read_csv("./output/marker_percent_malignant.csv", index_col=0)

# markers = ["CD20", "CD4", "CD8", "FOXP3", "CD68", "CD11c", "CD56", "CK"] 
# genes= ["BAP1", "NF2", "MTAP"]

# # corr_intensity = df_intensity.corr(method="spearman")
# # corr_intensity = corr_intensity.loc[genes, markers]
# # heatmap_plot(corr_intensity.round(2), "Intensity correlation in all malignant cores",-0.6,0.6)
# # plt.show()
# # corr_intensity.to_csv("./output/correlation/corr_intensity.csv")

# #select cores 
# cores_sel = df_percent.loc[df_percent['CK'] > 0.2,:].index.to_list()
# # cores_sel = df_intensity.loc[(df_intensity['CK'] > df_intensity['CK'].quantile(0.3)) ,:].index.to_list()

# studytype = "Intensity"
# df_type_sub = df_intensity.loc[cores_sel, :]

# df_intensity_sub = df_intensity.loc[cores_sel, markers + genes]
# corr_intensity_sub = df_intensity_sub.corr(method="spearman")
# corr_intensity_sub = corr_intensity_sub.loc[genes, markers]
# heatmap_plot(corr_intensity_sub.round(2), "Intensity correlation in CK intensive cores",-0.6,0.6)
# plt.show()
# corr_intensity_sub.to_csv("./output/correlation/corr_intensity_sub.csv")

# casetype = "Pleural"

# # cores_sel = df_intensity.loc[(df_intensity['casetype'] == casetype) ,:].index.to_list()

# # df_intensity_sub = df_intensity.loc[cores_sel, :]
# # corr_intensity_sub = df_intensity_sub.corr(method="spearman")
# # corr_intensity_sub = corr_intensity_sub.loc[genes, markers]
# # heatmap_plot(corr_intensity_sub, casetype + " Intensity correlation in all malignant cores",-0.6,0.6)
# # plt.show()

# cores_sel = df_percent.loc[(df_percent['CK'] > 0.2) & (df_percent['casetype'] == casetype) ,:].index.to_list()
# # cores_sel = df_intensity.loc[(df_intensity['CK'] > df_intensity['CK'].quantile(0.3)) & (df_intensity['casetype'] == casetype) ,:].index.to_list()

# df_intensity_sub = df_intensity.loc[cores_sel, markers + genes]
# corr_intensity_sub = df_intensity_sub.corr(method="spearman")
# corr_intensity_sub = corr_intensity_sub.loc[genes, markers]
# heatmap_plot(corr_intensity_sub.round(2), "Correlations for Pleural",-0.6,0.6)
# plt.show()


# casetype = "Peritoneal"

# cores_sel = df_percent.loc[(df_percent['CK'] > 0.2) & (df_percent['casetype'] == casetype) ,:].index.to_list()
# # cores_sel = df_intensity.loc[(df_intensity['CK'] > df_intensity['CK'].quantile(0.3)) & (df_intensity['casetype'] == casetype) ,:].index.to_list()

# df_intensity_sub = df_intensity.loc[cores_sel, markers + genes]
# corr_intensity_sub = df_intensity_sub.corr(method="spearman")
# corr_intensity_sub = corr_intensity_sub.loc[genes, markers]
# heatmap_plot(corr_intensity_sub.round(2),   "Correlation for Peritoneal",-0.45,0.45)
# plt.show()


# ################################################################
# # # scatter plot of correlation denoted on subtype
# # df_intensity = pd.read_csv("./output/marker_intensity_malignant.csv", index_col=0)
# # casetype = "Pleural"
# # cores_sel = df_percent.loc[(df_percent['CK'] > 0.2) & (df_percent['casetype'] == casetype) ,:].index.to_list()
# # df_intensity_sub = df_intensity.loc[cores_sel,:]
# # df_intensity_sub = df_intensity_sub.loc[df_intensity_sub['subtype'] != "Not specified", :]
# # marker = "CD56"
# # gene = "NF2"
# # corr = "0.67"

# # title = f"Pleural {gene}-{marker} corr={corr}"
# # scatter_plot(df_intensity_sub, marker, gene, 'subtype', title, save=True)


#################################################################
# correlation with the data genereated by r

# Panel 1
outpath = "../../Paper/Figure/Correlation/"

inte_type = "mean"
markers = ["CD20", "CD4", "CD8", "FOXP3", "CD68", "CD11c", "CD56", "CK"] 
genes= ["BAP1", "NF2", "MTAP", "LAG3"]

only3subs = False


grouptype = " intensity"
cols = [item + grouptype for item in markers+genes]

df_feature_ori = pd.read_csv("./output/cellContact/core_features_allMarkers_withIntensity_MVB_"+inte_type+".csv", index_col=0)
if only3subs:
    df_feature_ori = df_feature_ori.loc[df_feature_ori['subtype'].isin(["biphasic", "epithelioid", "sarcomatoid"]),:]
#104,42
# correlat_ion analysis based on CK patients
# df_feature_.replace(0, np.nan, inplace=True)
for CKcore in [False]:  # [False, True]:

    title = "Correlation for "
    df_feature = df_feature_ori
    
    if CKcore:
        df_feature= df_feature_ori.loc[df_feature_ori['CK percent'] > 0.2, :] #65,42
        title = "CK core " + title
        
    for casetype in ["Pleural", "Peritoneal"]:
        
        if only3subs:
            title =  title + casetype + " only3subs"
        else:
            title = title + casetype
       
        df_feature_sub = df_feature.loc[df_feature['casetype'] == casetype, :]
        
        
        df_feature_sel = df_feature_sub[cols]
        df_feature_sel.columns = markers+genes
        
        
        corr = df_feature_sel.corr(method="spearman")
        corr = corr.loc[genes, markers]
        
        if only3subs:
            corr.to_csv(f"{outpath}/{title}_only3subs.csv")
        else:
            corr.to_csv(f"{outpath}/{title}.csv")
            
        if casetype == "Pleural":
            heatmap_plot(corr.round(2), outpath,  title , -0.5,0.6, contain0=True)#-0.6,0.6
        else:
            heatmap_plot(corr.round(2), outpath,  title , -0.5,0.6, contain0=True)#-
            
        plt.show()
   
##################################################################
# scatter plot for all correaltions
CK_core = "CK core " #""
casetype = "Pleural"

subtypes={ "Pleural":  ["epithelioid", "biphasic" , "sarcomatoid"  ],
          "Peritoneal": ["epithelioid", "biphasic" ]
          }

grouptype = " intensity"
inte_type = "mean"
df_feature_ori = pd.read_csv("./output/cellContact/core_features_allMarkers_withIntensity_MVB_"+inte_type+".csv", index_col=0)

markers = ["CD20", "CD4", "CD8", "FOXP3", "CD68", "CD11c", "CD56", "CK"] 
genes= ["BAP1", "NF2", "MTAP", "LAG3"]

for CK_core in ["", "CK core "]:
    
    if CK_core != "":
        df_feature= df_feature_ori.loc[df_feature_ori['CK percent'] > 0.2, :]
        outpath = "../Paper/Figure/Correlation/ScatterPlot/CK cores"
    else:
        df_feature = df_feature_ori
        outpath = "../Paper/Figure/Correlation/ScatterPlot/All cores"
        
    for casetype in  ["Pleural", "Peritoneal"]:
        
        df_corr = pd.read_csv("../../Paper/Figure/Correlation/{CK_core}Correlation for {casetype}.csv", index_col=0)

        df_feature_sub = df_feature.loc[(df_feature['casetype'] == casetype) ,:]
        df_feature_sub = df_feature_sub.loc[df_feature_sub['subtype'] != "Not specified", :]
        
        for gene in genes:
            for marker in markers:
           
                corr = df_corr.loc[gene, marker]
                
                
                title = f"{CK_core}{casetype} {gene}-{marker} intensity correlation"
                
                
                
                scatter_plot(df_feature_sub , marker, gene, grouptype, 'subtype', subtypes[casetype], title,corr, outpath, save=True)


        


        

# correlation among panel2 genes    
inte_type = "mean"
df_feature_ori = pd.read_csv("./output/cellContact/core_features_allMarkers_withIntensity_MVB_"+inte_type+".csv", index_col=0)

# markers = ["CD20", "CD4", "CD8", "FOXP3", "CD68", "CD11c", "CD56", "CK"] 
genes= ["BAP1", "NF2", "MTAP", "LAG3"]


grouptype = " intensity"
cols = [item + grouptype for item in genes]


for CKcore in [True, False]:

    title = "Panel 2 marker correlation for "
    df_feature = df_feature_ori
    
    if CKcore:
        df_feature= df_feature_ori.loc[df_feature_ori['CK percent'] > 0.2, :] #65,42
        title = "CK core " + title
            
    for casetype in ["Pleural", "Peritoneal"]:
       
        df_feature_sub = df_feature.loc[df_feature['casetype'] == casetype, :]
        
        
        df_feature_sel = df_feature_sub[cols]
        df_feature_sel.columns = genes
        
        corr = df_feature_sel.corr(method="spearman")
        mask = np.zeros_like(corr, dtype=bool)
        mask[np.triu_indices_from(mask)] = True  #set lower left triangla as True
            
        if casetype == "Pleural":
            heatmap_plot(corr.round(2),   title+casetype ,-0,0.6, mask=mask)#-0.6,0.6
        else:
            heatmap_plot(corr.round(2),   title+casetype ,-0.2,0.6, mask=mask)#-

                    
        plt.show()


        
################################################################
# scatter plot of correlation denoted on subtype - data generated from R
inte_type = "mean"
df_feature_ori = pd.read_csv("./output/cellContact/core_features_allMarkers_withIntensity_MVB_"+inte_type+".csv", index_col=0)

grouptype = " percent" #" intensity"

subtypes={ "Pleural":  ["epithelioid", "biphasic" , "sarcomatoid"  ],
          "Peritoneal": ["epithelioid", "biphasic" ]
          }


casetype = "Peritoneal" #"Pleural"
df_feature_sub = df_feature_ori.loc[(df_feature_ori['casetype'] == casetype) ,:]
df_feature_sub = df_feature_sub.loc[df_feature_sub['subtype'] != "Not specified", :]

#pleural
marker = "CD56"
gene = "NF2"
corr = 0.61


marker = "CD8"
gene = "BAP1"
corr=-0.50

marker = "CD11c"
gene = "BAP1"
corr=-0.40

marker = "CD56"
gene = "MTAP"
corr=0.56

marker = "CD56"
gene = "LAG3"
corr=0.54   

#peritoneal
marker = "CD11c"
gene = "BAP1"
corr=-0.52

marker = "CD11c"
gene = "MTAP"
corr=-0.47

title = f"{casetype} {gene}-{marker} intensity correlation"
scatter_plot(df_feature_sub , marker, gene, grouptype, 'subtype', subtypes[casetype], title,corr, save=True)



################################################################################
# density compartion  between casetype, subtypes -- from R generated data
inte_type = "mean"
non_parametric=True

grouptype = "percent" # "intensity" # #"density" # # #"percent" # # # # #"density" #"percent"  # #
title = f'{grouptype.capitalize()} of Phenotype' 

# outpath = "../../Paper/Figure/BoxPlot_coreBased"
outpath = "../../Paper/Figure/BoxPlot"
# df_feature_ori = pd.read_csv("./output/cellContact/core_features_allMarkers_withIntensity_core.csv")

df_feature_ori = pd.read_csv("./output/cellContact/core_features_allMarkers_withIntensity_MVB_"+inte_type+".csv")
df_feature_ori.loc[(df_feature_ori['subtype']=='epithelial') | (df_feature_ori['subtype']=='epithelial or epithelioid'), 'subtype' ] = "epithelioid"

markers = ["CD20", "CD4", "CD8", "FOXP3", "CD68", "CD11c", "CD56", "CK"] # unidentified does not have intensity, "Unidentified"]

cols = [marker+" "+grouptype for marker in markers]
df_feature = df_feature_ori[['MVB','casetype','subtype'] + cols]
# change the column names (remove the grouptype)
df_feature.columns = ['MVB','casetype','subtype'] + markers
df_feature = df_feature.melt(id_vars=['MVB','casetype','subtype'], var_name='phenotype', value_name=grouptype)
df_feature.replace(0, np.nan, inplace=True) #do not plot the density of 0 cores
if grouptype == "density":
    df_feature['density'] = np.log(df_feature['density']+ 1e-16 )
    
# df_feature[grouptype] = min_max_norm(df_feature[grouptype])
# outpath = "../../Paper/Figure/BoxPlot/min_max_norm"  


#---- plot and pvalues for all tissue---------------------------------------------            
    
#Pleural and Peritoneal
casetype = None
hue_order=["Pleural", "Peritoneal"]
# generate_barplot_and_pvalues(df_feature,outpath, title, casetype, markers, hue_order, grouptype, non_parametric, size=(30, 10.5))
 
df_pvalue = generate_pvalues(df_feature, markers, "casetype", hue_order, grouptype, outpath, title, non_parametric=non_parametric)  
# bar_plot_plus(title,df_plot, 'phenotype', grouptype, 'casetype', markers, hue_order=hue_order)
box_plot_plus(outpath, title,df_feature, 'phenotype', grouptype, 'casetype', markers, hue_order=hue_order, size=(15, 10), short=True)
box_plot_plus_horizontal(outpath, title,df_feature, grouptype, 'phenotype',  'casetype', markers, hue_order=hue_order, size=(15,10), short=False)



# "Pleural
casetype = "Pleural"
title = f'{casetype} - {grouptype.capitalize()} of phenotype'
    
hue_order=["epithelioid", "biphasic", "sarcomatoid"]
# generate_barplot_and_pvalues(df_feature,outpath, title, casetype, markers, hue_order, grouptype, non_parametric, size=(17, 10.5))
df_plot = df_feature.loc[df_feature['casetype']==casetype,:]
df_plot = df_plot.loc[df_plot['subtype'].isin(hue_order), :]



df_pvalue = generate_pvalues(df_plot, markers, "subtype", hue_order, grouptype, outpath, title, non_parametric=non_parametric) 
# bar_plot_plus(title,df_plot, 'phenotype', grouptype, 'subtype', markers, hue_order=hue_order)
box_plot_plus(outpath, title,df_plot, 'phenotype', grouptype, 'subtype', markers, hue_order=hue_order, size=(17, 10.5), short=True)
box_plot_plus_horizontal(outpath, title+"_sw",df_plot, grouptype, 'phenotype', 'subtype', markers, hue_order=hue_order, size=(14,12), short=False)


casetype = "Peritoneal"
title = f'{casetype} - {grouptype.capitalize()} of phenotype'
hue_order=["epithelioid", "biphasic"]
# generate_barplot_and_pvalues(df_feature,outpath, title, casetype, markers, hue_order, grouptype, non_parametric, size=(17, 10.5))

df_plot = df_feature.loc[df_feature['casetype']==casetype,:]
df_plot = df_plot.loc[df_plot['subtype'].isin(hue_order), :]

generate_pvalues(df_plot, markers, "subtype", hue_order, grouptype, outpath, title, non_parametric=non_parametric) 
# bar_plot_plus(title,df_plot, 'phenotype', grouptype, 'subtype', markers, hue_order=hue_order)
box_plot_plus(outpath, title,df_plot, 'phenotype', grouptype, 'subtype', markers, hue_order=hue_order, size=(17, 10.5), short=True)
box_plot_plus_horizontal(outpath, title+"_sw",df_plot, grouptype, 'phenotype', 'subtype', markers, hue_order=hue_order, size=(14,12), short=False)


# compare the subtypes in between Pleural and Peritoneal
# eg. pl-"epithelioid" vs pe-"epithelioid"     
# eg. pl-"biphasic" vs pe-"biphasic"

from statsmodels.stats.multitest import multipletests

pvals_ttest_all = []
pvals_wilcox_all = []

pvals_ttest_all_adj = []
pvals_wilcox_all_adj = []

subtypes_list = []
for subtype in ["epithelioid", "biphasic"]:
    pvals_ttest = []
    pvals_wilcox = []
    for marker in markers:
    
        pl = df_feature.loc[(df_feature['casetype'] == "Pleural") & (df_feature['subtype']==subtype) &(df_feature['phenotype']==marker), grouptype ]
        pe = df_feature.loc[(df_feature['casetype'] == "Peritoneal") & (df_feature['subtype']==subtype)&(df_feature['phenotype']==marker), grouptype ]
        pval_ttest = get_pvalue(pl.dropna(), pe.dropna(), non_parametric=False) 
        pval_wilcox = get_pvalue(pl.dropna(), pe.dropna(), non_parametric=True)
        pvals_ttest.append(float(pval_ttest))
        pvals_wilcox.append(float(pval_wilcox))
        subtypes_list.append(subtype)
    
    # pvals_ttest_adj = 
    pvals_ttest_all = pvals_ttest_all + pvals_ttest
    pvals_wilcox_all = pvals_wilcox_all + pvals_wilcox
    
    pvals_ttest_all_adj = pvals_ttest_all_adj + list(multipletests(pvals_ttest, method="fdr_bh")[1])
    pvals_wilcox_all_adj = pvals_wilcox_all_adj + list(multipletests(pvals_wilcox, method="fdr_bh")[1])
    
df_pvalues = pd.DataFrame({
    "subtype": subtypes_list,
    "markers": markers + markers,   
    "ttest_pvals": pvals_ttest_all,
    "ttest_pvals_adj": pvals_ttest_all_adj,
    "wilcox_pvals": pvals_wilcox_all,
    "wilcox_pvals_adj": pvals_wilcox_all_adj}
    
    )
        
df_pvalues.to_csv(f"../../Paper/Figure/NEW_ORGANIZATION/Phenotype/MPM_vs_MPeM/subtype_difference_{grouptype}.csv", index=False)
    
#---- plot and pvalues for tunor and stroma separately---------------------------------------------  
# generat MVB feature, only execute once
# for tissue in ["Tumor", "Stroma"]:    
#     df_feature_core = pd.read_csv(f"./output/cellContact/core_features_panel1Marker_{tissue}_core.csv")
#     df_feature_MVB = df_feature_core.groupby(['MVB','casetype','subtype']).mean().reset_index()
#     df_feature_MVB.to_csv(f"./output/cellContact/core_features_panel1Marker_{tissue}_MVB_mean.csv")

#     df_feature_MVB = df_feature_core.groupby(['MVB','casetype','subtype']).median().reset_index()
#     df_feature_MVB.to_csv(f"./output/cellContact/core_features_panel1Marker_{tissue}_MVB_median.csv")
    
markers = ["CD20", "CD4", "CD8", "FOXP3", "CD68", "CK", "Unidentified"]
grouptype = " density"
inte_type = "mean"

for tissue in ["Tumor", "Stroma"]:    
    df_feature_ori = pd.read_csv(f"./output/cellContact/core_features_panel1Marker_{tissue}_MVB_{inte_type}.csv")

   
    
    cols = [marker+grouptype for marker in markers]
    df_feature = df_feature_ori[['MVB','casetype','subtype'] + cols]
    df_feature.columns = ['MVB','casetype','subtype'] + markers
    df_feature = df_feature.melt(id_vars=['MVB','casetype','subtype'], var_name='phenotype', value_name='density')
    df_feature.replace(0, np.nan, inplace=True)
    df_feature['density'] = np.log(df_feature['density'])
    
    
    title = 'Density of Phenotype for ' + tissue    
    non_parametric=True    
        
    #Pleural and Peritoneal
    casetype = None
    hue_order=["Pleural", "Peritoneal"]
    generate_barplot_and_pvalues(df_feature,title, casetype, markers, hue_order, non_parametric)
    
    # "Pleural
    casetype = "Pleural"
    title = f'{casetype} - Density of phenotype for {tissue}'
    hue_order=["epithelioid", "biphasic", "sarcomatoid"]
    generate_barplot_and_pvalues(df_feature,title, casetype, markers, hue_order, non_parametric)
    
    casetype = "Peritoneal"
    title = f'{casetype} - Density of phenotype for {tissue}'
    hue_order=["epithelioid", "biphasic"]
    generate_barplot_and_pvalues(df_feature,title, casetype, markers, hue_order, non_parametric)

        
    


# --------------------------------------------------------------------------- 
# compare between combined groups
## -------------------------------------------
# compare in total
df_merged = df_type_sub

BAP1_q50 = df_merged['BAP1'].quantile(0.5)
NF2_q50 = df_merged['NF2'].quantile(0.5)
MTAP_q50 = df_merged['MTAP'].quantile(0.5)

df_merged['BAP1status'] = np.where(df_merged['BAP1'] > BAP1_q50, 'BAP1+', 'BAP1-')
df_merged['NF2status'] = np.where(df_merged['NF2'] > NF2_q50, 'NF2+', 'NF2-')
df_merged['MTAPstatus'] = np.where(df_merged['MTAP'] > MTAP_q50, 'MTAP+', 'MTAP-')
df_merged['status'] = df_merged['BAP1status'] + df_merged['NF2status'] + df_merged['MTAPstatus']

''
# total
df_sub = df_merged.loc[(df_merged['status'] == 'BAP1-NF2-MTAP-') | (df_merged['status'] == 'BAP1-NF2+MTAP+'), ['CD11c', 'CD20', 'CD4', 'CD56', 'CD68', 'CD8',
       'CK', 'FOXP3', 'status']]

# df_sub = df_merged[['CD11c', 'CD20', 'CD4', 'CD56', 'CD68', 'CD8',
#        'CK', 'FOXP3', 'status']]

df_sub_plot = df_sub.melt(id_vars = 'status', var_name='phenotype', value_name=studyType)

bar_plot(studyType+' comparison ',df_sub_plot, 'phenotype', studyType, 'status')
plt.show()
# Pleural
df_sub_pleural = df_merged.loc[(df_merged['casetype']=="Pleural") & ((df_merged['status'] == 'BAP1-NF2-MTAP-') | (df_merged['status'] == 'BAP1-NF2+MTAP+')) , ['CD11c', 'CD20', 'CD4', 'CD56', 'CD68', 'CD8',
       'CK', 'FOXP3', 'status']]

# df_sub_pleural = df_merged.loc[(df_merged['casetype']=="Pleural")  , ['CD11c', 'CD20', 'CD4', 'CD56', 'CD68', 'CD8',
#        'CK', 'FOXP3', 'status']]

df_sub_plot = df_sub_pleural.melt(id_vars = 'status', var_name='phenotype', value_name=studyType)

bar_plot(studyType+' comparison in Pleural ',df_sub_plot, 'phenotype', studyType, 'status')
plt.show()
# Perttoneal
df_sub_peritoneal = df_merged.loc[(df_merged['casetype']=="Peritoneal") & ((df_merged['status'] == 'BAP1-NF2-MTAP-') | (df_merged['status'] == 'BAP1-NF2+MTAP+')) , ['CD11c', 'CD20', 'CD4', 'CD56', 'CD68', 'CD8',
       'CK', 'FOXP3', 'status']]

# df_sub_peritoneal = df_merged.loc[(df_merged['casetype']=="Peritoneal") , ['CD11c', 'CD20', 'CD4', 'CD56', 'CD68', 'CD8',
#        'CK', 'FOXP3', 'status']]

df_sub_plot = df_sub_peritoneal.melt(id_vars = 'status', var_name='phenotype', value_name=studyType)

bar_plot(studyType+' comparison in Peritoneal ',df_sub_plot, 'phenotype', studyType, 'status')
plt.show()

#####################################################################################
# check correlation of each pair for the subtypes
#add subtype information
pleural_cores = list(df_pleural.index)
pleural_cores_subtype = df_map.loc[df_map['Annotation ID'].isin(pleural_cores), ['Annotation ID','subtype']]
pleural_cores_subtype.set_index('Annotation ID', inplace=True)

peritoneal_cores = list(df_peritoneal.index)
peritoneal_cores_subtype = df_map.loc[df_map['Annotation ID'].isin(peritoneal_cores), ['Annotation ID','subtype']]
peritoneal_cores_subtype.set_index('Annotation ID', inplace=True)

df_plot_pleural = pd.concat([df_pleural,pleural_cores_subtype], axis=1)
df_plot_peritoneal = pd.concat([df_peritoneal,peritoneal_cores_subtype], axis=1) 


for marker in markers:
    for gene in genes:
        corr = round(df_plot_pleural[marker].corr(df_plot_pleural[gene]),2)
        title = f"Pleural {gene}-{marker} corr={corr}"
        scatter_plot(df_plot_pleural, marker, gene, 'subtype', title, save=True)
        plt.show()
        
        corr = round(df_plot_peritoneal[marker].corr(df_plot_peritoneal[gene]),2)
        title = f"Peritoneal {gene}-{marker} corr={corr}"
        scatter_plot(df_plot_peritoneal, marker, gene, 'subtype', title, save=True)
        plt.show()
        
        
#----------------------------------------------------
df2_malignant = df2.loc[df2['Classification']=='Malignant',:]
gene_percent = {}
for gene in genes:
    gene_percent[gene] =  sum( df2_malignant[f"Phenotype-{gene}"]==f"{gene}+") /df2_malignant.shape[0]

bar_plot_series(pd.Series(gene_percent), "Gene assignment percentage")


############################################################
# check the panel2 cores for BAP1, NF2, MTAP
check = "BAP1" #, "NF2", "MTAP"]

df2_tumor = df2.loc[df2['Annotation ID'].isin(cores_sel), :] #select the cores which tumor percentage > 0.25
df2_tumor = df2

# BAP1_percents = []
# BAP1_tumor_ratios = []

gb = df2_tumor[['Annotation2', f"Phenotype-{check}"]].groupby('Annotation2')

# percent_cores_pos = []
# percent_cores_neg = []

# ratio_cores_pos = []
# ratio_cores_neg = []
 
BAP1_densities = []
tumor_densities = []

density_cores_pos = []
density_cores_neg = []

for k, gp in gb:
    k1 = df_map.loc[df_map["Annotation2"]==k, "Annotation ID"].values[0]
    core1 = df.loc[df['Annotation ID']==k1 ,:]
    
    area2 = df2.loc[df2["Annotation2"]==k, "Tissue Area (mm2)"].values[0]
    area1_tumor = core1.loc[:, "tumorArea"].values[0]
    area1_stroma = core1.loc[:, "stromaArea"].values[0]
    area1 = area1_tumor +  area1_stroma
    
    nBAP1 = sum(gp[f"Phenotype-{check}"]==f"{check}+")
    nTumor = sum(core1['Tissue Category']=="Tumor")
    density_BAP1 = np.log(nBAP1/area2 + sys.float_info.epsilon)
    density_tumor = log(nTumor/area1 + sys.float_info.epsilon)
    
    BAP1_densities.append(density_BAP1)
    tumor_densities.append(density_tumor)
    
    
    if density_tumor < tumor_density_cut: continue
    
    if density_BAP1 > density_high:
        density_cores_pos.append(k1)
    if density_BAP1 < density_low:
        density_cores_neg.append(k1)
    
    # tumor_percent = round(sum(core1['Tissue Category'] == "Tumor")/core1.shape[0],3)
    # BAP1_percent = round(sum(gp.loc[:,f"Phenotype-{check}"]==f"{check}+")/gp.shape[0],3)
    # BAP1_ratio = round(BAP1_percent/tumor_percent,3)
    # BAP1_tumor_ratios.append(BAP1_ratio)
    # BAP1_percents.append(BAP1_percent)
   
    # if BAP1_percent > percent_high:
    #     percent_cores_pos.append(k1)
    # if BAP1_percent < percent_low:
    #     percent_cores_neg.append(k1)
        
    # if BAP1_ratio > ratio_high:
    #     ratio_cores_pos.append(k1)
    # if BAP1_ratio < ratio_low:
    #     ratio_cores_neg.append(k1)


tumor_density_cut = pd.Series(tumor_densities).quantile(0.4)   

density_low = pd.Series(BAP1_densities).quantile(0.5)
density_high = pd.Series(BAP1_densities).quantile(0.5)



path_out = './output/cellContact/data/'

np.savetxt(path_out + "density_cores_pos.csv", np.array(density_cores_pos), fmt="%s")
np.savetxt(path_out + "density_cores_neg.csv", np.array(density_cores_neg), fmt="%s")

# np.savetxt(path_out + "percent_cores_pos.csv", np.array(percent_cores_pos), fmt="%s")
# np.savetxt(path_out + "percent_cores_neg.csv", np.array(percent_cores_neg), fmt="%s")
# np.savetxt(path_out + "ratio_cores_pos.csv", np.array(ratio_cores_pos), fmt="%s")
# np.savetxt(path_out + "ratio_cores_neg.csv", np.array(ratio_cores_neg), fmt="%s")


sns.histplot(BAP1_percents, bins=100)
plt.title("BAP1 percents in the core")
plt.show()


sns.histplot(BAP1_tumor_ratios, bins=100)
plt.title("BAP1 tumor ratios in the core")
plt.show()

ratio_low = pd.Series(BAP1_tumor_ratios).quantile(0.3)
ratio_high = pd.Series(BAP1_tumor_ratios).quantile(0.7)

percent_low = pd.Series(BAP1_percents).quantile(0.3)
percent_high = pd.Series(BAP1_percents).quantile(0.7)


#########################################################################################
#check BAP1 density
# panel2 counts from statistics of df2
BAP1_density = []
tumor_density = []
for k in list(set(df2["Annotation2"])):

    area = df_density2.loc[df_density2["Annotation2"]==k, "Tissue Area (mm2)"].values[0]
    
    BAP1_density.append(np.log((sum(df2.loc[(df2['Annotation2']==k) ,f"Phenotype-{check}"] ==f"{check}+")/area + sys.float_info.epsilon))
    
    k1 = df_map.loc[df_map["Annotation2"]==k, "Annotation ID"].values[0]
    
    area_stroma = df_map.loc[(df_map['Annotation ID']==k1) , 'stromaArea'].values[0]
    area_tumor = df_map.loc[df_map['Annotation ID']==k1, 'tumorArea'].values[0]
    area = area_stroma + area_tumor
    
    tumor_density.append(np.log(sum(df.loc[df['Annotation ID']==k1 ,'Tissue Category'] =="Tumor")/area + sys.float_info.epsilon))
    
stats.pearsonr(tumor_density, BAP1_density)
#(0.16835083744745277, 0.00229013664267726)

#########################################################################################
# heatmap cell-cell contact on mean contacts volcano plot criteria
def normalize_df(df):
    df = (df-df.min().min())/(df.max().max()-df.min().min())
    return(df)

path_in = "../../Paper/Figure/CellContact/Heatmap/"

region = "tumor"
casetype = "allCore" 
file = f"cellContact_score_mean_{casetype}_{region}Region_5types.csv"

df = pd.read_csv(path_in + file, index_col=0)
df = normalize_df(df)

title = file[:-4]
heatmap_plot(df, path_in, title, 0,1)

region = "tumor"
casetype = "Pleural" 
file = f"cellContact_score_mean_{casetype}_{region}Region_5types.csv"
df = pd.read_csv(path_in + file, index_col=0)
df = normalize_df(df)
title = file[:-4]
heatmap_plot(df, path_in, title, 0,1)

region = "tumor"
casetype = "Peritoneal" 
file = f"cellContact_score_mean_{casetype}_{region}Region_5types.csv"
df = pd.read_csv(path_in + file, index_col=0)
df = normalize_df(df)
title = file[:-4]
heatmap_plot(df, path_in, title, 0,1)

#---------------------------------------------
region = "all"
casetype = "allCore" 
file = f"cellContact_score_mean_{casetype}_{region}Region_5types.csv"
df = pd.read_csv(path_in + file, index_col=0)
df = normalize_df(df)
title = file[:-4]
heatmap_plot(df, path_in, title, 0,1)

region = "all"
casetype = "Pleural" 
file = f"cellContact_score_mean_{casetype}_{region}Region_5types.csv"
df = pd.read_csv(path_in + file, index_col=0)
df = normalize_df(df)
title = file[:-4]
heatmap_plot(df, path_in, title, 0,1)

region = "all"
casetype = "Peritoneal" 
file = f"cellContact_score_mean_{casetype}_{region}Region_5types.csv"
df = pd.read_csv(path_in + file, index_col=0)
df = normalize_df(df)
title = file[:-4]
heatmap_plot(df, path_in, title, 0,1)




#########################################################################################
# heatmap cell-cell contact on median results 


path_in = "./output/cellContact/tumor_region_contactScore/"

file= "cellContact_score_median_pleural_tumorRegion.csv"
df = pd.read_csv(path_in + file, index_col=0)
title = "Pleural tumorRegion celltype contact score"
heatmap_plot(df, title, -10, 5) #, 0.9, 8.5

file= "cellContact_score_median_peritoneal_tumorRegion.csv"
df = pd.read_csv(path_in + file, index_col=0)
title = "Peritoneal tumorRegion celltype contact score"
heatmap_plot(df, title,  -10, 5) #, 0.9, 8.5

file= "cellContact_score_median_all_tumorRegion.csv"
df = pd.read_csv(path_in + file, index_col=0)
title = "Allcores tumorRegion celltype contact score"
heatmap_plot(df, title, -10, 5) #, 0.9, 8.5

#........................................................
# allregion

path_in = "./output/cellContact/all_region_contactScore/"

file= "cellContact_score_median_Pleural_allRegion.csv"
df = pd.read_csv(path_in + file, index_col=0)
title = "Pleural allRegion celltype contact score"
heatmap_plot(df, title,-10, 6) #, 0.9, 8.5

file= "cellContact_score_median_Peritoneal_allRegion.csv"
df = pd.read_csv(path_in + file, index_col=0)
title = "Peritoneal allRegion celltype contact score"
heatmap_plot(df, title,-10, 5) #, 0.9, 8.5

file= "cellContact_score_median_all_allRegion.csv"
df = pd.read_csv(path_in + file, index_col=0)
title = "Allcores allRegion celltype contact score"
heatmap_plot(df, title,-10, 5) #, 0.9, 8.5

#........................................................
# allregion type exist

path_in = "./output/cellContact/all_region_contactScore/"

file= "cellContact_score_median_Pleural_allRegion_typeExist.csv"
df = pd.read_csv(path_in + file, index_col=0)
title = "Pleural allRegion celltype contact score"
heatmap_plot(df, title,-1, 8, contain0=False) #, 0.9, 8.5

file= "cellContact_score_median_Peritoneal_allRegion_typeExist.csv"
df = pd.read_csv(path_in + file, index_col=0)
title = "Peritoneal allRegion celltype contact score"
heatmap_plot(df, title,-1, 8, contain0=False) #, 0.9, 8.5

file= "cellContact_score_median_all_allRegion.csv"
df = pd.read_csv(path_in + file, index_col=0)
title = "Allcores allRegion celltype contact score"
heatmap_plot(df, title,-10, 5) #, 0.9, 8.5

#........................................................
# allregion type exist

path_in = "./output/cellContact/all_region_contactScore/"

file= "cellContact_score_median_Pleural_allRegion_type0.01.csv"
df = pd.read_csv(path_in + file, index_col=0)
title = "Pleural allRegion celltype contact score"
heatmap_plot(df, title,0, 7, contain0=False) #, 0.9, 8.5

file= "cellContact_score_median_Peritoneal_allRegion_type0.01.csv"
df = pd.read_csv(path_in + file, index_col=0)
title = "Peritoneal allRegion celltype contact score"
heatmap_plot(df, title,0, 7, contain0=False) #, 0.9, 8.5


#........................................................
# all region, filtered cores (contact_score > -10)

cellContact_score_median_Peritoneal_allRegion_filtered.csv"
path_in = "./output/cellContact/all_core_contactScore/"

file= "cellContact_score_median_Pleural_allRegion_filtered1.csv"
df = pd.read_csv(path_in + file, index_col=0)
title = "Pleural allRegion celltype contact score vialid cores"
heatmap_plot(df, title,2,9.5, contain0=False) #, 0.9, 8.5

file= "cellContact_score_median_Peritoneal_allRegion_filtered1.csv"
df = pd.read_csv(path_in + file, index_col=0)
title = "Peritoneal allRegion celltype contact score valid cores"
heatmap_plot(df, title,2,9.5, contain0=False) #, 0.9, 8.5

file= "cellContact_score_median_all_allRegion_filtered1.csv"
df = pd.read_csv(path_in + file, index_col=0)
title = "Allcores allRegion celltype contact score valid cores"
heatmap_plot(df, title, 2,9.5, contain0=False) #, 0.9, 8.5
#########################################################################################
# read in cell-cell contact results and plot the heatmap

path_in = "./output/cellContact/"
file= "cellContact_scores_stratified_Pleural3.csv"
df = pd.read_csv(path_in + file, index_col=0)
title = "Pleural celltype contact score\n"
heatmap_plot(df, title, 0.9, 8.5)


path_in = "./output/cellContact/"
file= "cellContact_scores_stratified_Peritoneal3.csv"
df = pd.read_csv(path_in + file, index_col=0)
title = "Peritoneal celltype contact score\n"
heatmap_plot(df, title,  0.9, 8.5)

# reshape cell-cell contact matrix for circular plots
# https://app.flourish.studio/visualisation/13417798/edit
path_in = "./output/cellContact/"
file= "Pleural_Tumor_celltype_contact.csv"
df = pd.read_csv(path_in + file, index_col=0)
df_long = df.stack().reset_index()
df_long.columns = ["type1","type2", "score"]
df_long['score'] = df_long['score'] +5
df_long.to_csv(path_in + "Pleural_Tumor_celltype_contact_reshaped.csv", index=False)

path_in = "./output/cellContact/"
file= "Pleural_Stroma_celltype_contact.csv"
df = pd.read_csv(path_in + file, index_col=0)
df_long = df.stack().reset_index()
df_long.columns = ["type1","type2", "score"]
df_long['score'] = df_long['score'] +5
df_long.to_csv(path_in + "Pleural_Stroma_celltype_contact_reshaped.csv", index=False)

#-------------------------------------------------------------------
path_in = "./output/cellContact/"
file= "Peritoneal_Tumor_celltype_contact.csv"
df = pd.read_csv(path_in + file, index_col=0)
df_long = df.stack().reset_index()
df_long.columns = ["type1","type2", "score"]
df_long['score'] = df_long['score'] +5
df_long.to_csv(path_in + "Peritoneal_Tumor_celltype_contact_reshaped.csv", index=False)

path_in = "./output/cellContact/"
file= "Peritoneal_Stroma_celltype_contact.csv"
df = pd.read_csv(path_in + file, index_col=0)
df_long = df.stack().reset_index()
df_long.columns = ["type1","type2", "score"]
df_long['score'] = df_long['score'] +5
df_long.to_csv(path_in + "Peritoneal_Stroma_celltype_contact_reshaped.csv", index=False)

#-------------------------------------------------------------------
path_in = "./output/cellContact/"
file="cellContact_scores_stratified_Pleural3.csv"
df = pd.read_csv(path_in + file, index_col=0)
df_long = df.stack().reset_index()
df_long.columns = ["type1","type2", "score"]

df_long.to_csv(path_in + "Pleural_celltype_contact_reshaped.csv", index=False)

path_in = "./output/cellContact/"
file= "cellContact_scores_stratified_Peritoneal3.csv"
df = pd.read_csv(path_in + file, index_col=0)
df_long = df.stack().reset_index()
df_long.columns = ["type1","type2", "score"]

df_long.to_csv(path_in + "Peritoneal_celltype_contact_reshaped.csv", index=False)

#-----------------------------------------------------------------------
path_in = "./output/cellContact/"
file="cellContact_scores_stratified_BAP1-NF2+MTAP+_0.01_200_800.csv"
df = pd.read_csv(path_in + file, index_col=0)
title = "celltype contact score BAP1-NF2+MTAP+\n"
heatmap_plot(df, title, 1.5, 7)
plt.show()

path_in = "./output/cellContact/"
file="cellContact_scores_stratified_BAP1+NF2-MTAP+_0.01_200_800.csv"
df = pd.read_csv(path_in + file, index_col=0)
title = "celltype contact score BAP1+NF2-MTAP+\n"
heatmap_plot(df, title, 1.5, 7)
plt.show()

path_in = "./output/cellContact/"
file="cellContact_scores_stratified_BAP1+NF2+MTAP+_0.01_200_800.csv"
df = pd.read_csv(path_in + file, index_col=0)
title = "celltype contact score BAP1+NF2+MTAP+\n"
heatmap_plot(df, title, 1.5, 7)
plt.show()
#-------------------------------------------------------------------
# Compare Pleural and Peritoneal
path_in = "./output/cellContact/"

file1="cellContact_scores_CK0.01_Pleural_Tumor.csv"
df1 = pd.read_csv(path_in + file1, index_col=0)
df_long1 = df1.stack().reset_index()
df_long1.columns = ["type1","type2", "score"]
df_long1['casetype'] = "Pleural"

# create a mask to identify rows with duplicate features as mentioned above
mask_dups = (df_long1[['type1', 'type2']].apply(frozenset, axis=1).duplicated()) 
df_long1 = df_long1[~mask_dups]



file2= "cellContact_scores_CK0.01_Peritoneal_Tumor.csv"
df2 = pd.read_csv(path_in + file2, index_col=0)
df_long2 = df2.stack().reset_index()
df_long2.columns = ["type1","type2", "score"]
df_long2['casetype'] = "Peritoneal"

mask_dups = (df_long2[['type1', 'type2']].apply(frozenset, axis=1).duplicated()) 
df_long2 = df_long2[~mask_dups]



df_long = pd.concat([df_long1, df_long2], axis = 0)
df_long['type'] = df_long['type1'] + "_" + df_long['type2']

mask_same = df_long['type1']==df_long['type2']

df_sameType = df_long[mask_same]
df_notSameType = df_long[~mask_same]

df_immune2CK = df_notSameType[df_notSameType['type2']=='CK']

df_immune2immune = df_notSameType[df_notSameType['type2']!='CK']

#
title = "Immune to immune cell type contacts"
bar_plot(title, df_immune2immune, 'type','score', 'casetype')
plt.show()
#
title = "Immune to CK cell type contacts"
bar_plot(title, df_immune2CK, 'type','score', 'casetype')
plt.show()
#
title = "Celltype vs itself contacts"
bar_plot(title, df_sameType, 'type','score', 'casetype')
plt.show()



###############################################################
# Composition waterfall plots

casetype = "Pleural" #"Peritoneal" #"Pleural" #"Peritoneal"
# casetype = "Pleural"
# casetype = None
#"sarcomatoid"# "epithelioid"   #"sarcomatoid" "biphasic"
subtype_dict = {}
subtype_dict['Pleural'] = ["epithelioid", "biphasic", "sarcomatoid"]
subtype_dict['Peritoneal'] = ["epithelioid", "biphasic"]



for casetype in [None ,"Pleural", "Peritoneal"]:
        if casetype is None:
            data1 = df.loc[  (df['Classification']=='Malignant') ,:]
            data2 = df2.loc[  (df2['Classification']=='Malignant'),:]
        else:
            data1 = df.loc[  (df['Classification']=='Malignant') & (df['CaseType']==casetype) ,:]
            data2 = df2.loc[  (df2['Classification']=='Malignant') & (df2['CaseType']==casetype) ,:]
        
        waterfall_plot(data1, data2, casetype, by="Annotation ID")

for casetype in ["Pleural", "Peritoneal"]:    
    for subtype in subtype_dict[casetype]:
        
        data1 = df.loc[  (df['Classification']=='Malignant') & (df['CaseType']==casetype) &(df['subtype']==subtype) ,:]
        data2 = df2.loc[  (df2['Classification']=='Malignant') & (df2['CaseType']==casetype) ,:]

        waterfall_plot(data1, data2, casetype, subtype, by="Annotation ID")
    
    
colours = ['tab:brown', 'tab:olive', 'tab:purple', 'tab:red', 'tab:green', 'tab:blue', 'tab:pink', 'tab:orange', 'tab:gray']
celltypes = ['CD20',     'CD4', '       CD8',      'FOXP3',      'CD68',     'CD11c',     'CD56',     'CK',       'other']


################################################################################
# generate data for waterfall plot in R
outpath = "./output/waterfall_plot/"

data1 = df.loc[  (df['Classification']=='Malignant') ,:]
data2 = df2.loc[  (df2['Classification']=='Malignant'),:]

by = "MVB"

fc_p1 = data1.groupby([by]).apply(lambda x: x['phenotype_combined'].value_counts(sort=False, normalize=True).rename_axis('Type').reset_index(name='Percent'))
fc_p2 = data2.groupby([by]).apply(lambda x: x['phenotype_combined'].value_counts(sort=False, normalize=True).rename_axis('Type').reset_index(name='Percent'))


df_fc1 = fc_p1.reset_index()
df_fc2 = fc_p2.reset_index()

# pivot celltype to combine panel1 and panel2, panel1[other] -  panel2[CD56+CD11c]
df_percent1 = df_fc1.set_index([by, 'Type'])['Percent'].unstack().reset_index() #pivot dataframe unstack used to reshape dataframe
df_percent2 = df_fc2.set_index([by, 'Type'])['Percent'].unstack().reset_index() #colnames ['MVB', 'CD11c', 'CD56', 'Unidentified']
df_percent2.drop(['other'], axis=1, inplace=True) # delete column other
# combine df_percent and df_percent2, deduct CD56 and CD11c percents from other in df_percent

df_percent = pd.merge(df_percent1, df_percent2, on=by)

df_percent['Unidentified'] = df_percent['Unidentified'] - df_percent['CD56'] - df_percent['CD11c']
# remove MVB that Unidentified < 0
df_percent = df_percent.loc[df_percent['Unidentified'] >= 0, :]

# return back the long form and add on casetype and subtype
df_full = pd.merge(df_percent, data1[['MVB', 'CaseType', 'subtype']].drop_duplicates(), on="MVB", how="left" )

# only select MVB where Pleural have subtypes = ["epithelioid", "biphasic", "sarcomatoid"] and
# Pertoneal =   ["epithelioid", "biphasic"]

case_cond = df_full['CaseType'].isin(['Pleural', 'Peritoneal'])
pleural_cond = (df_full['CaseType'] == "Pleural") & (df_full['subtype'].isin(["epithelioid", "biphasic", "sarcomatoid"]))
peritoneal_cond =  (df_full['CaseType'] == "Peritoneal") &  (df_full['subtype'].isin(["epithelioid", "biphasic"])) 

df_sub = df_full.loc[case_cond,:]
df_sub = df_sub.loc[pleural_cond | peritoneal_cond, :]      #77 12                                                     
                                                             

df_out_wide = df_sub.sort_values(['CaseType','subtype','Unidentified'],ascending=False).groupby(['CaseType','subtype']).head(100)

df_out_long = df_sub.melt(id_vars=['MVB', 'CaseType', 'subtype'], var_name='Phenotype', value_name='Percent')

df_out_wide.to_csv(outpath + "waterfall_compostion_by_" + by + "_wide.csv", index=False)
df_out_long.to_csv(outpath + "waterfall_compostion_by_" + by + "_long.csv", index=True)

#############################################################################
# celltype composition differences between subtypes

subtype_dict = {}
subtype_dict['Pleural'] = ["epithelioid", "biphasic", "sarcomatoid"]
subtype_dict['Peritoneal'] = ["epithelioid", "biphasic"]



data1 = df.loc[ (df['Classification']=='Malignant') ,:]
data2 = df2.loc[ (df2['Classification']=='Malignant'),:]

by = "Annotation ID"
fc_p1 = data1.groupby([by]).apply(lambda x: x['phenotype_combined'].value_counts(sort=False, normalize=True).rename_axis('Type').reset_index(name='Percent'))
fc_p2 = data2.groupby([by]).apply(lambda x: x['phenotype_combined'].value_counts(sort=False, normalize=True).rename_axis('Type').reset_index(name='Percent'))


df_fc1 = fc_p1.reset_index()
df_fc2 = fc_p2.reset_index()

df_percent1 = df_fc1.set_index([by, 'Type'])['Percent'].unstack().reset_index()
df_percent2 = df_fc2.set_index([by, 'Type'])['Percent'].unstack().reset_index() #colnames ['MVB', 'CD11c', 'CD56', 'other']
df_percent2.drop(['other'], axis=1, inplace=True) # delete column other

# combine df_percent and df_percent2, deduct CD56 and CD11c percents from other in df_percent
df_percent = pd.merge(df_percent1, df_percent2, on=by)
df_percent['other'] = df_percent['other'] - df_percent['CD56'] - df_percent['CD11c']
# remove MVB that other < 0
df_percent = df_percent.loc[df_percent['other'] >= 0, :]

df_percent.index = df_percent[by]
df_percent.drop([by], axis=1, inplace=True)
df_percent.index.name= None
df_percent.columns.name = None


import itertools
type_list = []
for casetype in ['Pleural', 'Peritoneal']:
    for subtype in subtype_dict[casetype]:
        type_list.append(f"{casetype}-{subtype}")
pair_order_list =   list(itertools.combinations(type_list,2))


markers  = [ 'CD20', 'CD4', 'CD8', 'FOXP3', 'CD68', 'CD11c', 'CD56', 'CK' ]
pvalues_dict = {}
marker_comp_dict = {}
for marker in markers:
    comp_dict = {}
    for casetype in ['Pleural', 'Peritoneal']:
        for subtype in subtype_dict[casetype]:
            cores_sub = df_map.loc[(df_map['CaseType']==casetype) & (df_map['subtype']==subtype), 'Annotation ID' ]
            df_percent_sub = df_percent.loc[df_percent.index.isin(cores_sub),:]
            comp_dict[f"{casetype}-{subtype}"] = pd.Series(df_percent_sub[marker].dropna()) # use Series help to create dataframe from dictionary of different length
    
    marker_comp_dict[marker] = comp_dict   
    pvalues = []    
    for pair in pair_order_list:
        x = comp_dict[pair[0]]
        y = comp_dict[pair[1]]
        pvalue = stats.ranksums(x, y).pvalue
        pvalues.append(pvalue)
    
    pvalues_dict[marker] = pvalues
 
#--- plot pvalue for sig pairs

markers_sel = []  # ['CD4', 'CD8', 'CD68', 'CK']
for marker in markers:
    if any(ele < 0.05 for ele in pvalues_dict[marker]):
        markers_sel.append(marker)
 
yloc = {}
yloc['CD4'] = 1.05
yloc['CD8'] = 1.5
yloc['CD68'] = 1.35
yloc ['CK'] = 1.35
for marker in markers_sel:  

    df_comp = pd.DataFrame.from_dict(marker_comp_dict[marker])
    
    df_comp = df_comp.reset_index()  # convert index to 'index' columns
    df_comp = df_comp.melt(id_vars=['index'])
    
    pairs = [pair_order_list[i] for i, x in  enumerate(pvalues_dict[marker]) if x < 0.05 ]
    pvalues = [x for i, x in  enumerate(pvalues_dict[marker]) if x < 0.05 ]
    
  
    x="variable"
    y="value"
    ylim= [0,1]
    plot_with_anno (marker, df_comp, x, y, pairs, pvalues, ylim, yloc[marker])
    
    
#----plot pvalue of anova


for marker in markers:
    comp_dict = marker_comp_dict[marker]
    res = stats.kruskal(*list(comp_dict.values())) # with* python pass each list in the list, without *, it only pass the whole list of list
    if res.pvalue > 0.05:
        print (marker, res.pvalue)
        continue
    
    
    pvalues = [res.pvalue]
    
    
    df_comp = pd.DataFrame.from_dict(marker_comp_dict[marker])
    
    df_comp = df_comp.reset_index()  # convert index to 'index' columns
    df_comp = df_comp.melt(id_vars=['index'])
    
    pairs = [(type_list [0], type_list[len(type_list )-1])]
    
    x="variable"
    y="value"
    ylim = [0,1]
    plot_with_anno (marker, df_comp, x, y, pairs, pvalues, ylim)
    
    
######################################################################
# celltype composition difference between casetype

data1 = df.loc[ (df['Classification']=='Malignant') ,:]
data2 = df2.loc[ (df2['Classification']=='Malignant'),:]

by = "MVB"
fc_p1 = data1.groupby([by]).apply(lambda x: x['phenotype_combined'].value_counts(sort=False, normalize=True).rename_axis('Type').reset_index(name='Percent'))
fc_p2 = data2.groupby([by]).apply(lambda x: x['phenotype_combined'].value_counts(sort=False, normalize=True).rename_axis('Type').reset_index(name='Percent'))


df_fc1 = fc_p1.reset_index()
df_fc2 = fc_p2.reset_index()

df_percent1 = df_fc1.set_index([by, 'Type'])['Percent'].unstack().reset_index()
df_percent2 = df_fc2.set_index([by, 'Type'])['Percent'].unstack().reset_index() #colnames ['MVB', 'CD11c', 'CD56', 'other']
df_percent2.drop(['other'], axis=1, inplace=True) # delete column other

# combine df_percent and df_percent2, deduct CD56 and CD11c percents from other in df_percent
df_percent = pd.merge(df_percent1, df_percent2, on=by)
df_percent['other'] = df_percent['other'] - df_percent['CD56'] - df_percent['CD11c']
# remove MVB that other < 0
df_percent = df_percent.loc[df_percent['other'] >= 0, :]

df_percent.index = df_percent[by]
df_percent.drop([by], axis=1, inplace=True)
df_percent.index.name= None
df_percent.columns.name = None



markers  = [ 'CD20', 'CD4', 'CD8', 'FOXP3', 'CD68', 'CD11c', 'CD56', 'CK' ]
pvalues_dict = {}
marker_comp_dict = {}
for marker in markers:
    comp_dict = {}
    for casetype in ['Pleural', 'Peritoneal']:
            cores_sub = df_map.loc[(df_map['CaseType']==casetype) , 'MVB' ]
            df_percent_sub = df_percent.loc[df_percent.index.isin(cores_sub),:]
            comp_dict[f"{casetype}"] = pd.Series(df_percent_sub[marker].dropna()) # use Series help to create dataframe from dictionary of different length
    
    marker_comp_dict[marker] = comp_dict   
    pvalues = []    
   
    x = comp_dict['Pleural']
    y = comp_dict['Peritoneal']
    pvalue = stats.ranksums(x, y).pvalue
    
    pvalues_dict[marker] = pvalue
    
markers_sel = []  # ['CD4', 'CD8', 'CD68', 'CK']
for marker in markers:
    if pvalues_dict[marker]<0.05:
        markers_sel.append(marker)
 
 # only have CK with pavelu = 0.02788
