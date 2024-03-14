
<img src="https://github.com/osmanbeyoglulab/MesotheliomaSpatialAtlas_analyses/blob/main/assets/diagram.jpg" alt="drawing" width="600"/>


## Introduction

The repository contains the data analyses code for the paper paper: Ma X, Lembersky D, Kim ES, Bruno TC, Testa JR, Osmanbeyoglu HU. [Spatial landscape of malignant pleural and peritoneal mesothelioma tumor immune microenvironment](https://www.biorxiv.org/content/10.1101/2021.09.07.459263v3.full). bioRxiv 2021:2021.2009.2007.459263.

We conducted multiplex immunofluorescence (mIF) analyses on tissue microarrays (n=3) from malignant peritoneal (MPeM, n=23) and pleural (MPM, n=79) mesothelioma patients. Our study aimed to elucidate spatial distributions of key immune cell populations and their association with LAG3, BAP1, NF2, and MTAP, with MTAP serving as a CDKN2A/B surrogate marker. Additionally, we examined the relationship between the spatial distribution of major immune cell types with MM patient prognosis and clinical characteristics. We observed a higher degree of interaction between immune cells and tumor cells in MPM compared to MPeM. Notably, within MPM tumors, we detected a significantly increased interaction between tumor cells and CD8+ T cells in tumors with low BAP1 expression compared to those with high BAP1 expression


## TMA data

The TMA_data from Akoya software which is the inital input of our data analysis pipeline. The unrecognized cells by Akoya have been removed. Data for each panel and patient clinic information can be downloaded from https://sites.pitt.edu/~xim33/Mesothelioma

## Intermidiate data
Some itermidiate data generated in the pipe line is under /TMA_data folder in the repository 

## Analysis code

Out analysis pipeline involve both Python and R prgramming. The codes are located in the root folder. We also integrated relavent code to generate the manuscript figure in the /notebook folder

## Citation
If you find the data or code from this repository helpful, please cite this paper:
```
@inproceedings{Ma2023git,
  title = {Spatial landscape of malignant pleural and peritoneal mesothelioma tumor immune microenvironment},
  author = {Xiaojun Ma and 
            David Lemberskya and 
            Elena S Kimaf and 
            Tullia C Brunoa and 
            Michael J. Becicha and 
            Joseph R. Testad and 
            Hatice U Osmanbeyoglua},
  url = {https://www.biorxiv.org/content/10.1101/2023.09.06.556559v1},
  year = {2023},
}
