# IMPPAT 2.0
 
This repository contains the codes associated with the following manuscript:<br>
<i>IMPPAT 2.0: an enhanced and expanded phytochemical atlas of Indian medicinal plants.</i><br>
R. P. Vivek-Ananth, Karthikeyan Mohanraj, Ajaya Kumar Sahoo and Areejit Samal*, Manuscript in submission.

## Brief description of IMPPAT 2.0
![DatabaseOverview](https://github.com/asamallab/IMPPAT2/blob/main/DatabaseOverview.png)

Indian Medicinal Plants, Phytochemistry And Therapeutics 2.0 (IMPPAT 2.0) is a manually curated database which has been constructed via digitalization of information from more than 100 books on traditional Indian medicine, 7000+ published research articles and other existing resources. IMPPAT 2.0 is the largest digital database on phytochemicals of Indian medicinal plants to date, and is a significant enhancement and expansion over IMPPAT 1.0. It can be freely acccessed at https://cb.imsc.res.in/imppat.

IMPPAT 2.0 captures three different types of associations:
(a) Indian medicinal plants - Plant parts- Phytochemical constituents
(b) Indian medicinal plants - Plant parts - Therapeutic uses
(c) Indian medicinal plants - Plant parts - Traditional Indian medicinal formulations

The current version 2.0 released on June 1, 2022 of the IMPPAT database captures 4010 Indian medicinal plants, 17967 phytochemicals, 1095 therapeutic uses and 1133 traditional Indian medicinal formulations. This is more than 2-fold increase in Indian medicinal plants coverage and nearly 2-fold increase in the number of phytochemical constituents compared to IMPPAT 1.0. More importantly, the phytochemical constituents and the therapeutic uses of the Indian medicinal plants are now provided at the level of plant part, while IMPPAT 1.0 did not capture the associations at the plant part level. Another significant update is the 5-fold increase in the number of Indian medicinal plants - Phytochemical constituents and Indian medicinal plants - Therapeutic use associations in comparision to IMPPAT 1.0. Note that the earlier version 1.0 of IMPPAT was released on January 25, 2018.

Importantly, for the 17967 phytochemicals in this database, we provide the two-dimensional (2D) and three-dimensional (3D) chemical structures, and have employed cheminformatics tools to compute their physicochemical properties, drug-likeness based on multiple scoring schemes and predicted Absorption, distribution, metabolism, excretion and toxicity (ADMET) properties. Apart from viewing phytochemicals based on physicochemical and drug-likeness properties, and chemical similarity search, the website now enables viewing the phytochemicals at the level of molecular scaffolds.

In a nutshell, IMPPAT is the largest database on phytochemicals of Indian medicinal plants to date, and this resource is a culmination of our ongoing efforts to digitize the wealth of information contained within traditional Indian medicine. IMPPAT provides an integrated platform to apply cheminformatic approaches to accelerate natural product based drug discovery. IMPPAT is also expected to enable application of system-level approaches towards future elucidation of mechanistic links between phytochemicals of Indian medicinal plants and their therapeutic action.

## Code details

The following python scripts can be used to analyze the chemical structures:
1) ChemicalSimilarityNetwork.py : Calculate Tanimoto coefficient for quantifying chemical structure similarity between molecules in a library
2) ChemicalStructureImages.py : Create SVG or PNG images for chemical structures in a library
3) DruglikenessProperties.py : Evaluate drug-likeness of chemical structures in a library
4) MolecularProperties.py : Compute physicochemical properties for chemical structures in a library 
5) MolecularScaffolds.py : Compute molecular scaffolds for chemical structures in a library
6) MurckoScaffold.py : Edited MurckoScaffold.py code of RDKit package to compute Scaffold at Graph/Node level

The above scripts have been provided with information on required packages and input files for their execution.

## Citation
In case you use the codes herein, please cite the manuscript:<br/>
<i>IMPPAT 2.0: an enhanced and expanded phytochemical atlas of Indian medicinal plants</i>, R. P. Vivek-Ananth, Karthikeyan Mohanraj, Ajaya Kumar Sahoo and Areejit Samal*. Manuscript in submission.
