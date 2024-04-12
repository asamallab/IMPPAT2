# IMPPAT 2.0
 
This repository contains the codes for cheminformatics analysis performed in the following research article:<br>

R. P. Vivek-Ananth, Karthikeyan Mohanraj, Ajaya Kumar Sahoo, and Areejit Samal*, [<i>IMPPAT 2.0: an enhanced and expanded phytochemical atlas of Indian medicinal plants</i>](https://pubs.acs.org/doi/full/10.1021/acsomega.3c00156), ACS Omega, 8(9):8827–8845, 2023.<br>
(* Corresponding author)

## Brief description of IMPPAT 2.0
![DatabaseOverview](https://github.com/asamallab/IMPPAT2/blob/main/DatabaseOverviewUpdated.png)

Indian Medicinal Plants, Phytochemistry And Therapeutics 2.0 (IMPPAT 2.0) is a manually curated database that has been constructed via digitalization of information from more than 100 books on traditional Indian medicine, 7000+ published research articles, and other existing resources. IMPPAT 2.0 is a significant enhancement and expansion over IMPPAT 1.0. It can be accessed at: https://cb.imsc.res.in/imppat.

IMPPAT 2.0 captures the following associations:<br>
(a) Indian medicinal plant - Plant part- Phytochemical <br>
(b) Indian medicinal plant - Plant part - Therapeutic use<br>

The current version 2.0 released on June 17, 2022, of the IMPPAT database, captures 4010 Indian medicinal plants, 17967 phytochemicals, and 1095 therapeutic uses. This is a more than 2-fold increase in the coverage of Indian medicinal plants and a nearly 2-fold increase in the number of phytochemicals compared to IMPPAT 1.0. More importantly, the phytochemicals and the therapeutic uses of the Indian medicinal plants are now provided at the level of plant part. Another significant update is the 5-fold increase in the number of Indian medicinal plant-Phytochemical associations and Indian medicinal plant-Therapeutic use associations in comparison to IMPPAT 1.0. Note that the earlier version 1.0 of IMPPAT was released on January 25, 2018.

Importantly, for the 17967 phytochemicals in this database, we provide the two-dimensional (2D) and three-dimensional (3D) chemical structures, and employed cheminformatics tools to compute their physicochemical properties, drug-likeness based on multiple scoring schemes, and predicted Absorption, distribution, metabolism, excretion and toxicity (ADMET) properties. Apart from viewing phytochemicals based on physicochemical properties, drug-likeness properties, and chemical similarity, the website now enables viewing the phytochemicals based on molecular scaffolds.

In a nutshell, IMPPAT is a manually curated comprehensive database on phytochemicals of Indian medicinal plants to date, and this resource is a culmination of our ongoing efforts to digitize the wealth of information contained within traditional Indian medicine. 

## Code details

The following Python scripts can be used to analyze the chemical structures:
1) [ChemicalSimilarityNetwork.py](https://github.com/asamallab/IMPPAT2/blob/main/CODES/ChemicalSimilarityNetwork.py): Calculate Tanimoto coefficient for quantifying chemical structure similarity between molecules in a library
2) [ChemicalStructureImages.py](https://github.com/asamallab/IMPPAT2/blob/main/CODES/ChemicalStructureImages.py): Create SVG or PNG images for chemical structures in a library
3) [DruglikenessProperties.py](https://github.com/asamallab/IMPPAT2/blob/main/CODES/DruglikenessProperties.py) : Evaluate drug-likeness of chemical structures in a library
4) [MolecularProperties.py](https://github.com/asamallab/IMPPAT2/blob/main/CODES/MolecularProperties.py): Compute physicochemical properties for chemical structures in a library 
5) [MolecularScaffolds.py](https://github.com/asamallab/IMPPAT2/blob/main/CODES/MolecularScaffolds.py): Compute molecular scaffolds for chemical structures in a library
6) [MurckoScaffold.py](https://github.com/asamallab/IMPPAT2/blob/main/CODES/MurckoScaffold.py): Edited MurckoScaffold.py code of RDKit package to compute Scaffold at Graph/Node level

The above scripts have been provided with information on required packages and input files for their execution.

## Citation
In case you use the codes herein, please cite the research article:<br>

R. P. Vivek-Ananth, Karthikeyan Mohanraj, Ajaya Kumar Sahoo, and Areejit Samal*, [<i>IMPPAT 2.0: an enhanced and expanded phytochemical atlas of Indian medicinal plants</i>](https://pubs.acs.org/doi/full/10.1021/acsomega.3c00156), ACS Omega, 8(9):8827–8845, 2023.
