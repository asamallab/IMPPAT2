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

The current version 2.0 released on June 1, 2022 of the IMPPAT database captures 4010 Indian medicinal plants, 17967 phytochemicals, 1095 therapeutic uses and 1133 traditional Indian medicinal formulations. This is more than 2-fold increase in Indian medicinal plants coverage and nearly 2-fold increase in the number of phytochemical constituents compared to IMPPAT 1.0. More importantly, the phytochemical constituents and the therapeutic uses of the Indian medicinal plants are now provided at the level of plant part, while IMPPAT 1.0 did not capture the associations at the plant part level. Another significant update is the 5-fold increase in the number of Indian medicinal plants - Phytochemical constituents and Indian medicinal plants - Therapeutic use associations in comparision to IMPPAT 1.0. Note that the earlier version 1.0 of IMPPAT ws released on January 25, 2018.

Importantly, for the 17967 phytochemicals in this database, we provide the two-dimensional (2D) and three-dimensional (3D) chemical structures, and have employed cheminformatics tools to compute their physicochemical properties, drug-likeliness based on multiple scoring schemes and predicted Absorption, distribution, metabolism, excretion and toxicity (ADMET) properties. Apart from viewing phytochemicals based on physicochemical and drug-likeness properties, and chemical similarity search, the website now enables viewing the phytochemicals at the level of molecular scaffolds.

In a nutshell, IMPPAT is the largest database on phytochemicals of Indian medicinal plants to date, and this resource is a culmination of our ongoing efforts to digitize the wealth of information contained within traditional Indian medicine. IMPPAT provides an integrated platform to apply cheminformatic approaches to accelerate natural product based drug discovery. IMPPAT is also expected to enable application of system-level approaches towards future elucidation of mechanistic links between phytochemicals of Indian medicinal plants and their therapeutic action.

Code/Script Details:
=================================================

The following two scripts can be used to filter the cross-correlation matrices and generate edge files and node files of the filtered networks:
1) mst_wt.py : Python script to generate a weighted or unweighted filtered minimum spanning tree + thresholded network from the weighted network of cross-correlation values. The weights are interpreted as distances (costs).
2) PMFG_wt.py: Python script to generate a weighted PMFG from the weighted network of cross-correlation values. The weights are interpreted as distances (costs).

The following scripts can be used to compute the different network measures for the filtered networks:
1) clique_number.py : Clique number
2) diameter_wt.py : Diameter of a weighted network
3) eigenvector_centality.py : Eigenvector centrality for all the nodes of a weighted network
4) FormanUndirected.cpp : Forman-Ricci curvature for all the edges of a weighted/unweighted network 
5) ga_wt: Global assortativity of a weighted network
6) graph_measures.py : Number of edges, Average degree, Average Weighted Degree, Edge Density, Average Clustering coefficient
7) grc_wt_undir.py : Global Reaching Centrality of a weighted network
8) network_entropy.py : Entropy of an unweighted network
9) comm_eff.py : Communication efficiency of a weighted network
10) OR-UnDir.py : Ollivier-Ricci curvature for all the edges of a weighted/unweighted but undirected network
11) MengerHaantjesUnweighted.py : Menger-curvature and Haantjes-curvature for all the edges of an unweighted network

### CITATION
In case you use the codes or data herein, please cite the manuscript:<br/>
<i>IMPPAT 2.0: an enhanced and expanded phytochemical atlas of Indian medicinal plants</i>, R. P. Vivek-Ananth, Karthikeyan Mohanraj, Ajaya Kumar Sahoo and Areejit Samal*. Manuscript in submission.

We welcome suggestions from researchers in academia and industry for future improvements of our database.
