# FCM_cluster_annotation

This script is meant for the detection and visualization of heterogenous and homoegenous expression clusters within Flow Cytometry data. The original purpose of this algorithm was the identification of leukemic clusters that were both subtype and patient specific. In a way, this is a form of automated gating for FCM data. It was primarely developed for CLL and B-ALL leukemias, but can further edited for other leukemia types on request basis ( contact: aolmomarchal@ucsd.edu). The program takes the input of Flow Cytometry data that has been converted into a csv/txt file format. Additionally, it makes use of the FlowSOM algorithm, with is cited below (some edits had to be made due to errors in the initial implementation). 

<img src="https://github.com/user-attachments/assets/a6d87bea-4f34-45a1-a147-11a0ddfbaa6f" width="200" /> 
  

To run:


FlowSOM was developed by:

A. Couckuyt, B. Rombaut, Y. Saeys, and S. Van Gassen, “Efficient cytometry analysis with FlowSOM in Python boosts interoperability with other single-cell tools,” Bioinformatics, vol. 40, no. 4, p. btae179, Apr. 2024, doi: 10.1093/bioinformatics/btae179.

S. Van Gassen et al., “FlowSOM: Using self-organizing maps for visualization and interpretation of cytometry data,” Cytometry Part A, vol. 87, no. 7, pp. 636–645, 2015, doi: 10.1002/cyto.a.22625.
