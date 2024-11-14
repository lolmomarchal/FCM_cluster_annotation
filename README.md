# FCM_cluster_annotation

This script is meant for the detection and visualization of heterogenous and homoegenous expression clusters within Flow Cytometry data. The original purpose of this algorithm was the identification of leukemic clusters that were both subtype and patient specific. In a way, this is a form of automated gating for FCM data. It was primarely developed for CLL and B-ALL leukemias, but can further edited for other leukemia types on request basis ( contact: aolmomarchal@ucsd.edu). The program takes the input of Flow Cytometry data that has been converted into a csv/txt file format. Additionally, it makes use of the FlowSOM algorithm, with is cited below (some edits had to be made due to errors in the initial implementation). 

<img src="https://github.com/user-attachments/assets/a6d87bea-4f34-45a1-a147-11a0ddfbaa6f" width="200" /> 
  

# Installation:
1. pip install git+https://github.com/saeyslab/FlowSOM_Python
2. git clone https://github.com/lolmomarchal/FCM_cluster_annotation.git
3. Use conda environment

# Running

--significance_value (default: 0.05): Sets the p-value threshold for statistical significance.
--input_directory: Directory containing the FCM data files to be processed.
--metadata_csv: Path to the CSV file with metadata associated with the FCM data.
--percentile (default: 0.9): Sets the percentile threshold for filtering data points.
--standard_dev_threshold (default: 3): Threshold for filtering based on standard deviations.
--neighbors (default: 30): Number of neighbors to consider in UMAP analysis.
--min_dist (default: 0.5): Minimum distance parameter for UMAP clustering.
--leukemia_type (default: "CLL"): Specifies the leukemia type (e.g., CLL, ALL) for customized processing.
--percent_sampling (default: 1): Percentage of cells to sample for analysis (e.g., 0.5 for 50%).
--times_sampled (default: 1): Number of times to resample the data.
--output_directory: Directory where processed data and results will be saved.
--processes (default: os.cpu_count()): Number of CPU cores to use for parallel processing.


python FCM_preprocessing.py --input_directory /data/fcm_input --metadata_csv /data/metadata.csv --output_directory /data/fcm_output --significance_value 0.01 --plot_UMAP --percentile 0.95 --neighbors 15 --min_dist 0.3 --leukemia_type ALL --processes 4

FlowSOM was developed by:

A. Couckuyt, B. Rombaut, Y. Saeys, and S. Van Gassen, “Efficient cytometry analysis with FlowSOM in Python boosts interoperability with other single-cell tools,” Bioinformatics, vol. 40, no. 4, p. btae179, Apr. 2024, doi: 10.1093/bioinformatics/btae179.

S. Van Gassen et al., “FlowSOM: Using self-organizing maps for visualization and interpretation of cytometry data,” Cytometry Part A, vol. 87, no. 7, pp. 636–645, 2015, doi: 10.1002/cyto.a.22625.
