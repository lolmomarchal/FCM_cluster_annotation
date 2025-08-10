### FCM Cluster Annotation 🔬

This script is for the detection and visualization of **heterogeneous** and **homogeneous** expression clusters in Flow Cytometry (FCM) data. It was originally developed to identify patient-specific and subtype-specific leukemic clusters, serving as a form of automated gating for FCM data. While primarily created for **CLL** and **B-ALL** leukemias, it can be adapted for other leukemia types upon request. Any questions should be directed to lorenzo.olmomarchal@gmail.com

![FCM Cluster Annotation Image](https://github.com/user-attachments/assets/a6d87bea-4f34-45a1-a147-11a0ddfbaa6f)

---

### **Installation** ⚙️

1.  **Clone the repository:**

    ```bash
    git clone [https://github.com/lolmomarchal/FCM_cluster_annotation.git](https://github.com/lolmomarchal/FCM_cluster_annotation.git)
    ```

2.  **Install the FlowSOM dependency:** The script uses the FlowSOM algorithm for its clustering analysis.

    ```bash
    pip install git+[https://github.com/saeyslab/FlowSOM_Python](https://github.com/saeyslab/FlowSOM_Python)
    ```

3.  **Set up the environment:** It's recommended to use a Conda environment to manage dependencies.

    ```bash
    # Create and activate a new conda environment
    conda create -n fcm_env python=3.8
    conda activate fcm_env

    # Install the required packages
    pip install -r requirements.txt
    ```
    (Note: You will need a `requirements.txt` file listing all necessary packages.)

---

### **How to Run** ▶️

You can run the script from the command line using the following parameters:

* `--significance_value` (default: 0.05): Sets the **p-value threshold** for statistical significance.
* `--input_directory`: Path to the directory containing your **FCM data files** in `.csv` or `.txt` format.
* `--metadata_csv`: Path to the CSV file with **metadata** for the FCM data.
* `--percentile` (default: 0.9): Threshold for filtering data points based on percentile.
* `--standard_dev_threshold` (default: 3): Threshold for filtering based on standard deviations.
* `--neighbors` (default: 30): Number of neighbors to consider in **UMAP analysis**.
* `--min_dist` (default: 0.5): Minimum distance parameter for **UMAP clustering**.
* `--leukemia_type` (default: "CLL"): Specifies the leukemia type, such as **CLL** or **ALL**, for specific processing.
* `--percent_sampling` (default: 1): The percentage of cells to sample (e.g., `0.5` for 50%).
* `--times_sampled` (default: 1): The number of times to resample the data.
* `--output_directory`: The directory where processed data and results will be saved.
* `--processes` (default: `os.cpu_count()`): The number of **CPU cores** to use for parallel processing.

#### **Example Command**

```bash
python FCM_preprocessing.py --input_directory /data/fcm_input --metadata_csv /data/metadata.csv --output_directory /data/fcm_output --significance_value 0.01 --plot_UMAP --percentile 0.95 --neighbors 15 --min_dist 0.3 --leukemia_type ALL --processes 4
```


#### Citations 📖
This project uses the FlowSOM algorithm. Please cite the following papers if you use this tool in your work:

A. Couckuyt, B. Rombaut, Y. Saeys, and S. Van Gassen, “Efficient cytometry analysis with FlowSOM in Python boosts interoperability with other single-cell tools,” Bioinformatics, vol. 40, no. 4, p. btae179, Apr. 2024.

S. Van Gassen et al., “FlowSOM: Using self-organizing maps for visualization and interpretation of cytometry data,” Cytometry Part A, vol. 87, no. 7, pp. 636–645, 2015.
