# Standard Library imports

import gc
import os
import pickle
import random
import multiprocessing as mp
import logging
import datetime
import tqdm
import argparse

# Data Science Library imports
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import umap.umap_ as umap
import pandas as pd
import numpy as np
from scipy.stats import mannwhitneyu

import FlowSOM_Python.src.FlowSOM as fs

# Constants

global UMAP_files
global d05n30
global d04n30
global d03n30
global d02n30
global leukemia_type
global count
global cell_proportion
global prepro_path_

# logging
logging.disable(logging.DEBUG)
logging.disable(logging.ERROR)
logging.disable(logging.WARNING)


# ok first thing need to unite all patient info

def upload_patient_csvs(directory, df):
    # open dataframe with all sample info and make sure to put directory in file path
    now = datetime.datetime.now()
    print(f"Uploading all Patient information {now.time()}")
    file_data = pd.read_csv(df)
    file_data["File"] = file_data["File"].apply(lambda x: os.path.join(directory, x))
    output_path = os.path.join(directory, "appended.csv")
    patient_index = os.path.join(directory, "patient_index.pickle")
    patient_dict = {}
    # need to append all patient data together
    if not os.path.isfile(os.path.join(directory, "appended.csv")):
        index = 0
        all_dfs = []
        for i, row in tqdm.tqdm(file_data.iterrows(), total=file_data.shape[0]):
            label = row["Label"]
            file_name = row["File"]
            temp = pd.read_csv(file_name, sep="\t")
            if "Time" in temp.columns:
                temp = temp.drop(columns=["Time"])
            patient_dict[row["patient_id"]] = [index, index + len(temp)]  # save index of when patient starts/ends
            index += len(temp)

            all_dfs.append(temp)
        df = pd.concat(all_dfs, ignore_index=True)
        df.to_csv(output_path, index=False)
        with open(patient_index, "wb") as handle:
            pickle.dump(patient_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)
    else:
        with open(patient_index, "rb") as f:
            patient_dict = pickle.load(f)
    now = datetime.datetime.now()
    print(f"Done with appending all patients {now.time()}")
    return output_path, file_data, patient_dict


def three_std_away(positive, negative, standard_dev_threshold=3):
    mean = np.mean(negative)
    std_dev = np.std(negative)
    return (positive - mean) > standard_dev_threshold * std_dev


def assigningClusters(clusters, patient_dict, df, significance_threshold=0.05, percentile=90, standard_dev_threshold=3):
    labels = pd.read_csv(df)
    # initiate cluster count
    positive_count = 0
    negative_count = 0
    positive_cluster = {}
    negative_cluster = {}
    for i in range(0, 40):
        positive_cluster[i] = []
        negative_cluster[i] = []
    # get patient level proportions for clusters
    for key in patient_dict.keys():
        start = patient_dict[key][0]
        end = patient_dict[key][1]
        patient_cells = clusters.iloc[start:end]
        label = labels[labels["patient_id"] == key].iloc[0]["Label"]
        if label == 1:
            for cluster in positive_cluster.keys():
                positive_cluster[cluster].append(
                    len(patient_cells[patient_cells["metaclusters"] == cluster]) / len(patient_cells))
        else:
            for cluster in positive_cluster.keys():
                negative_cluster[cluster].append(
                    len(patient_cells[patient_cells["metaclusters"] == cluster]) / len(patient_cells))
        # getting homogenous clusters -> want to get most siginificant clusters
    t_values = {}
    for i in range(0, 40):
        t_values[i] = float("inf")
    for key in t_values.keys():
        t_statistic, p_value = mannwhitneyu(positive_cluster[key], negative_cluster[key])
        t_values[key] = p_value
    # filter out for significant clusters
    significant_items = {k: v for k, v in t_values.items() if v < significance_threshold}
    if significant_items:
        if leukemia_type == "CLL":
            # only want the top 3 clusters(the ones with most significant p values)
            largest_items = sorted(significant_items.items(), key=lambda x: x[1])[:3]
            significant_clusters = [item[0] for item in largest_items]
            clusters["cell_label"] = clusters.apply(
                lambda row: 1 if row["metaclusters"] in significant_clusters else 0, axis=1)

        else:
            # want all significant clusters
            largest_items = sorted(significant_items.items(), key=lambda x: x[1])
            significant_clusters = [item[0] for item in largest_items]
            clusters["cell_label"] = clusters.apply(
                lambda row: 1 if row["metaclusters"] in significant_clusters else 0, axis=1)
    else:
        print(f"There were 0 homogenous clusters found with a significance threshold of {significance_threshold}")
        clusters["cell_label"] = 0

    # print(f"homogenous:clusters {significant_clusters}")
    # second step get the patient-specific cluster
    for key in patient_dict.keys():
        cluster_patient = {}
        for k in range(0, 40):
            cluster_patient[k] = 0
        # get patient cells
        start = patient_dict[key][0]
        end = patient_dict[key][1]
        patient_cells = clusters.iloc[start:end]
        label = labels[labels["patient_id"] == key].iloc[0]["Label"]
        if label == 1:
            # get proportion of each cluster for the patient
            for cluster in cluster_patient.keys():
                cluster_patient[cluster] = len(patient_cells[patient_cells["metaclusters"] == cluster]) / len(
                    patient_cells)
            # now we are going to check the two conditions: 3 standard deviations away from mean of out of class and within the 90th percentile of expression in the in class
            patient_specific_clusters = []
            for cluster in cluster_patient.keys():
                threshold_upper = np.percentile(np.array(positive_cluster[cluster]), percentile)
                is_3std = three_std_away(cluster_patient[cluster], negative_cluster[cluster],
                                         standard_dev_threshold=standard_dev_threshold)
                if is_3std and cluster_patient[cluster] >= threshold_upper:
                    patient_specific_clusters.append(cluster)
            patient_specific_clusters = np.concatenate(
                (np.array(significant_clusters), np.setdiff1d(patient_specific_clusters, significant_clusters)))
            clusters.loc[start:end - 1, "cell_label"] = clusters.loc[start:end - 1].apply(
                lambda row: 1 if row["metaclusters"] in patient_specific_clusters else 0, axis=1)


def metaclustering(directory, appended_csv, patient_dict, df, significance_threshold=0.05, percentile=90,
                   standard_dev_threshold=3):
    skip = True
    if not skip:
        labels = pd.read_csv(df)
        now = datetime.datetime.now()
        print(f"clustering patient data {now.time()}")
        cluster_path = os.path.join(directory, "ClusterInformation.csv")
        if not os.path.isfile(cluster_path):
            clusters = pd.read_csv(appended_csv)
            ff = fs.io.read_csv(appended_csv)
            ff.uns['meta'] = {}
            if leukemia_type == "CLL":
                cols = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]
            else:
                cols = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]

            fsom = fs.FlowSOM(ff, cols, xdim=20, ydim=20, n_clus=40)
            clusters["metaclusters"] = np.array(fsom.get_cell_data().obs["metaclustering"])
            clusters.to_csv(os.path.join(directory, "ClusterInformation.csv"), index=False)
        else:
            clusters = pd.read_csv(cluster_path)

            now = datetime.datetime.now()
            print(f"Done with Clustering, analyzing clusters {now.time()}")
            assigningClusters(clusters, patient_dict, df,
                              significance_threshold=significance_threshold, percentile=percentile,
                              standard_dev_threshold=standard_dev_threshold)
            now = datetime.datetime.now()
            print(f"Finishined analyzing clusters, saving patient single cell labels {now.time()}")
            for key in tqdm.tqdm(patient_dict.keys(), total=len(patient_dict.keys())):
                start = patient_dict[key][0]
                end = patient_dict[key][1]
                patient_cells = clusters.iloc[start:end]
                patient_path = labels[labels["patient_id"] == key].iloc[0]["File"]
                patient_cells.to_csv(os.path.join(directory, patient_path), index=False)
                del patient_cells
                now = datetime.datetime.now()
            print(f"Finishined assigning single cell labels {now.time()}")


def randomized_cells(path, proportion):
    """
  randomized_cells, reads the data from a patient CSV file
   and then randomly samples rows.
  :param path: path of sample
  :param proportion: proportion of cells we are sampling
  :return: sampled cell/marker dataframe from patient
  """
    if os.path.exists(path):
        # File exists, so proceed to read or process it
        data = pd.read_csv(path)
        if "Time" in data.columns:
            data = data.drop(columns="Time")
        if proportion == 1:
            return data
        random_state = random.randint(1, 900000)  # Generate a random integer between 1 and 100000
        random_data = data.sample(frac=proportion, random_state=random_state)
        random_data = random_data.reset_index(drop=True)
        del data
        return random_data
    else:
        # File does not exist, so skip processing
        print("File does not exist. Skipping. " + path)
        return pd.DataFrame()


def UMAP_initialize(path, template_path, min_dist, neighbors):
    """
    :param template_path: path to template
 :param path: of folder to save UMAP models
 :param proportion: to get appropriate of cells to be used
 :return: UMAP model paths
 """
    model_filename = os.path.join(path, f"umap_d{min_dist}_n{neighbors}.pkl")
    if not os.path.isfile(model_filename):
        template = pd.read_csv(template_path)
        model = umap.UMAP(n_neighbors=neighbors, min_dist=min_dist).fit(template)
        with open(model_filename, "wb") as file:
            pickle.dump(model, file)
        del template
        gc.collect()
    else:
        with open(model_filename, "rb") as file:
            model = pickle.load(file)
    return model


def UMAP_generation(neighbors, distance, data, path, label, patient_id, umap_model, run_count, type):
    """
    :param neighbors: amount of neighbors used
    :param distance: distance used
    :param data: what data is being transformed
    :param path: to save outputs
    :param label: to determine where to save outputs (Classification folders)
    :param patient_id: For naming purposes (i.e. sample 5567)
    :param umap_model: to see what pre-made model is being used.
    :param count: Used for iterative runs, for what run the output belongs to
    :param type: Leukemia type, affects the type of 2d plots used
    """

    matplotlib.use('Agg')
    # Creating directory name and filenames based on label
    if type == "CLL":
        if label == 1:
            directory = os.path.join(path, "CLL", str(patient_id))
            filename = str(type) + "_" + str(patient_id) + "_d" + str(distance) + "_n" + str(
                neighbors) + "_" + "run_" + str(run_count)
            filenamefcm = "FCMDATA_" + str(type) + "_" + str(patient_id) + "_d" + str(distance) + "_n" + str(
                neighbors) + "_" + "run_" + str(run_count)
        else:
            directory = os.path.join(path, "non-CLL", str(patient_id))
            filename = "NO" + str(type) + "_" + str(patient_id) + "_d" + str(distance) + "_n" + str(
                neighbors) + "_" + "run_" + str(
                run_count)
            filenamefcm = "FCMDATA_NO" + str(type) + "_" + str(patient_id) + "_d" + str(distance) + "_n" + str(
                neighbors) + "_" + "run_" + str(run_count)
    elif type == "B-ALL":
        if label == 1:
            directory = os.path.join(path, "B-ALL", str(patient_id))
            filename = str(type) + "_" + str(patient_id) + "_d" + str(distance) + "_n" + str(
                neighbors) + "_" + "run_" + str(run_count)
            filenamefcm = "FCMDATA_" + str(type) + "_" + str(patient_id) + "_d" + str(distance) + "_n" + str(
                neighbors) + "_" + "run_" + str(run_count)
        else:
            directory = os.path.join(path, "non-B-ALL", str(patient_id))
            filename = "NO" + str(type) + "_" + str(patient_id) + "_d" + str(distance) + "_n" + str(
                neighbors) + "_" + "run_" + str(
                run_count)
            filenamefcm = "FCMDATA_NO" + str(type) + "_" + str(patient_id) + "_d" + str(distance) + "_n" + str(
                neighbors) + "_" + "run_" + str(run_count)

    mask = data["cell_label"] == 0

    if not os.path.isfile(os.path.join(directory, filenamefcm + ".csv")):
        sns.set(context='poster', style='white', rc={'figure.figsize': (30, 30)})

    # UMAP layout generated from template
    trans = umap_model.transform(data.drop(columns=["cell_label", "metaclusters"]))
    template = pd.DataFrame(trans, columns=["umap1", "umap2"])
    data = pd.concat([data, template], axis=1)

    if type == "CLL":
        fig, axes = plt.subplots(3, 3, figsize=(30, 30))
        plots = [
            (0, 0, "umap1", "umap2"),
            (0, 1, "FSC-A", "SSC-A"),
            (0, 2, "SSC-H", "CD45"),
            (1, 0, "CD5", "CD19"),
            (1, 1, "CD10", "CD79b"),
            (1, 2, "CD38", "CD22"),
            (2, 0, "CD43", "CD81"),
            (2, 1, "CD3", "CD19"),
            (2, 2, "CD79b", "CD81")
        ]

        for row, col, x, y in plots:
            # Plot negative cells (background) first
            sns.scatterplot(
                x=x, y=y, data=data[mask], color="blue", ax=axes[row, col], s=1, label = "negative")
            # Overlay positive cells
            sns.scatterplot(
                x=x, y=y, data=data[~mask], color="red", ax=axes[row, col], s=1, label = "positive"
            )

            # Set labels and margins
            axes[row, col].set_xlabel(x)
            axes[row, col].set_ylabel(y)
            axes[row, col].margins(x=0, y=0)

            # Remove legend for all but the first subplot
            if not (row == 0 and col == 0):
                legend = axes[row, col].get_legend()
                if legend is not None:  # Check if legend exists before trying to remove it
                    legend.remove()

            # Adjust layout for spacing
        axes[0, 0].legend(loc='upper right')
        plt.subplots_adjust(hspace=0.5, wspace=0.5)
        plt.tight_layout(rect=[0, 0, 0.85, 1])  # Adjusts for space on right side for legend
        # Save the figure
        os.makedirs(directory, exist_ok=True)
        plt.savefig(os.path.join(directory, filename + "_cancerous_cell_labels.png"), bbox_inches="tight")
        plt.close(fig)
        gc.collect()

        # Plot metaclusters for each subplot

        fig, axes = plt.subplots(3, 3, figsize=(30, 30))
        for row, col, x, y in plots:
            sns.scatterplot(
                x=x, y=y, data=data, hue="metaclusters", ax=axes[row, col], s=1, palette="viridis"
            )
            axes[row, col].set_xlabel(x)
            axes[row, col].set_ylabel(y)
            axes[row, col].margins(x=0, y=0)

            if not (row == 0 and col == 0):
                axes[row, col].get_legend().remove()
        plt.subplots_adjust(hspace=0.5, wspace=0.5)
        plt.tight_layout(rect=[0, 0, 0.85, 1])  # Leave space on the right for legend
        plt.savefig(os.path.join(directory, filename + "_metacluster_cell_labels.png"), bbox_inches="tight")
        plt.close(fig)
        gc.collect()
    if type == "B-ALL":
        fig, axes = plt.subplots(3, 3, figsize=(15, 15))
        plots = [
            (0, 0, "umap1", "umap2"),
            (0, 1, "FSC", "SSC"),
            (0, 2, "CD19", "CD34"),
            (1, 0, "CD66b", "CD22"),
            (1, 1, "CD20", "CD38"),
            (1, 2, "CD20", "CD10"),
            (2, 0, "CD24", "CD45"),
            (2, 1, "CD10", "CD38"),
            (2, 2, "SSC", "CD45")
        ]
        for row, col, x, y in plots:
            sns.scatterplot(x=x, y=y, data=data, hue="cell_label", ax=axes[row, col], s=10)
        # Set labels and margins
        axes[row, col].set_xlabel(x)
        axes[row, col].set_ylabel(y)
        axes[row, col].margins(x=0, y=0)
        if not (row == 0 and col == 0):
            axes[row, col].get_legend().remove()
        os.makedirs(directory, exist_ok=True)
        # saving full-image
        plt.savefig(os.path.join(directory, filename + "_cancerous_cell_labels" + ".png"))
        plt.close(fig)
        gc.collect()
        # get metacluster specific
        for row, col, x, y in plots:
            fig, axes = plt.subplots(3, 3, figsize=(15, 15))
            sns.scatterplot(x=x, y=y, data=data, hue="metaclusters", ax=axes[row, col], s=10)
            # Set labels and margins
            axes[row, col].set_xlabel(x)
            axes[row, col].set_ylabel(y)
            axes[row, col].margins(x=0, y=0)
            if not (row == 0 and col == 0):
                axes[row, col].get_legend().remove()
        plt.savefig(os.path.join(directory, filename + "_metacluster_cell_labels" + ".png"))
        plt.close(fig)
        gc.collect()

    # saving full dataframes
    pd.DataFrame.to_csv(template, os.path.join(directory, "UMAP_" + filename + ".csv"))
    pd.DataFrame.to_csv(data, os.path.join(directory, filenamefcm + ".csv"))

    plt.cla()
    plt.clf()
    plt.close('all')


def process_file(args):
    UMAP_model, args, row = args
    # get sample information
    path = row['File']
    label = row['Label']
    patient_id = row['patient_id']
    count = args.times_sampled
    cell_proportion = args.percent_sampling
    path = os.path.join(args.input_directory, path)
    try:
        if os.path.exists(path):
            for i in range(count):
                data = randomized_cells(path, cell_proportion)
                UMAP_generation(args.neighbors, args.min_dist, data, args.output_directory,
                                label, patient_id, UMAP_model, i, args.leukemia_type)
                del data
        else:
            print("File does not exist. Skipping. " + path)
        gc.collect()
    except Exception as e:
        print(f"Error occurred while processing {path}: {e}")


def template(directory, label_data):
    if not os.path.isfile(os.path.join(directory, "template.csv")):
        # get 10 random patients (5 healthy 5 cancerous)
        healthy_samples = label_data[label_data['Label'] == 0].sample(n=5, random_state=42)
        cancer_samples = label_data[label_data['Label'] == 1].sample(n=5, random_state=42)
        patients = pd.concat([healthy_samples, cancer_samples])
        template = pd.DataFrame()

        for i, row in patients.iterrows():
            path = row["File"]
            data = pd.read_csv(os.path.join(directory, path))
            if "Time" in data.columns:
                data = data.drop(columns="Time")
            data = data.drop(columns=["cell_label", "metaclusters"])
            template = pd.concat([template, data])
        template.to_csv(os.path.join(directory, "template.csv"), index=False)
    now = datetime.datetime.now()
    print(f"Finished Generating Template {now.time()}")
    return os.path.join(directory, "template.csv")


def get_args():
    parser = argparse.ArgumentParser(
        prog='FCM_preprocessing',
        description='Preprocessing FCM data, obtains single-cell labels and represents them in visual form')
    parser.add_argument("--significance_value", default=0.05)
    parser.add_argument("--input_directory")
    parser.add_argument("--metadata_csv")
    parser.add_argument("--percentile", default=0.9)
    parser.add_argument("--standard_dev_threshold", default=3)
    parser.add_argument("--plot_UMAP", action="store_true")
    parser.add_argument("--neighbors", default=30)
    parser.add_argument("--min_dist", default=0.5)
    parser.add_argument("--leukemia_type", default="CLL")
    parser.add_argument("--percent_sampling", default=1)
    parser.add_argument("--times_sampled", default=1)
    parser.add_argument("--output_directory")
    parser.add_argument("--processes", default=os.cpu_count())
    return parser.parse_args()


def main():
    args = get_args()
    directory = args.input_directory
    df = args.metadata_csv
    percentile = args.percentile
    standard_dev_threshold = args.standard_dev_threshold
    plot_UMAP = args.standard_dev_threshold
    neighbors = args.neighbors
    min_dist = args.min_dist

    # set global variables
    gc.enable()
    matplotlib.use('Agg')
    global leukemia_type
    leukemia_type = args.leukemia_type
    global count
    count = args.times_sampled
    global cell_proportion
    cell_proportion = args.percent_sampling
    global prepro_path_
    prepro_path_ = args.output_directory
    max_processes = args.processes

    # make output directory
    os.makedirs(args.output_directory, exist_ok=True)

    # first cluster patients
    output_path, file_data, patient_dict = upload_patient_csvs(directory, df)
    metaclustering(directory, output_path, patient_dict, df, significance_threshold=args.significance_value
                   , percentile=percentile, standard_dev_threshold=standard_dev_threshold)

    # create UMAP template
    template_path = template(directory, pd.read_csv(df))

    # create subcategories for in class vs out of class
    if args.leukemia_type == "CLL":
        os.makedirs(os.path.join(args.output_directory, "CLL"), exist_ok=True)
        os.makedirs(os.path.join(args.output_directory, "non-CLL"), exist_ok=True)
    if args.leukemia_type == "B-ALL":
        os.makedirs(os.path.join(args.output_directory, "B-ALL"), exist_ok=True)
        os.makedirs(os.path.join(args.output_directory, "non-B-ALL"), exist_ok=True)

    # initialize UMAP model

    UMAP_model = UMAP_initialize(args.input_directory, template_path, min_dist=min_dist, neighbors=neighbors)

    now = datetime.datetime.now()
    print(f"Visualizing samples and obtaining UMAP representation: {now.time()}")
    label_data = pd.read_csv(df)

    mult_args = [[UMAP_model, args, row] for _, row in label_data.iterrows()]

    with mp.Pool(processes=max_processes) as pool:
        results = pool.map(process_file, mult_args)

    # After intializing model


if __name__ == "__main__":
    main()
