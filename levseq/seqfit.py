"""
A script for visualzing the sequnece-fitness relationship
"""

from typing import Optional, Dict, Union, List

import re
import os
import html
import psutil
import socket
import param

import warnings

from pathlib import Path
from copy import deepcopy
from glob import glob

from ast import literal_eval
import urllib.parse

import numpy as np
import pandas as pd
import biopandas as Bio

import panel as pn
import holoviews as hv

from bokeh.io import output_notebook, show
from bokeh.plotting import figure

import ninetysix as ns

try:
    from rdkit import Chem

    import torch
    import torch.nn as nn

    from chai_lab.chai1 import run_inference

    gen_struct = True

except:
    pass

    gen_struct = False


# Get them w.r.t to a mutation
from scipy.stats import mannwhitneyu
from tqdm import tqdm
import pandas as pd
import numpy as np
from collections import defaultdict

# Enable Bokeh to display plots in the notebook
hv.extension("bokeh")
pn.extension()
output_notebook()

warnings.simplefilter(action="ignore", category=FutureWarning)

# Amino acid code conversion
AA_DICT = {
    "Ala": "A",
    "Cys": "C",
    "Asp": "D",
    "Glu": "E",
    "Phe": "F",
    "Gly": "G",
    "His": "H",
    "Ile": "I",
    "Lys": "K",
    "Leu": "L",
    "Met": "M",
    "Asn": "N",
    "Pro": "P",
    "Gln": "Q",
    "Arg": "R",
    "Ser": "S",
    "Thr": "T",
    "Val": "V",
    "Trp": "W",
    "Tyr": "Y",
    "Ter": "*",
}

def calculate_mutation_combinations(stats_df):
    mutation_dict = defaultdict(list)
    for mutation in stats_df['mutation'].values:
        mutations = mutation.split('_')
        for m in mutations:
            mutation_dict[m].append(mutation)

    rows = []
    with pd.ExcelWriter('mutations.xlsx', engine='xlsxwriter') as writer:
        for mutation, mutations in mutation_dict.items():
            # Here we want to now get the values for each of these i.e. the stats values for each one and summarize it maybe for now we'll just make a excel file
            df1 = stats_df[stats_df['mutation'].isin(mutations)]
            mutation = mutation.replace('*', '.')
            df1.to_excel(writer, sheet_name=mutation)
            # Also just take the mean of the mean lol and the sum of the number of the wells
            rows.append([mutation, np.sum(df1['number of wells with mutation'].values), '|'.join(set(list(mutations))),
                         np.mean(df1['mean'].values),
                         np.median(df1['median'].values), np.mean(df1['amount greater than parent mean'].values),
                         np.max(df1['amount greater than parent mean'].values)])

    df = pd.DataFrame(rows, columns=['mutation', 'number of wells with mutation',
                                     'other-mutations', 'mean', 'median',
                                     'mean amount greater than parent', 'max amount greater than parent'])
    df.sort_values(by='mean amount greater than parent', ascending=False)
    return df


def normalise_calculate_stats(processed_plate_df, value_columns, normalise='standard', stats_method='mannwhitneyu',
                              parent_label='#PARENT#', normalise_method='median'):
    parent = parent_label
    # if nomrliase normalize with standard normalisation
    normalised_value_columns = []
    normalised_df = pd.DataFrame()
    if normalise:
        for plate in set(processed_plate_df['Plate'].values):
            for value_column in value_columns:
                sub_df = processed_plate_df[processed_plate_df['Plate'] == plate]
                parent_values = sub_df[sub_df['amino-acid_substitutions'] == parent][value_column].values
                # By default use the median
                if normalise_method == 'median':
                    parent_mean = np.median(parent_values)
                else:
                    parent_mean = np.mean(parent_values)
                parent_sd = np.std(parent_values)

                # For each plate we normalise to the parent of that plate
                sub_df[f'{value_column} plate standard norm'] = (sub_df[value_column].values - parent_mean) / parent_sd
                normalised_value_columns.append(f'{value_column} plate standard norm')
                normalised_df = pd.concat([normalised_df, sub_df])
    else:
        normalised_df = processed_plate_df

    normalised_value_columns = list(set(normalised_value_columns))
    processed_plate_df = normalised_df

    sd_cutoff = 1.5  # The number of standard deviations we want above the parent values
    # Now for all the other mutations we want to look if they are significant, first we'll look at combinations and then individually
    grouped_by_mutations = processed_plate_df.groupby('amino-acid_substitutions')

    rows = []
    for mutation, grp in tqdm(grouped_by_mutations):
        # Get the values and then do a ranksum test
        if mutation != parent:
            for value_column in normalised_value_columns:
                parent_values = list(processed_plate_df[processed_plate_df['amino-acid_substitutions'] == parent][value_column].values)
                if normalise_method == 'median':
                    parent_mean = np.median(parent_values)
                else:
                    parent_mean = np.mean(parent_values)
                parent_sd = np.std(parent_values)

                vals = list(grp[value_column].values)
                U1, p = None, None
                # Now check if there are 3 otherwise we just do > X S.D over - won't be sig anyway.
                if len(grp) > 2:
                    # Do stats
                    U1, p = mannwhitneyu(parent_values, vals, method="exact")
                mean_vals = np.mean(vals)
                std_vals = np.std(vals)
                median_vals = np.median(vals)
                sig = mean_vals > ((sd_cutoff * parent_sd) + parent_mean)
                rows.append(
                    [value_column, mutation, len(grp), mean_vals, std_vals, median_vals, mean_vals - parent_mean, sig,
                     U1, p])
    stats_df = pd.DataFrame(rows, columns=['value_column', 'amino-acid_substitutions', 'number of wells with amino-acid substitutions', 'mean', 'std',
                                           'median', 'amount greater than parent mean',
                                           f'greater than > {sd_cutoff} parent', 'man whitney U stat', 'p-value'])
    return stats_df


def checkNgen_folder(folder_path: str) -> str:

    """
    Check if the folder and its subfolder exists
    create a new directory if not
    Args:
    - folder_path: str, the folder path
    """
    # get rid of the very first / if it exists
    if folder_path[0] == "/":
        folder_path = folder_path[1:]

    # if input path is file
    if bool(os.path.splitext(folder_path)[1]):
        folder_path = os.path.dirname(folder_path)

    split_list = os.path.normpath(folder_path).split("/")
    for p, _ in enumerate(split_list):
        subfolder_path = "/".join(split_list[: p + 1])
        if not os.path.exists(subfolder_path):
            print(f"Making {subfolder_path} ...")
            os.mkdir(subfolder_path)
    return folder_path


def free_port(port):
    """
    Kill any processes using the specified port.
    If the port is not free after this, return a new available port.
    """
    # Attempt to free the specified port by killing any process using it
    for conn in psutil.net_connections(kind="inet"):
        if conn.laddr.port == port:
            try:
                proc = psutil.Process(conn.pid)
                proc.terminate()  # Attempt to terminate the process
                proc.wait(timeout=3)  # Wait for the process to terminate
                print(f"Terminated process {conn.pid} using port {port}.")
            except (psutil.NoSuchProcess, psutil.AccessDenied, psutil.TimeoutExpired):
                print(f"Could not terminate process {conn.pid} using port {port}.")

    # Check if the specified port is still in use
    if is_port_in_use(port):
        # Find and return the next available port
        port = find_free_port()
        print(f"Port {port} is in use. Using new port {port}.")
    return port


def is_port_in_use(port):
    """Check if a specific port is currently in use."""
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
        return s.connect_ex(("localhost", port)) == 0


def find_free_port():
    """Find an available port."""
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
        s.bind(("", 0))
        s.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
        return s.getsockname()[1]


def work_up_lcms(
    file,
    products,
    substrates=None,
    drop_string=None,
):
    """Works up a standard csv file from Revali.
    Parameters:
    -----------
    file: string
        Path to the csv file
    products: list of strings
        Name of the peaks that correspond to the product
    substrates: list of strings
        Name of the peaks that correspond to the substrate
    drop_string: string, default 'burn_in'
        Name of the wells to drop, e.g., for the wash/burn-in period that are not samples.
    Returns:
    --------
    plate: ns.Plate object (DataFrame-like)
    """
    if isinstance(file, str):
        # Read in the data
        df = pd.read_csv(file, header=[1])
    else:
        # Change to handling both
        df = file
    # Convert nans to 0
    df = df.fillna(0)
    # Only grab the Sample Acq Order No.s that have a numeric value
    index = [True for _ in df["Sample Acq Order No"]]
    for i, value in enumerate(df["Sample Acq Order No"]):
        try:
            int(value)
        except ValueError:
            index[i] = False
    # Index on this
    df = df[index]

    def fill_vial_number(series):
        for i, row in enumerate(series):
            if pd.isna(row):
                series[i] = series[i - 1]
        return series

    df["Sample Vial Number"] = fill_vial_number(df["Sample Vial Number"].copy())
    # Drop empty ones!
    df = df[df["Sample Vial Number"] != 0]
    # Remove unwanted wells
    df = df[df["Sample Name"] != drop_string]
    # Get wells

    df.insert(0, "Well", df["Sample Vial Number"].apply(lambda x: str(x).split("-")[-1]))
    # Rename
    df = df.rename({"Sample Name": "Plate"}, axis="columns")
    # Create minimal DataFrame
    df = df[["Well", "Plate", "Compound Name", "Area"]].reset_index(drop=True)
    # Pivot table; drop redundant values by only taking 'max' with aggfunc
    # (i.e., a row is (value, NaN, NaN) and df is 1728 rows long;
    # taking max to aggregate duplicates gives only (value) and 576 rows long)
    df = df.pivot_table(
        index=["Well", "Plate"], columns="Compound Name", values="Area", aggfunc="max"
    ).reset_index()
    # Get rows and columns
    df.insert(1, "Column", df["Well"].apply(lambda x: int(x[1:]) if x[1:].isdigit() else None))
    df.insert(1, "Row", df["Well"].apply(lambda x: x[0]))
    # Set values as floats
    cols = products + substrates if substrates is not None else products
    for col in cols:
        df[col] = df[col].astype(float)
    plate = ns.Plate(df, value_name=products[-1]).set_as_location("Plate", idx=3)
    plate.values = products
    return plate


def process_files(results_df, plate_df, plate: str, product: list) -> pd.DataFrame:
    """
    Process and combine a single plate file

    Args:
    - product : str
        The name of the product to be analyzed. ie pdt
    - plate : str, ie 'HMC0225_HMC0226.csv'
        The name of the input CSV file containing the plate data.

    Returns:
    - pd.DataFrame
        A pandas DataFrame containing the processed data.
    - str
        The path of the output CSV file containing the processed data.
    """
    filtered_df = results_df[["Plate", "Well", "amino-acid_substitutions", "nt_sequence", "aa_sequence"]]
    filtered_df = filtered_df[(filtered_df["amino-acid_substitutions"] != "#N.A.#")].dropna()

    # Extract the unique entries of Plate
    unique_plates = filtered_df["Plate"].unique()

    # Create an empty list to store the processed plate data
    processed_data = []

    # Iterate over unique Plates and search for corresponding CSV files in the current directory
    plate_object = work_up_lcms(plate_df, product)

    # Extract attributes from plate_object as needed for downstream processes
    if hasattr(plate_object, "df"):
        # Assuming plate_object has a dataframe-like attribute 'df' that we can work with
        plate_df = plate_object.df
        plate_df["Plate"] = plate  # Add the plate identifier for reference

        # Merge filtered_df with plate_df to retain amino-acid_substitutionss and nt_sequence columns
        merged_df = pd.merge(
            plate_df, filtered_df, on=["Plate", "Well"], how="left"
        )
        columns_order = (
                ["Plate", "Well", "Row", "Column", "amino-acid_substitutions"]
                + product
                + ["nt_sequence", "aa_sequence"]
        )
        merged_df = merged_df[columns_order]
        processed_data.append(merged_df)

    # Concatenate all dataframes if available
    if processed_data:
        processed_df = pd.concat(processed_data, ignore_index=True)
    else:
        processed_df = pd.DataFrame(
            columns=["Plate", "Well", "Row", "Column", "amino-acid_substitutions"]
                    + product
                    + ["nt_sequence", "aa_sequence"]
        )

    # Ensure all entries in 'Mutations' are treated as strings
    processed_df["amino-acid_substitutions"] = processed_df["amino-acid_substitutions"].astype(str)

    # Remove any rows with empty values
    processed_df = processed_df.dropna()

    # Return the processed DataFrame for downstream processes
    return processed_df


# Function to process the plate files
def process_plate_files(product: str, input_csv: str) -> pd.DataFrame:

    """
    Process the plate files to extract relevant data for downstream analysis.
    Assume the same directory contains the plate files with the expected names.
    The expected filenames are constructed based on the Plate values in the input CSV file.
    The output DataFrame contains the processed data for the specified product
    and is saved to a CSV file named 'seqfit.csv' in the same dirctory.

    Args:
    - product : str
        The name of the product to be analyzed. ie pdt
    - input_csv : str, ie 'HMC0225_HMC0226.csv'
        The name of the input CSV file containing the plate data.

    Returns:
    - pd.DataFrame
        A pandas DataFrame containing the processed data.
    - str
        The path of the output CSV file containing the processed data.
    """

    dir_path = os.path.dirname(input_csv)
    print(f"Processing data from '{dir_path}'")

    # Load the provided CSV file
    results_df = pd.read_csv(input_csv)

    # Extract the required columns: Plate, Well, amino-acid_substitutionss, and nt_sequence, and remove rows with '#N.A.#' and NaN values
    # barcode_plate	Plate	Well	Alignment Count	nucleotide_amino-acid_substitutions	amino-acid_substitutions	Alignment Probability	Average amino-acid_substitutions frequency	P value	P adj. value	nt_sequence	aa_sequence
    filtered_df = results_df[["Plate", "Well", "amino-acid_substitutions", "nt_sequence", "aa_sequence"]]
    filtered_df = filtered_df[(filtered_df["amino-acid_substitutions"] != "#N.A.#")].dropna()

    # Extract the unique entries of Plate
    unique_plates = filtered_df["Plate"].unique()

    # Create an empty list to store the processed plate data
    processed_data = []

    # Iterate over unique Plates and search for corresponding CSV files in the current directory
    for plate in unique_plates:
        # Construct the expected filename based on the Plate value
        filename = os.path.join(dir_path, f"{plate}.csv")

        # Check if the file exists in the current directory
        if os.path.isfile(filename):
            print(f"Processing data for Plate: {plate}")
            # Work up data to plate object
            plate_object = work_up_lcms(filename, product)

            # Extract attributes from plate_object as needed for downstream processes
            if hasattr(plate_object, "df"):
                # Assuming plate_object has a dataframe-like attribute 'df' that we can work with
                plate_df = plate_object.df
                plate_df["Plate"] = plate  # Add the plate identifier for reference

                # Merge filtered_df with plate_df to retain amino-acid_substitutionss and nt_sequence columns
                merged_df = pd.merge(
                    plate_df, filtered_df, on=["Plate", "Well"], how="left"
                )
                columns_order = (
                    ["Plate", "Well", "Row", "Column", "amino-acid_substitutions"]
                    + product
                    + ["nt_sequence", "aa_sequence"]
                )
                merged_df = merged_df[columns_order]
                processed_data.append(merged_df)

    # Concatenate all dataframes if available
    if processed_data:
        processed_df = pd.concat(processed_data, ignore_index=True)
    else:
        processed_df = pd.DataFrame(
            columns=["Plate", "Well", "Row", "Column", "amino-acid_substitutions"]
            + product
            + ["nt_sequence", "aa_sequence"]
        )

    # Ensure all entries in 'Mutations' are treated as strings
    processed_df["amino-acid_substitutions"] = processed_df["amino-acid_substitutions"].astype(str)

    # Remove any rows with empty values
    processed_df = processed_df.dropna()

    seqfit_path = os.path.join(dir_path, "seqfit.csv")

    # Optionally, save the processed DataFrame to a CSV file
    processed_df.to_csv(seqfit_path, index=False)
    print(f"Processed data saved to {seqfit_path} in the same directory")

    # Return the processed DataFrame for downstream processes
    return processed_df, seqfit_path


def match_plate2parent(df: pd.DataFrame, parent_dict: Optional[Dict] = None) -> dict:

    """
    Find plate names correpsonding to each parent sequence.

    Args:
    - df : pd.DataFrame
        A pandas DataFrame containing the data for a single plate.
        The DataFrame should have the following columns:
        - "Plate" : str
            The plate identifier.
        - "Well" : str
            The well identifier.
        - "Mutations" : str
            The mutations in the well.
    - parent_dict : dict
        A dictionary containing the parent name for each aa_varient.

    Returns:
    - dict
        A dictionary containing the plate names for each parent sequence.
    """

    if parent_dict is None:

        # add aa_sequence column if not present by translating from the nt_sequence column
        if "aa_sequence" not in df.columns:
            df["aa_sequence"] = df["nt_sequence"].apply(
                Bio.sequence.Sequence(df["nt_sequence"]).translate
            )

        # get all the parents from the df
        parents = df[df["amino-acid_substitutions"] == "#PARENT#"].reset_index(drop=True).copy()

        # get the parent nt_sequence
        parent_aas = (
            df[df["amino-acid_substitutions"] == "#PARENT#"][["amino-acid_substitutions", "aa_sequence"]]
            .drop_duplicates()["aa_sequence"]
            .tolist()
        )

        parent_dict = {f"Parent-{i+1}": parent for i, parent in enumerate(parent_aas)}

    # get the plate names for each parent
    parent2plate = {
        p_name: df[df["aa_sequence"] == p_seq]["Plate"].unique().tolist()
        for p_name, p_seq in parent_dict.items()
    }

    # reverse the dictionary to have plate names as keys and rasie flag if there are multiple parents for a plate
    plate2parent = {}
    for parent, plates in parent2plate.items():
        for plate in plates:
            if plate in plate2parent:
                raise ValueError(f"Multiple parents found for plate {plate}")
            else:
                plate2parent[plate] = parent

    return parent_dict, plate2parent


def detect_outliers_iqr(series: pd.Series) -> pd.Index:

    """
    Calculate the Interquartile Range (IQR) and
    determine the lower and upper bounds for outlier detection.

    The IQR is a measure of statistical dispersion and
    is calculated as the difference between the third quartile (Q3)
    and the first quartile (Q1) of the data

    Args:
    - series : pandas.Series
        A pandas Series containing the data for which the IQR and bounds are to be calculated.

    Returns:
    - tuple
        A tuple containing the lower bound and upper bound for outlier detection.

    Example:
    --------
    >>> import pandas as pd
    >>> data = pd.Series([10, 12, 14, 15, 18, 20, 22, 23, 24, 25, 100])
    >>> calculate_iqr_bounds(data)
    (-1.0, 39.0)
    """

    Q1 = series.quantile(0.25)
    Q3 = series.quantile(0.75)
    IQR = Q3 - Q1
    lower_bound = Q1 - 1.5 * IQR
    upper_bound = Q3 + 1.5 * IQR

    return series[(series < lower_bound) | (series > upper_bound)].index


def norm2parent(plate_df: pd.DataFrame) -> pd.DataFrame:

    """
    For each given plate,
    normalize the pdt values of a plate to the mean of the parent
    without the outliers.

    Args:
    - plate_df : pd.DataFrame
        A pandas DataFrame containing the data for a single plate.
        The DataFrame should have the following columns:
        - "Plate" : str
            The plate identifier.
        - "Mutations" : str
            The mutations in the well.
        - "pdt" : float
            The pdt value for the well.

    Returns:
    - pd.DataFrame
        A pandas DataFrame containing the normalized pdt values.
    """

    # get all the parents from the df
    parents = (
        plate_df[plate_df["amino-acid_substitutions"] == "#PARENT#"].reset_index(drop=True).copy()
    )
    filtered_parents = (
        parents.drop(index=detect_outliers_iqr(parents["pdt"]))
        .reset_index(drop=True)
        .copy()
    )

    # normalize the whole plate to the mean of the filtered parent
    plate_df["pdt_norm"] = plate_df["pdt"] / filtered_parents["pdt"].mean()

    return plate_df


def process_mutation(mutation: str) -> pd.Series:
    # Check if mutation is #PARENT#
    if mutation == "#PARENT#":
        return pd.Series([0, [(None, None, None)]])  # Return 0 sites and NaN details

    # Split by "_" to get number of sites
    sites = mutation.split("_")
    num_sites = len(sites)

    # Extract details if it matches the pattern
    details = []
    for site in sites:
        match = re.match(r"^([A-Z])(\d+)([A-Z*])$", site)
        if match:
            parent_aa, site_number, mutated_aa = match.groups()
            details.append((parent_aa, site_number, mutated_aa))
        else:
            details.append((None, None, None))

    return pd.Series([num_sites, details])


def prep_single_ssm(df: pd.DataFrame) -> pd.DataFrame:
    """
    Prepare the data for a single sitessm summary plot.

    Args:
    - df: pd.DataFrame, input full dataframe

    Returns:
    - pd.DataFrame, output dataframe
    """

    # slice out single site SSM and add in parentAA, site, and mutAA columns
    single_ssm_df = df[df["num_sites"] <= 1].copy()

    # Expand the single entry in Details for these rows into three columns
    single_ssm_df[["parent_aa", "site_numb", "mut_aa"]] = pd.DataFrame(
        single_ssm_df["mut_dets"].apply(lambda x: x[0]).tolist(),
        index=single_ssm_df.index,
    )

    single_ssm_df["parent_aa_loc"] = (
        single_ssm_df["parent_aa"] + single_ssm_df["site_numb"]
    )

    # fill nan site numbers with 0 and convert to int
    single_ssm_df["site_numb"] = single_ssm_df["site_numb"].fillna(0).astype(int)

    return single_ssm_df


def get_single_ssm_site_df(
    single_ssm_df: pd.DataFrame, parent: str, site: str
) -> pd.DataFrame:
    """
    Get the single site SSM data for a given site with appended parent data.

    Args:
    - single_ssm_df: pd.DataFrame, input single site SSM dataframe
    - parent: str, parent to filter the data on
    - site: str, site to filter the data on

    Returns:
    - pd.DataFrame, output dataframe
    """

    # get the site data
    site_df = (
        single_ssm_df[
            (single_ssm_df["Parent_Name"] == parent)
            & (single_ssm_df["parent_aa_loc"] == site)
        ]
        .reset_index(drop=True)
        .copy()
    )

    # get parents from those plates
    site_parent_df = (
        single_ssm_df[
            (single_ssm_df["amino-acid_substitutions"] == "#PARENT#")
            & (single_ssm_df["Plate"].isin(site_df["Plate"].unique()))
        ]
        .reset_index(drop=True)
        .copy()
    )

    # rename those site_numb, mut_aa, parent_aa_loc None or NaN to corresponding parent values
    site_parent_df["mut_aa"] = site_parent_df["mut_aa"].fillna(
        site_df["parent_aa"].values[0]
    )
    site_parent_df["site_numb"] = site_parent_df["site_numb"].fillna(
        site_df["site_numb"].values[0]
    )
    site_parent_df["parent_aa_loc"] = site_parent_df["parent_aa_loc"].fillna(
        site_df["parent_aa_loc"].values[0]
    )

    # now merge the two dataframes
    return pd.concat([site_parent_df, site_df]).reset_index(drop=True).copy()


def prep_aa_order(df: pd.DataFrame, add_na: bool = False) -> pd.DataFrame:
    """
    Prepare the data for a single sitessm summary plot.

    Args:
    - df: pd.DataFrame, input full dataframe

    Returns:
    - pd.DataFrame, output dataframe
    """

    # Define the order of x-axis categories
    x_order = list(AA_DICT.values())

    if add_na:
        x_order += ["#N.A.#"]

    # Convert `Mutations` to a categorical column with specified order
    df["mut_aa"] = pd.Categorical(df["mut_aa"], categories=x_order, ordered=True)

    # Sort by the `x_order`, filling missing values
    return (
        df.sort_values("mut_aa", key=lambda x: x.cat.codes)
        .reset_index(drop=True)
        .copy()
    )


def get_parent2sitedict(df: pd.DataFrame) -> dict:

    """
    Get a dictionary of parent to site mapping for single site mutants.

    Args:
    - df : pd.DataFrame

    Returns:
    - dict
        A dictionary containing the parent sequence and site number for each parent.
    """

    site_dict = deepcopy(
        df[["Parent_Name", "parent_aa_loc"]]
        .drop_duplicates()
        .dropna()
        .groupby("Parent_Name")["parent_aa_loc"]
        .apply(list)
        .to_dict()
    )

    # Sort the site list for each parent as an integer
    for parent, sites in site_dict.items():
        # Ensure each site is processed as a string and sorted by the integer part
        site_dict[parent] = sorted(sites, key=lambda site: int(str(site)[1:]))

    return site_dict


def get_x_label(x: str):
    
    """
    Function to return the x-axis label based on the input string.
    """

    if "mut_aa" in x.lower():
        clean_x = x.replace("mut_aa", "Amino acid substitutions")
    else:
        clean_x = x.replace("_", " ").capitalize()

    return clean_x


def get_y_label(y: str):

    """
    Function to return the y-axis label based on the input string.
    """
    clean_y = ""
    if "pdt" in y.lower():
        clean_y = "Product"
    elif "area" in y.lower():
        clean_y = "Yield"
    elif y == "fitness_ee2/(ee1+ee2)":
        clean_y = "ee2/(ee1+ee2)"
    elif y == "fitness_ee1/(ee1+ee2)":
        clean_y = "ee1/(ee1+ee2)"
    else:
        clean_y = y

    # normalize the y label
    if "norm" in y.lower():
        clean_y = f"Normalized {clean_y.lower()}"
    return clean_y


def plot_bar_point(
    df: pd.DataFrame,
    x: str,
    y: str,
    x_label: str = None,
    y_label: str = None,
    title: str = None,
    if_max: bool = False,
) -> hv.Layout:

    # Create Bars plot
    bars = hv.Bars(
        df[[y, x]].sort_values(x).groupby(x).mean(),
        kdims=x,
        vdims=y,
    )

    # Display the plot
    bars.opts(
        title=title,
        xlabel=x_label or get_x_label(x),
        ylabel=y_label or get_y_label(y),
        color=y,
        cmap="coolwarm",
        width=600,
        height=400,
        xrotation=45,
    )

    # Create Scatter chart
    points = hv.Scatter(df, x, [y, "Plate", "Well"]).opts(
        color=y, cmap="gray", size=8, alpha=0.5, tools=["hover"]
    )

    # create another scatter plot to highlight the max value
    if if_max:
        max_points = hv.Scatter(
            df.loc[df.groupby(x)[y].idxmax()],
            x,
            [y, "Plate", "Well"],
        ).opts(color="orange", size=10, alpha=1, tools=["hover"])
        return bars * points * max_points

    else:
        return bars * points


def get_parent_plot(df: pd.DataFrame, y: str = "pdt_norm") -> hv.Bars:

    """
    Function to plot the max value by parent.

    Args:
    - df : pd.DataFrame
        A pandas DataFrame containing the data for all parents.
        The DataFrame should have the Parent_Name columns
    - y : str
        The column name for which the max value is to be calculated.

    Returns:
    - hv.Bars
        A holoviews Bars object containing the plot.
    """

    parent_summary = df.groupby("Parent_Name")[y].max().reset_index()
    return hv.Bars(parent_summary, kdims="Parent_Name", vdims=y).opts(
        title="Max Value by Parent", width=600, height=400
    )


def agg_parent_plot(df: pd.DataFrame, ys: list = ["pdt_norm"]) -> pn.Row:

    """
    Function to plot the max value by parent for different y metrics.

    Args:
    - df : pd.DataFrame
        A pandas DataFrame containing the data for all parents.
        The DataFrame should have the Parent_Name columns
    - ys : list
        The list of column name for which the max value is to be calculated.

    Returns:
    - hv.Bars
    """

    # find single site mutations
    # avg_parnet_plots = [get_parent_plot(y=y) for y in ys if y in df.columns]
    avg_parnet_plots = [
        plot_bar_point(
            df,
            x="Parent_Name",
            y=y,
            title=f"{get_y_label(y)} across parents",
            if_max=True,
        )
        for y in ys
        if y in df.columns
    ]

    if len(avg_parnet_plots) == 0:
        return None
    # elif len(avg_ssm_plots) == 1:
    #     return avg_ssm_plots[0]
    else:
        return pn.Row(*avg_parnet_plots)


def plot_single_ssm_avg(
    single_ssm_df: pd.DataFrame,
    parent_name: str,
    y: str = "pdt_norm",
    width: int = 600,
):
    """
    Function to plot single site mutations with average values.

    Parameters:
    - df: DataFrame containing mutation data.
    """

    sliced_df = prep_aa_order(
        single_ssm_df[single_ssm_df["Parent_Name"] == parent_name].copy()
    )

    height = max(20 * sliced_df["site_numb"].nunique() + 60, 160)

    return hv.HeatMap(
        data=sliced_df[["parent_aa_loc", "mut_aa", y]]
        .dropna()
        .groupby(by=["parent_aa_loc", "mut_aa"])
        .mean()
        .sort_values(
            ["parent_aa_loc", "mut_aa"],
            key=lambda col: col.str.extract(r"(\d+)$").fillna(0).astype(int).iloc[:, 0]
            if col.name == "parent_aa_loc"
            else col
        )
        .reset_index(),
        kdims=["mut_aa", "parent_aa_loc"],
        vdims=[y],
    ).opts(
        height=height,
        width=width,
        cmap="coolwarm",
        colorbar=True,
        colorbar_opts=dict(title=get_y_label(y), width=8),
        xrotation=45,
        title=f"Average single site substitution for {parent_name}",
        xlabel="Amino acid substitutions",
        ylabel="Position",
        invert_yaxis=True,
        tools=["hover"],
    )


def agg_single_ssm_exp_avg(
    single_ssm_df: pd.DataFrame,
    parent_name: str,
    ys: list = ["pdt_norm"],
):

    # find single site mutations
    avg_ssm_plots = [
        plot_single_ssm_avg(single_ssm_df=single_ssm_df, parent_name=parent_name, y=y)
        for y in ys
        if y in single_ssm_df.columns
    ]

    if len(avg_ssm_plots) == 0:
        return None
    # elif len(avg_ssm_plots) == 1:
    #     return avg_ssm_plots[0]
    else:
        return pn.Row(*avg_ssm_plots)


def canonicalize_smiles(smiles_string: str) -> str:

    """
    A function to canonicalize a SMILES string.

    Args:
    - smiles_string (str): The input SMILES string.

    Returns:
    - str: The canonicalized SMILES string.
    """

    molecule = Chem.MolFromSmiles(smiles_string)
    if molecule:
        canonical_smiles = Chem.MolToSmiles(molecule, canonical=True)
        return canonical_smiles


def get_chaistructure(
    chai_dir: str,
    seq: str,
    seq_name: str,
    smiles: str = "",
    smiles_name: str = "",
    cofactor_smiles: Union[List, str] = "",
    cofactor_name: str = "",
    joinsubcofactor: bool = True,
    torch_device: str = "cuda",
    ifrerun: bool = False,
):

    """
    A function for getting chai structure for a gvien sequence.

    Args:
    - seq (str): sequence of the protein
    - seq_name (str): label for the protein-ligand pair
    - smiles (str): SMILES string of the ligand
    - smiles_name (str): label for the ligand
    - cofactor_smiles (list or str): list of SMILES strings of cofactors, default is ""
    - cofactor_name (str): label for the cofactor, default is ""
    - joinsubcofactor (bool): whether to join the substrate and cofactor in the same fasta file, default is True

    Returns:
    - list: list of the output files
    """

    chai_dir = checkNgen_folder(os.path.normpath(chai_dir))

    # make sure output dir is dir
    output_subdir = os.path.join(chai_dir, seq_name)

    # Need to clean up the sequence
    seq = seq.strip().replace("*", "").replace(" ", "").upper()

    input_fasta = f">protein|{seq_name}\n{seq}\n"

    if cofactor_smiles != "":
        # convert cofactor_smiles to a list if it is a string
        if isinstance(cofactor_smiles, str):
            # use ast.literal_eval to convert string to list
            try:
                cofactor_smiles = literal_eval(cofactor_smiles)
            except:
                cofactor_smiles = [cofactor_smiles]

        # add cofactor SMILES to the fasta
        for cofactor_smile in cofactor_smiles:
            input_fasta += f">ligand|{seq_name}-cofactor\n{cofactor_smile}\n"

    if smiles:
        smiles = canonicalize_smiles(smiles)
        # now add substrate
        input_fasta += f">ligand|{seq_name}-substrate\n{smiles}\n"

    # only rerun if the flag is set and the output folder doies not exists
    if ifrerun or not os.path.exists(output_subdir):

        output_subdir = Path(checkNgen_folder(output_subdir))

        fasta_path = Path(f"{output_subdir}/{seq_name}.fasta")
        fasta_path.write_text(input_fasta)

        output_paths = run_inference(
            fasta_file=fasta_path,
            output_dir=output_subdir,
            # 'default' setup
            num_trunk_recycles=3,
            num_diffn_timesteps=200,
            seed=42,
            device=torch.device(torch_device),
            use_esm_embeddings=True,
        )

        renamed_output_files = []

        # get name of the output cif or pdb files
        output_strcut_files = sorted(
            glob(f"{output_subdir}/*.cif") + glob(f"{output_subdir}/*.pdb")
        )

        # rename the output files cif or pdb files
        for output_strcut_file in output_strcut_files:
            renamed_output_file = output_strcut_file.replace("pred.model_idx", seq_name)
            os.rename(
                output_strcut_file,
                renamed_output_file,  # output_strcut_file.replace("pred.model_idx", seq_name)
            )
            renamed_output_files.append(renamed_output_file)

        renamed_scores_files = []

        # for npz files do the same
        output_scores_files = sorted(glob(f"{output_subdir}/*.npz"))

        for output_scores_file in output_scores_files:
            renamed_output_file = output_scores_file.replace(
                "scores.model_idx", seq_name
            )
            os.rename(
                output_scores_file,
                renamed_output_file,  # output_scores_file.replace("scores.model_idx", seq_name)
            )
            renamed_scores_files.append(renamed_output_file)

    else:
        renamed_output_files = glob(f"{output_subdir}/*.cif") + glob(
            f"{output_subdir}/*.pdb"
        )
        renamed_scores_files = glob(f"{output_subdir}/*.npz")

    return renamed_output_files, renamed_scores_files


def export_structure_as_html(
    parent_name, file_path, output_dir="", highlight_residues=None
):
    """
    Exports the 3D structure as an interactive HTML file with highlighted residues.

    Parameters:
    - parent_name: str, the title/name for the viewer.
    - file_path: str, path to the PDB or CIF structure file.
    - output_dir: str, directory to save the HTML file (creates the directory if it doesn't exist).
    - highlight_residues: list of dicts, each with 'chain' and 'resi' keys to specify residues to highlight.
    """

    output_filename = f"{parent_name}_structure.html"

    # Read the structure file content
    with open(file_path, "r") as file:
        file_content = file.read()

    # Determine the output path
    if output_dir:
        os.makedirs(
            output_dir, exist_ok=True
        )  # Create the directory if it doesn't exist
        output_path = os.path.join(output_dir, output_filename)
    else:
        output_path = output_filename  # Save in the current directory by default

    # Convert highlight_residues list to JavaScript array format
    residues_js = (
        ", ".join(
            f"{{chain: '{res['chain']}', resi: {res['resi']}}}"
            for res in highlight_residues
        )
        if highlight_residues
        else ""
    )

    # Define the HTML content
    html_content = f"""
    <!DOCTYPE html>
    <html lang="en">
    <head>
        <meta charset="UTF-8">
        <meta name="viewport" content="width=device-width, initial-scale=1.0">
        <title>3Dmol Structure Viewer</title>
        <script src="https://3Dmol.org/build/3Dmol-min.js"></script>
    </head>
    <body>
        <div id="viewer" style="width: 600px; height: 360px; position: relative;"></div>
        <script>
            function loadStructure(content, format) {{
                var viewer = $3Dmol.createViewer("viewer", {{ backgroundColor: "white" }});
                viewer.addModel(content, format);
                viewer.setStyle({{}}, {{ stick: {{}} }});
                viewer.setStyle({{ chain: 'A' }}, {{ cartoon: {{}} }});
                
                // Highlight specified residues
                const residuesToHighlight = [{residues_js}];
                residuesToHighlight.forEach(res => {{
                    viewer.addStyle(
                        {{chain: res.chain, resi: res.resi, atom: "CA"}}, 
                        {{sphere: {{radius: 1.5, color: "orange"}}}}
                    );
                }});

                viewer.zoomTo();
                viewer.render();
            }}

            const structureContent = `{file_content}`;
            const modelType = "{'pdb' if file_path.endswith('.pdb') else 'cif'}";
            loadStructure(structureContent, modelType);
        </script>
    </body>
    </html>
    """

    # Write the HTML content to the file
    with open(output_path, "w") as output_file:
        output_file.write(html_content)

    print(f"Exported structure view to {output_path}")
    return output_path


def gen_seqfitvis(
    seqfit_path: str,
    output_dir: str = "",
    protein_chain: str = "A",
    chai_meta_data={
        "smiles": "",
        "smiles_name": "",
        "cofactor_smiles": "",
        "cofactor_name": "",
        "joinsubcofactor": True,
        "torch_device": "cuda",
        "ifrerun": False,
    },
    gen_struct: bool = False,
    # port=8000,
):

    # normalized per plate to parent
    if isinstance(seqfit_path, str):
        df = pd.read_csv(seqfit_path)
    else:
        df = seqfit_path
    # ignore deletion meaning "Mutations" == "-"
    df = df[df["amino-acid_substitutions"] != "-"].copy()
    # count number of sites mutated and append mutation details
    # df["num_sites"] = df['Mutations'].apply(lambda x: 0 if x == "#PARENT#" else len(x.split("_")))

    # Apply function to the column
    df[["num_sites", "mut_dets"]] = df["amino-acid_substitutions"].apply(process_mutation)

    # apply the norm function to all plates
    df = df.groupby("Plate").apply(norm2parent).reset_index(drop=True).copy()

    # add a new column called parent name to the df
    # using the dict out put from match_plate2parent
    # that matches the plate to the parent
    parent_dict, plate2parent = match_plate2parent(df, parent_dict=None)
    df["Parent_Name"] = df["Plate"].map(plate2parent)

    parents = df["Parent_Name"].unique().tolist()
    single_ssm_df = prep_single_ssm(df)
    sites_dict = get_parent2sitedict(single_ssm_df)

    struct_dict = {}

    if gen_struct:
        # get structures for all parents
        for parent_name, parent_seq in parent_dict.items():
            # get parent chai structure
            chai_files, chai_scores_files = get_chaistructure(
                chai_dir=os.path.join(os.path.dirname(seqfit_path), "chai"),
                seq=parent_seq,
                seq_name=parent_name,
                **chai_meta_data,
            )

            # Export the 3D structure as an interactive HTML file
            html_path = export_structure_as_html(
                parent_name=parent_name,
                file_path=chai_files[0],
                output_dir=os.path.dirname(seqfit_path),
                highlight_residues=[
                    {"chain": protein_chain, "resi": int(i[1:])}
                    for i in sites_dict.get(parent_name, [])
                ],
            )

            struct_dict[parent_name] = html_path

    def get_subplots(
        parent,
    ):

        if parent not in struct_dict:
            html_pane = None
        else:
            # Load HTML content from a file
            with open(struct_dict[parent], "r") as file:
                html_content = file.read()

            # URL-encode the HTML content
            data_url = "data:text/html;charset=utf-8," + urllib.parse.quote(
                html_content
            )

            # Embed in an iframe
            iframe_code = f'<iframe id="viewerFrame" width="600" height="400" src="{data_url}"></iframe>'

            # Create an HTML pane with the loaded content
            html_pane = pn.pane.HTML(iframe_code)

        # Get the list of sites for the selected parent
        site_options = sites_dict.get(parent, [])

        # Set the initial site to the first item in the list if it exists
        initial_site = site_options[0] if site_options else None

        # Create a site dropdown with the initial site as the default
        site_dropdown = pn.widgets.Select(
            name="Sites", options=site_options, value=initial_site
        )

        def update_site_plot(site):

            site_df = prep_aa_order(
                get_single_ssm_site_df(single_ssm_df, parent=parent, site=site)
            )

            if site_df.empty:
                return pn.pane.Markdown("### No data available for the selected site")

            site_info = (
                site_df["parent_aa_loc"].unique()[0] if not site_df.empty else "Unknown"
            )

            return plot_bar_point(
                df=site_df,
                x="mut_aa",
                y="pdt_norm",
                # y_label: str = None,
                title=f"{site_info} for {parent}",
                if_max=False,
            )

        site_plot = pn.Column(pn.bind(update_site_plot, site=site_dropdown))

        return pn.Column(
            html_pane,
            agg_single_ssm_exp_avg(
                single_ssm_df=single_ssm_df,
                parent_name=parent,
                # ys: list,
            ),
            site_dropdown,
            site_plot,
        )

    # Dropdown for parent selection
    parent_dropdown = pn.widgets.Select(name="Parent", options=parents)

    # Initial parent plots
    initial_subplots = get_subplots(parents[0])

    # Panel layout
    dashboard = pn.Column(
        agg_parent_plot(df),
        parent_dropdown,
        pn.Column(pn.bind(get_subplots, parent=parent_dropdown)),
    )

    port = find_free_port()

    if not output_dir:
        output_dir = os.path.dirname(seqfit_path)
    else:
        output_dir = checkNgen_folder(os.path.normpath(output_dir))

    # Serve the dashboard
    pn.serve(dashboard, port=port, open=True, show=True, start=True)

    print(f"Access the interactive dashboard via: http://localhost:{port}")