"""
A script for visualzing the sequnece-fitness relationship
"""

import re
import os
import html

from pathlib import Path
from copy import deepcopy
from glob import glob

from ast import literal_eval
import urllib.parse

import numpy as np
import pandas as pd
import biopandas as Bio

from rdkit import Chem

import panel as pn
import holoviews as hv

from bokeh.io import output_notebook, show
from bokeh.plotting import figure

import torch
import torch.nn as nn

from chai_lab.chai1 import run_inference

# Enable Bokeh to display plots in the notebook
hv.extension('bokeh')
pn.extension()
output_notebook()

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


def match_plate2parent(df: pd.DataFrame, parent_dict: dict | None = None) -> dict:

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

        # add aa_variant column if not present by translating from the nc_variant column
        if "aa_variant" not in df.columns:
            df["aa_variant"] = df["nc_variant"].apply(
                Bio.sequence.Sequence(df["nc_variant"]).translate
            )

        # get all the parents from the df
        parents = df[df["Mutations"] == "#PARENT#"].reset_index(drop=True).copy()

        # get the parent nc_variant
        parent_aas = (
            df[df["Mutations"] == "#PARENT#"][["Mutations", "aa_variant"]]
            .drop_duplicates()["aa_variant"]
            .tolist()
        )

        parent_dict = {f"Parent-{i+1}": parent for i, parent in enumerate(parent_aas)}

    # get the plate names for each parent
    parent2plate = {
        p_name: df[df["aa_variant"] == p_seq]["Plate"].unique().tolist()
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
        plate_df[plate_df["Mutations"] == "#PARENT#"].reset_index(drop=True).copy()
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

def get_single_ssm_site_df(single_ssm_df: pd.DataFrame, parent: str, site: str) -> pd.DataFrame:
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
    site_df = single_ssm_df[
        (single_ssm_df["Parent_Name"] == parent)
        & (single_ssm_df["parent_aa_loc"] == site)
    ].copy()

    # get parents from those plates
    site_parent_df = single_ssm_df[
        (single_ssm_df["Mutations"] == "#PARENT#")
        & (single_ssm_df["Plate"].isin(site_df["Plate"].unique()))
    ].copy()

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
    df["mut_aa"] = pd.Categorical(
        df["mut_aa"], categories=x_order, ordered=True
    )

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
        .drop_duplicates().dropna()
        .groupby("Parent_Name")["parent_aa_loc"]
        .apply(list)
        .to_dict()
    )

    # Sort the site list for each parent as an integer
    for parent, sites in site_dict.items():
        # Ensure each site is processed as a string and sorted by the integer part
        site_dict[parent] = sorted(sites, key=lambda site: int(str(site)[1:]))

    return site_dict


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
        clean_y=  y

    # normalize the y label
    if "norm" in y.lower():
        clean_y = f"Normalized {clean_y.lower()}"
    return clean_y


def plot_bar_point(
    df: pd.DataFrame,
    x: str,
    y: str,
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

    sliced_df = prep_aa_order(single_ssm_df[single_ssm_df["Parent_Name"] == parent_name].copy())

    height = max(30 * sliced_df["site_numb"].nunique(), 160)

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
            # key=lambda col: col[1:].astype(int)
            # if col.name == "single_mutated_sites_w_parent"
            # else col,
        )
        .reset_index(),
        kdims=["mut_aa", "parent_aa_loc"],
        vdims=[y],
    ).opts(
        height=height,
        width=width,
        cmap="coolwarm",
        # color_levels=color_levels,
        colorbar=True,
        colorbar_opts=dict(title=y, width=8),
        xrotation=45,
        title=f"Average single site mutations for {parent_name}",
        xlabel="Residue",
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
    cofactor_smiles: list | str = "",
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
        output_strcut_files = sorted(glob(f"{output_subdir}/*.cif") + glob(
            f"{output_subdir}/*.pdb"
        ))

        # rename the output files cif or pdb files
        for output_strcut_file in output_strcut_files:
            renamed_output_file = output_strcut_file.replace(
                "pred.model_idx", seq_name
            )
            os.rename(
                output_strcut_file, renamed_output_file # output_strcut_file.replace("pred.model_idx", seq_name)
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
                output_scores_file, renamed_output_file # output_scores_file.replace("scores.model_idx", seq_name)
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



class SeqFitVis:
    def __init__(self, seqfit_path: str, protein_chain: str = "A", output_dir: str = "", chai_meta_data: dict = {
            "smiles": "",
            "smiles_name": "",
            "cofactor_smiles": "",
            "cofactor_name": "",
            "joinsubcofactor": True,
            "torch_device": "cuda",
            "ifrerun": False,
        }):

        """
        Initialize the SeqFitVis class.
        
        Args:
        - seqfit_path (str): Path to the sequence fitness data file.
        - protein_chain (str): Chain identifier for the protein structure.
        
        """

        self._seqfit_path = seqfit_path
        self._protein_chain = protein_chain

        if not output_dir:
            self._output_dir = os.path.dirname(seqfit_path)
        else:
            self._output_dir = checkNgen_folder(output_dir)

        self._chai_meta_data = chai_meta_data

        self._df, self._parent_dict = self._preprocess_data()

        self._parents = self._df["Parent_Name"].unique().tolist()
        self._single_ssm_df = prep_single_ssm(self._df)
        self._sites_dict = get_parent2sitedict(self._single_ssm_df)

        print("Generating parents structures...")
        self._struct_dict = self._gen_structure()

        self._save_serve_dashboard()

    def _preprocess_data(self):

        """
        Preprocess the sequence fitness data.
        """

        df = pd.read_csv(self._seqfit_path)
        # ignore deletion meaning "Mutations" == "-"
        df = df[df["Mutations"] != "-"].copy()
        # count number of sites mutated and append mutation details

        # Apply function to the column
        df[["num_sites", "mut_dets"]] = df["Mutations"].apply(process_mutation)

        # apply the norm function to all plates
        df = df.groupby("Plate").apply(norm2parent).reset_index(drop=True).copy()

        # add a new column called parent name to the df
        # using the dict out put from match_plate2parent
        # that matches the plate to the parent
        parent_dict, plate2parent = match_plate2parent(df, parent_dict=None)
        df["Parent_Name"] = df["Plate"].map(plate2parent)

        return df.copy(), parent_dict


    def _gen_structure(self):

        struct_dict = {}

        # get structures for all parents
        for parent_name, parent_seq in self._parent_dict.items():
        # get parent chai structure
            chai_files, chai_scores_files = get_chaistructure(
                chai_dir=os.path.join(os.path.dirname(self._seqfit_path), "chai"),
                seq=parent_seq,
                seq_name=parent_name,
                **self._chai_meta_data,
            )

            # Export the 3D structure as an interactive HTML file
            html_path = export_structure_as_html(
                parent_name=parent_name,
                file_path=chai_files[0],
                output_dir=os.path.dirname(self._seqfit_path),
                highlight_residues=[
                    {"chain": self._protein_chain, "resi": int(i[1:])} for i in self._sites_dict.get(parent_name, [])
                ],
            )

            struct_dict[parent_name] = html_path

        return struct_dict


    def _get_subplots(
        self,
        parent,
    ):

        # Load HTML content from a file
        with open(self._struct_dict[parent], "r") as file:
            html_content = file.read()
        
        # URL-encode the HTML content
        data_url = "data:text/html;charset=utf-8," + urllib.parse.quote(html_content)

        # Embed in an iframe
        iframe_code = f'<iframe id="viewerFrame" width="600" height="400" src="{data_url}"></iframe>'

        # Create an HTML pane with the loaded content
        html_pane = pn.pane.HTML(iframe_code)

        site_dropdown = pn.widgets.Select(name="Sites", options=self._sites_dict.get(parent, []))

        def update_site_plot(site):

            site_df = prep_aa_order(
                get_single_ssm_site_df(self._single_ssm_df, parent=parent, site=site)
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
                single_ssm_df=self._single_ssm_df,
                parent_name=parent,
                # ys: list,
            ),
            site_dropdown,
            site_plot,
        )

    def _save_serve_dashboard(self):

        # Dropdown for parent selection
        parent_dropdown = pn.widgets.Select(name="Parent", options=self._parents)

        # Initial parent plots
        initial_subplots = self._get_subplots(self._parents[0])

        # Panel layout
        dashboard = pn.Column(
            agg_parent_plot(self._df),
            parent_dropdown,
            pn.Column(pn.bind(self._get_subplots, parent=parent_dropdown)),
        )

        # Serve the dashboard
        port = 8000  # Specify the port you want to use
        server = pn.serve(dashboard, port=port, show=False, start=False)

        # Generate HTML wrapper file with iframe pointing to the served Panel app
        wrapper_html_content = f"""
        <!DOCTYPE html>
        <html lang="en">
        <head>
            <meta charset="UTF-8">
            <meta name="viewport" content="width=device-width, initial-scale=1.0">
            <title>Interactive Panel Dashboard</title>
        </head>
        <body>
            <iframe src="http://localhost:{port}" style="width: 100%; height: 100vh; border: none;"></iframe>
        </body>
        </html>
        """

        # Define the path where you want to save the HTML file
        output_path = os.path.join(checkNgen_folder(self._output_dir), "seqfit.html")
        with open(output_path, "w") as f:
            f.write(wrapper_html_content)

        print(f"Wrapper HTML saved at {output_path}")
        print(f"Access the interactive dashboard via: http://localhost:{port}")

        pn.serve(dashboard, open=True, browser="chrome")

