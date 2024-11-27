###############################################################################
#                                                                             #
#    This program is free software: you can redistribute it and/or modify     #
#    it under the terms of the GNU General Public License as published by     #
#    the Free Software Foundation, either version 3 of the License, or        #
#    (at your option) any later version.                                      #
#                                                                             #
#    This program is distributed in the hope that it will be useful,          #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of           #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            #
#    GNU General Public License for more details.                             #
#                                                                             #
#    You should have received a copy of the GNU General Public License        #
#    along with this program. If not, see <http://www.gnu.org/licenses/>.     #
#                                                                             #
###############################################################################

# Import all packages
from __future__ import annotations

from collections import Counter

import warnings

import os
import sys
from copy import deepcopy

import pandas as pd
import numpy as np

from Bio import AlignIO
from Bio.motifs import Motif
from Bio.PDB.Polypeptide import aa1
from Bio.Align import MultipleSeqAlignment

import matplotlib.pyplot as plt
import matplotlib as mpl

import holoviews as hv
from holoviews.streams import Tap

import ninetysix as ns
import colorcet as cc

from bokeh.plotting import figure
from bokeh.models import (
    ColumnDataSource,
    Range1d,
    CustomJS,
    RangeSlider,
    TapTool,
    HoverTool,
    Label,
    Div,
    FactorRange,
    Legend,
    LegendItem,
)
from bokeh.models.glyphs import Text, Rect
from bokeh.layouts import column, gridplot, row, Spacer
from bokeh.events import Tap
from bokeh.io import save, show, output_file, output_notebook

import panel as pn

from levseq.utils import *

output_notebook()

pn.extension()
pn.config.comms = "vscode"

hv.extension("bokeh")
hv.renderer("bokeh").webgl = True

# warnings.filterwarnings("ignore")
#warnings.filterwarnings("ignore", category=Warning)

#with warnings.catch_warnings():
#    warnings.simplefilter("ignore")

# Redirect stderr to devnull to suppress any remaining output
#sys.stderr = open(os.devnull, "w")

######## Define constants for MSA alignments ########

# Set light gray for seq matched with the reference
match_color = "#d9d9d9"

# Set nuecleotide colors for the MSA alignment plot
NUC_COLOR_DICT = {
    "A": "green",
    "T": "red",
    "G": "black",
    "C": "blue",
    "-": "white",
    "N": "gray",
}


def get_well_ids():
    """
    Generate a list of well IDs for a 96-well plate.
    """
    # Initialize an empty list to store
    well_ids = []

    # Loop through the rows (A to H)
    for row in "ABCDEFGH":
        # Loop through the columns (1 to 12)
        for col in range(1, 13):
            # Combine the row and column to form the well position
            well_ids.append(f"{row}{col}")
    return deepcopy(well_ids)


WELL_IDS = get_well_ids()


def well2nb(well_id):
    """
    Given a well ID, return the row and column.
    """
    row = ord(well_id[0]) - 64
    col = int(well_id[1:])
    nb = (row - 1) * 12 + col
    return f"NB{'0' if nb < 10 else ''}{nb}"


# Function for making plate maps
def _make_platemap(df, title, cmap=None):
    """Generates a plate heatmap from LevSeq data using Holoviews with
    bokeh backend.

    Called via `generate_platemaps`; see docs there.
    """
    # Handle empty dataframe case
    if df.empty:
        # Create a dummy plot with a message
        empty_df = pd.DataFrame({
            'Row': list('ABCDEFGH'),
            'Column': [str(i) for i in range(1, 13)],
            'logseqdepth': [0] * 96,
            'amino-acid_substitutions': [''] * 96,
            'Alignment Count': [0] * 96,
            'Alignment Probability': [0] * 96
        })
        df = empty_df

    # Convert SeqDepth to log for easier visualization.
    df["logseqdepth"] = np.log(
        df["Alignment Count"],
        out=np.zeros_like(df["Alignment Count"], dtype=float),
        where=(df["Alignment Count"] != 0),
    )

    # Create necessary Row and Column values and sort
    df = df.sort_values(["Column", "Row"])
    df["Column"] = df["Column"].astype("str")

    # Set some base opts
    opts = dict(invert_yaxis=True, title=title, show_legend=True)

    # Set the center and handle empty or zero-only data
    if df["logseqdepth"].max() <= 0:
        # For empty data, create minimal valid range
        df["logseqdepth"] = 1  # Set to constant value
        color_levels = [0.9, 1, 1.1]  # Create minimal valid range
        seq_depth_cmap = ['#f7f7f7', '#f7f7f7']  # Use same color for uniform appearance
    else:
        # Regular case with actual data
        seq_depth_cmap = list(reversed(cc.CET_D9))
        center = np.log(10)
        add_min = False
        
        if df["logseqdepth"].min() >= center:
            add_min = True

        if df["logseqdepth"].max() <= center:
            center = df["logseqdepth"].median()

        if not add_min:
            color_levels = ns.viz._center_colormap(df["logseqdepth"], center)
        else:
            color_levels = ns.viz._center_colormap(
                list(df["logseqdepth"]) + [np.log(1)], center
            )

    # Get heights
    n_rows = len(df["Row"].unique())
    n_cols = len(df["Column"].unique())
    height = int(50 * n_rows)
    width = height * n_cols // n_rows

    # add tooltips
    tooltips = [
        ("amino-acid_substitutions", "@amino-acid_substitutions"),
        ("Alignment Count", "@Alignment Count"),
        ("Alignment Probability", "@Alignment Probability"),
    ]

    def hook(plot, element):
        plot.handles["y_range"].factors = list("HGFEDCBA")
        plot.handles["x_range"].factors = [str(value) for value in range(1, 13)]

    # generate the heatmap
    hm = (
        hv.HeatMap(
            df,
            kdims=["Column", "Row"],
            vdims=[
                "logseqdepth",
                "amino-acid_substitutions",
                "Alignment Count",
                "Alignment Probability",
            ],
        )
        .redim.values(row=np.unique(df["Row"]), Column=np.unique(df["Column"]))
        .opts(
            **opts,
            colorbar=True,
            cmap=seq_depth_cmap,
            height=height,
            width=width,
            line_width=4,
            clipping_colors={"NaN": "#DCDCDC"},
            color_levels=color_levels,
            tools=["hover"],
            colorbar_opts=dict(title="LogSeqDepth", background_fill_alpha=0),
            hooks=[hook],
        )
    )
    # function to bin the alignment frequencies into more relevant groupings
    def bin_align_freq(value):
        if value > 0.95:
            bin_vals = "0.95+"
        elif value <= 0.95 and value > 0.9:
            bin_vals = "0.90-0.95"
        elif value <= 0.9 and value > 0.8:
            bin_vals = "0.80-0.90"
        # anything below 0.8 should really be discarded
        else:
            bin_vals = "<0.80"
        return bin_vals

    # Bin alignment frequencies for easier viz
    bins = ["0.95+", "0.90-0.95", "0.80-0.90", "<0.80"]
    if cmap is None:
        cmap = [cc.bmy[int((1.1 - i) * len(cc.bmy))] for i in [0.95, 0.9, 0.8, 0.4]]
    if "stoplight" in cmap:
        cmap = ["#337D1F", "#94CD35", "#FFC300", "#C62C20"]
    else:
        # Validate colormap
        if not isinstance(cmap, (list, tuple)):
            raise ValueError("cmap argument must be a list or tuple")
        if len(cmap) > 4:
            raise ValueError(
                "cmap argument has too many entries; only 4 should be passed"
            )
    cmap = {bin: color for bin, color in zip(bins, cmap)}

    # apply binning function to the AlignmentFrequency
    df["AlignmentProbabilityBinned"] = df["Alignment Probability"].apply(bin_align_freq)

    # Set up size of the outline boxes
    box_size = height // n_rows * 1.2

    # alignment frequency heatmap for edges around wells
    boxes = hv.Points(
        df.sort_values(["Alignment Probability"], ascending=False),
        ["Column", "Row"],
        "AlignmentProbabilityBinned",
    ).opts(
        **opts,
        marker="square",
        line_color="AlignmentProbabilityBinned",
        line_join="miter",
        cmap=cmap,
        line_width=6,
        fill_alpha=0,
        line_alpha=1,
        legend_position="top",
        size=box_size,
    )

    # Use in apply statement for residue labels
    def split_variant_labels(mutation_string):
        if pd.isna(mutation_string) or mutation_string == '':
            return ''
            
        num_mutations = len(mutation_string.split("_"))

        if num_mutations > 4:
            return str(num_mutations) + " muts"

        mutation_string = mutation_string.replace("?", "")
        new_line_mutations = mutation_string.replace("_", "\n")

        return new_line_mutations

    _df = df.copy()
    _df["Labels"] = _df["amino-acid_substitutions"].apply(split_variant_labels)

    # Set the font size based on if #PARENT# is in a well and num of mutations
    max_num_mutations = _df["Labels"].apply(lambda x: len(x.split("\n"))).max()
    has_parent = "#PARENT#" in _df["Labels"].values

    if max_num_mutations > 3 or has_parent:
        label_fontsize = "8pt"
    else:
        label_fontsize = "8pt"

    labels = hv.Labels(
        _df,
        ["Column", "Row"],
        "Labels",
    ).opts(text_font_size=label_fontsize, **opts, text_color="#000000")
    
    # return formatted final plot
    return (hm * boxes * labels).opts(
        frame_height=550, frame_width=550 * 3 // 2, border=50, show_legend=True
    )



# Main function to return heatmap with or without alignment
def generate_platemaps(
    max_combo_data,
    result_folder,
    cmap=None,
    show_msa=False,
    widget_location="top_left",
):
    """Saves a plate heatmap html generated from from evSeq data.

    Input:
    ------
    max_combo_data: path (str) or DartaFrame
        Path to 'variants.csv' from an LevSeq experiment or
        a pandas DataFrame of that file.
    cmap: list-like or str, default None
        The colormap to use for the well outline indicating alignment
        frequency. If None, defaults to a Plasma-like (colorcet.bmy)
        colormap. If 'stoplight', uses a green-yellow-red colormap (not
        the most colorblind friendly, but highly intuitive). Otherwise
        you may pass any list -like object containing four colors (e.g.,
        ['#337D1F', '#94CD35', '#FFC300', '#C62C20'] for 'stoplight').
    show_msa: bool, default False
    widget_location: string, default 'top_left'
        Location of the widget for navigating plots. Must be one of:
        ['left', 'bottom', 'right', 'top', 'top_left', 'top_right',
        'bottom_left', 'bottom_right', 'left_top', 'left_bottom',
        'right_top', 'right_bottom'].

    Returns:
    --------
    hm_holomap: an interactive Platemap
    unique_plates: list of unique plates in the data,
    plate2barcode: dictionary mapping plate to barcode_plate
    """

    # Convert to dataframe if necessary
    if isinstance(max_combo_data, str):
        max_combo_df = pd.read_csv(max_combo_data)
    else:
        max_combo_df = max_combo_data.copy()

    # Identify unique plates
    unique_plates = max_combo_df.Plate.unique()

    # Create a new DataFrame to modify without affecting the original
    temp_df = max_combo_df.copy()

    # Convert barcode_plate to string and modify its format
    temp_df["barcode_plate"] = temp_df["barcode_plate"].apply(
        lambda x: "RB0" + str(x) if x < 10 else "RB" + str(x)
    )

    # Create a dictionary with unique Plate to modified barcode_plate mapping
    plate2barcode = (
        temp_df[["Plate", "barcode_plate"]]
        .drop_duplicates("Plate")
        .set_index("Plate")["barcode_plate"]
        .to_dict()
    )

    # make logseqdepth column
    max_combo_df["logseqdepth"] = np.log(
        max_combo_df["Alignment Count"],
        out=np.zeros_like(max_combo_df["Alignment Count"], dtype=float),
        where=max_combo_df["Alignment Count"] != 0,
    )
    # Set the center
    center = np.log(10)

    add_min = False
    if max_combo_df["logseqdepth"].min() >= center:
        add_min = True

    # Adjust if it is greater than max of data (avoids ValueError)
    if max_combo_df["logseqdepth"].max() <= center:
        # Adjust the center
        center = max_combo_df["logseqdepth"].median()

    # center colormap
    if not add_min:
        color_levels = ns.viz._center_colormap(max_combo_df["logseqdepth"], center)
    else:
        color_levels = ns.viz._center_colormap(
            list(max_combo_df["logseqdepth"]) + [np.log(1)], center
        )

    # dictionary for storing plots
    hm_dict = {}
    aln_dict = {}

    # Uniform color levels
    for _hm in hm_dict.values():
        _hm.opts({"HeatMap": {"color_levels": color_levels}})

    # Create dropdowns
    plate_selector = pn.widgets.Select(name="Plate", options=list(unique_plates))

    if show_msa:

        well_selector = pn.widgets.Select(name="Well", options=WELL_IDS)

        for plate in unique_plates:

            # Split to just the information of interest
            df = max_combo_df.loc[max_combo_df.Plate == plate].copy()

            # generate a holoviews plot
            hm_dict[plate] = _make_platemap(df, title=plate, cmap=cmap)

        def get_plate_well(plate, well_id):

            """
            Get the platemap and alignment plot for a given well

            Args:
            - plate: str, plate name
            - well_id: str, well ID, ie. 'A1'
            """

            # Split to just the information of interest
            df = max_combo_df.loc[max_combo_df.Plate == plate].copy()

            hm_bokeh = hv.render(hm_dict[plate], backend="bokeh")
            hm_bokeh.toolbar_location = "right"
            hm_bokeh.toolbar.active_drag = None
            hm_bokeh.toolbar.active_scroll = None

            # hm_bokeh = hv.render(
            #         _make_platemap(df, title=plate, cmap=cmap), backend="bokeh"
            #     )

            # hm_bokeh.toolbar_location = 'right'
            # hm_bokeh.toolbar.active_drag = None
            # hm_bokeh.toolbar.active_scroll = None
            # hm_bokeh = hm_dict.get(plate, pn.pane.Markdown("No platemap available for this plate"))

            # Get the row and column
            aln_path = os.path.join(
                result_folder,
                plate,
                plate2barcode[plate],
                well2nb(well_id),
                f"msa_{plate}_{well_id}.fa",
            )

            # plot the alignment using the nc_variant sequence
            aln = plot_sequence_alignment(
                aln_path,
                parent_name=plate,
                well_seq=df[(df["Row"] == well_id[0]) & (df["Column"] == well_id[1:])][
                    "nc_variant"
                ].values[0],
                markdown_title=f"{result_folder} {plate} {plate2barcode[plate]} {well2nb(well_id)} {well_id}",
            )

            # generate a holoviews plot
            return gridplot(
                [[hm_bokeh], [aln]],
                toolbar_location="right",
                sizing_mode="fixed",  # "stretch_width",
            )

        # Function to update the plots based on dropdown selection
        @pn.depends(plate=plate_selector.param.value, well_id=well_selector.param.value)
        def update_plot(plate, well_id):
            """
            Update the plot based on the dropdown selection

            Args:
            - plate: str, plate name
            - well_id: str, well ID, ie. 'A1'
            """
            return get_plate_well(plate, well_id)

        # Layout the dropdowns and the plot
        return pn.Column(pn.Row(plate_selector, well_selector), pn.Column(update_plot))

    else:

        # Generate plots for each plate
        for plate in unique_plates:

            # Split to just the information of interest
            df = max_combo_df.loc[max_combo_df.Plate == plate].copy()

            # generate a holoviews plot
            hm_dict[plate] = _make_platemap(df, title=plate, cmap=cmap)

        # plot from the dictionary
        hm_holomap = hv.HoloMap(hm_dict, kdims=["Plate"])
        # Update widget location
        hv.output(widget_location=widget_location)

        return hm_holomap


########### Functions for the MSA alignment plot ###########
def get_sequence_colors(seqs: list, palette="viridis") -> list[str]:

    """
    Get colors for a sequence without parent seq highlighting differences

    Args:
    - seqs: list of sequences
    - palette: str, name of the color palette

    Returns:
    - list: colors for each nucleotide
    """

    aas = deepcopy(ALL_AAS)
    aas.append("-")
    aas.append("X")

    pal = plt.colormaps[palette]
    pal = [mpl.colors.to_hex(i) for i in pal(np.linspace(0, 1, 20))]
    pal.append("white")
    pal.append("gray")

    pcolors = {i: j for i, j in zip(aas, pal)}
    nuc = [i for s in list(seqs) for i in s]

    try:
        colors = [NUC_COLOR_DICT[i] for i in nuc]
    except:
        colors = [pcolors[i] for i in nuc]

    return colors


def get_sequence_diff_colorNseq(seqs: list, seq_ids: list, parent_seq: str) -> tuple:

    """
    Get colors and nucleotides for input sequences highlighting differences from parent

    Args:
    - seqs: str, list of sequences
    - seq_ids: str, list of sequence ids
    - parent_seq: str, parent sequence to compare against

    Returns:
    - block_colors: list of colors for each nucleotide highlighting differences
    - nuc_colors: list of colors for each nucleotide
    - nucs: list of nucleotides highlighting differences
    - spacers: list of spacers (for plotting)
    """
    # color for the highlighted nuc block over text
    block_colors = []
    # color for the nuc text to be annotated
    nuc_textcolors = []
    # init nuc to annotate
    diff_nucs = []
    # parent nuc or spacer annotation
    text_annot = []

    for seq, seq_id in zip(seqs, seq_ids):
        if seq_id == "parent":
            for p in list(parent_seq):
                block_colors.append(match_color)
                nuc_textcolors.append(NUC_COLOR_DICT[p])
                diff_nucs.append(" ")
                text_annot.append(" ")
        else:
            for n, p in zip(list(seq), list(parent_seq)):
                if n != p:
                    block_colors.append(NUC_COLOR_DICT[n])
                    diff_nucs.append(n)
                    if n == "-":
                        text_annot.append("-")
                    else:
                        text_annot.append(" ")
                else:
                    block_colors.append(match_color)
                    diff_nucs.append(" ")
                    text_annot.append(" ")

                nuc_textcolors.append("gray")

    return block_colors, nuc_textcolors, diff_nucs, text_annot


def get_cons(aln: MultipleSeqAlignment) -> list[float]:

    """
    Get conservation values from alignment

    Args:
    - aln: MultipleSeqAlignment, input alignment

    Returns:
    - list: conservation values
    """

    x = []
    l = len(aln)
    for i in range(aln.get_alignment_length()):
        a = aln[:, i]
        res = Counter(a)
        del res["-"]
        x.append(max(res.values()) / l)
    return x


def get_cons_seq(aln: MultipleSeqAlignment, ifdeg: bool = True) -> str:

    """
    Ger consensus sequence from alignment

    Args:
    - aln: MultipleSeqAlignment, input alignment
    - ifdeg: bool, if True, return degenerate consensus

    Returns:
    - str: consensus sequences
    """

    alignment = aln.alignment
    motif = Motif("ACGT", alignment)

    if ifdeg:
        return motif.degenerate_consensus
    else:
        return motif.consensus


def get_cons_diff_colorNseq(cons_seq: str, parent_seq: str) -> tuple:

    """
    Get consensus sequence highlighting differences from parent

    Args:
    - cons_seq: str, consensus sequence
    - parent_seq: str, parent sequence

    Returns:
    - colors: list, colors for each nucleotide highlighting differences
    - cons_seq_diff: list, nucleotides highlighting differences
    """

    colors = []
    cons_seq_diff = []
    for n, p in zip(list(cons_seq), list(parent_seq)):
        if n != p:
            cons_seq_diff.append(n)
            if n in NUC_COLOR_DICT.keys():
                colors.append(NUC_COLOR_DICT[n])
            else:
                colors.append("#f2f2f2")
        else:
            cons_seq_diff.append(" ")
            colors.append(match_color)

    return colors, cons_seq_diff


def aggregate_gray_blocks(x_vals: list, y_vals: list, colors: list, text: list):

    """
    Aggregate gray blocks in the MSA alignment plot
    to reduce the number of elements for plotting

    Args:
    - x_vals: list, x values
    - y_vals: list, y values
    - colors: list, colors
    - text: list, text

    Returns:
    - aggregated_x: list, aggregated x values
    - aggregated_y: list, aggregated y values
    - aggregated_width: list, aggregated width values
    - aggregated_height: list, aggregated height values
    - aggregated_colors: list, aggregated colors
    - aggregated_text: list, aggregated text
    """

    aggregated_x = []
    aggregated_y = []
    aggregated_width = []
    aggregated_height = []
    aggregated_colors = []
    aggregated_text = []

    current_x_start = None
    current_y = None
    current_width = 0

    for i, (x, y, color, t) in enumerate(zip(x_vals, y_vals, colors, text)):
        if color == "gray":
            # Start or continue aggregating gray blocks
            if current_x_start is None:
                current_x_start = x
                current_y = y
                current_width = 1
            else:
                if y == current_y:
                    # Continue aggregating in the same row
                    current_width += 1
                else:
                    # Row changed, finalize the current block
                    aggregated_x.append(current_x_start + current_width / 2)
                    aggregated_y.append(current_y)
                    aggregated_width.append(current_width)
                    aggregated_height.append(1)
                    aggregated_colors.append("gray")
                    aggregated_text.append("")

                    # Start a new gray block for the new row
                    current_x_start = x
                    current_y = y
                    current_width = 1
        else:
            # Add the current gray block if it exists
            if current_x_start is not None:
                aggregated_x.append(current_x_start + current_width / 2)
                aggregated_y.append(current_y)
                aggregated_width.append(current_width)
                aggregated_height.append(1)
                aggregated_colors.append("gray")
                aggregated_text.append("")

                # Reset aggregation variables
                current_x_start = None
                current_width = 0

            # Add the non-gray block as it is
            aggregated_x.append(x)
            aggregated_y.append(y)
            aggregated_width.append(1)
            aggregated_height.append(1)
            aggregated_colors.append(color)
            aggregated_text.append(t)

    # Add any remaining aggregated gray block
    if current_x_start is not None:
        aggregated_x.append(current_x_start + current_width / 2)
        aggregated_y.append(current_y)
        aggregated_width.append(current_width)
        aggregated_height.append(1)
        aggregated_colors.append("gray")
        aggregated_text.append("")

    return (
        aggregated_x,
        aggregated_y,
        aggregated_width,
        aggregated_height,
        aggregated_colors,
        aggregated_text,
    )


def aggregate_conservation(x_vals, heights, colors):

    """
    Aggregate gray blocks in the conservation plot
    to reduce the number of elements for plotting

    Args:
    - x_vals: list, x values
    - heights: list, heights
    - colors: list, colors

    Returns:
    - aggregated_x: list, aggregated x values
    - aggregated_height: list, aggregated heights
    - aggregated_colors: list, aggregated colors
    """

    aggregated_x = []
    aggregated_height = []
    aggregated_colors = []
    current_x_start = None
    current_width = 0

    for i, (x, height, color) in enumerate(zip(x_vals, heights, colors)):
        if color == "gray" and height == 2:
            # Start or continue aggregating gray blocks with height = 2
            if current_x_start is None:
                current_x_start = x
                current_width = 1
            else:
                current_width += 1
        else:
            # Add the current gray block if it exists
            if current_x_start is not None:
                aggregated_x.append(current_x_start + current_width / 2)
                aggregated_height.append(2)
                aggregated_colors.append("gray")
                current_x_start = None
                current_width = 0

            # Add the non-gray or non-height-2 block as it is
            aggregated_x.append(x)
            aggregated_height.append(height)
            aggregated_colors.append(color)

    # Append the final aggregated gray block if any
    if current_x_start is not None:
        aggregated_x.append(current_x_start + current_width / 2)
        aggregated_height.append(2)
        aggregated_colors.append("gray")

    return aggregated_x, aggregated_height, aggregated_colors


def plot_empty(msg="", plot_width=1000, plot_height=200) -> figure:
    """
    Return an empty bokeh plot with optional text displayed

    Args:
    - msg: str, message to display
    - plot_width: int, width of the plot
    - plot_height: int, height of the plot

    Returns:
    - figure: bokeh plot
    """

    p = figure(
        width=plot_width,
        height=plot_height,
        x_range=(0, 1),
        y_range=(0, 2),
        sizing_mode="fixed",  # "stretch_width",
        # output_backend="webgl"
    )

    text = Label(x=0.3, y=1, text=msg)
    p.add_layout(text)
    p.grid.visible = False
    p.xaxis.visible = False
    p.yaxis.visible = False
    return p


def plot_sequence_alignment(
    aln_path: str,
    parent_name: str = "parent",
    well_seq: str = "",
    markdown_title: str = "Multiple sequence alignment",
    fontsize: str = "4pt",
    plot_width: int = 1000,
    sizing_mode: str = "fixed",  # "stretch_width",
    palette: str = "viridis",
    row_height: float = 8,
) -> figure:

    """
    Plot sequence alignment

    Args:
    - aln_path: str, path to the alignment file
    - parent_name: str, name of the parent sequence
    - well_seq: str, sequence of the well
    - markdown_title: str, title of the plot
    - fontsize: str, fontsize of the text
    - plot_width: int, width of the plot
    - sizing_mode: str, sizing mode of the plot
    - palette: str, color palette
    - row_height: float, height of the row

    Returns:
    - figure: bokeh plot
    """

    # get text from markdown
    msa_title = Div(
        text=f"""
    {markdown_title}
    """
    )

    # check if alignment file exists
    if not os.path.exists(aln_path):
        p = plot_empty("Alignment file not found", plot_width)
        return gridplot(
            [[msa_title], [p]],
            toolbar_location=None,
            sizing_mode=sizing_mode,
        )

    # read in alignment
    aln = AlignIO.read(aln_path, "fasta")

    seqs = [rec.seq for rec in (aln)]
    ids = [rec.id for rec in aln]
    rev_ids = list(reversed(ids))

    seq_len = len(seqs[0])
    numb_seq = len(seqs)

    # get parent sequence
    parent_seq = None
    for rec in aln:
        if rec.id in ["parent", "Parent", parent_name]:
            parent_seq = rec.seq
            break

    # check if alignment has at least two sequences
    if len(seqs) <= 1:
        p = plot_empty("Alignment plot needs at least two sequences", plot_width)
        return gridplot(
            [[msa_title], [p]],
            toolbar_location=None,
            sizing_mode=sizing_mode,
        )

    seq_nucs = [i for s in list(seqs) for i in s]

    # get colors for the alignment
    if parent_seq == None:
        block_colors = get_sequence_colors(seqs=seqs, palette=palette)
        text = seq_nucs
        text_colors = "black"
    else:
        parent_nucs = [i for i in parent_seq] * numb_seq
        (
            block_colors,
            text_colors,
            diff_nucs,
            text,
        ) = get_sequence_diff_colorNseq(seqs=seqs, seq_ids=ids, parent_seq=parent_seq)

    # get conservation values
    cons = get_cons(aln)

    # get consensus sequence
    if well_seq == "":
        cons_seq = get_cons_seq(aln)
    else:
        cons_seq = well_seq

    cons_nucs = [i for i in cons_seq] * numb_seq

    # coords of the plot
    x = np.arange(1, seq_len + 1)
    y = np.arange(0, numb_seq, 1)
    xx, yy = np.meshgrid(x, y)
    gx = xx.ravel()
    gy = yy.flatten()
    recty = gy + 0.5

    # Apply transformation to IDs
    y_flipped = numb_seq - gy - 1
    ids_repeated = [ids[yi] for yi in y_flipped]
    rev_ids_repeated = [rev_ids[yi] for yi in y_flipped]

    # set up msa source
    msa_source = ColumnDataSource(
        dict(
            x=gx,
            y=y_flipped,
            ids=ids_repeated,
            rev_ids=rev_ids_repeated,
            seq_nucs=seq_nucs,
            parent_nucs=parent_nucs,
            cons_nucs=cons_nucs,
            recty=numb_seq - recty,
            block_colors=block_colors,
            text=text,
            text_colors=text_colors,
        )
    )

    # define the plot size and range
    plot_height = len(seqs) * row_height + 20
    x_range = Range1d(0, seq_len + 1, bounds="auto")

    # Aggregating gray blocks
    agg_x, agg_y, agg_width, agg_height, agg_colors, agg_text = aggregate_gray_blocks(
        msa_source.data["x"],
        msa_source.data["y"],
        msa_source.data["block_colors"],
        msa_source.data["text"],
    )

    # Create the updated ColumnDataSource with aggregated gray blocks
    msa_source_aggregated = ColumnDataSource(
        dict(
            x=agg_x,
            y=agg_y,
            width=agg_width,
            height=agg_height,
            fill_color=agg_colors,
            text=agg_text,
        )
    )

    # Create a new figure
    p_aln = figure(
        title=None,
        width=plot_width,
        height=plot_height,
        x_range=(0, max(agg_x) + 1),
        y_range=(-0.5, max(agg_y) + 0.5),
    )

    # Add aggregated rectangles (blocks of sequences)
    rect_glyph = Rect(
        x="x",
        y="y",
        width="width",
        height="height",
        fill_color="fill_color",
        line_color=None,
    )

    # Add the rectangles to the plot
    p_aln.add_glyph(msa_source_aggregated, rect_glyph)

    p_aln.grid.visible = False
    p_aln.yaxis.visible = False
    p_aln.yaxis.major_label_text_font_size = "0pt"
    p_aln.yaxis.minor_tick_line_width = 0
    p_aln.yaxis.major_tick_line_width = 0

    # conservation plot
    cons_colors, cons_text = get_cons_diff_colorNseq(
        cons_seq=cons_seq, parent_seq=parent_seq
    )

    cons_source = ColumnDataSource(
        dict(
            x=x,
            cons=cons,
            cons_height=[2 * c for c in cons],
            cons_colors=cons_colors,
            parent_nucs=list(parent_seq),
            cons_nucs=list(cons_seq),
            # cons_text=cons_text,
        )
    )

    # Example usage with the conservation data
    agg_cons_x, agg_cons_height, agg_cons_colors = aggregate_conservation(
        cons_source.data["x"],
        cons_source.data["cons_height"],
        cons_source.data["cons_colors"],
    )

    # Create the updated ColumnDataSource with aggregated gray blocks
    cons_source_aggregated = ColumnDataSource(
        dict(
            x=agg_cons_x,
            height=agg_cons_height,
            fill_color=agg_cons_colors,
        )
    )

    # Create a new figure for the aggregated conservation plot
    p_cons = figure(
        title=None,
        width=plot_width,
        height=60,
        x_range=p_aln.x_range,
        y_range=(Range1d(0, 1)),
        # tools=[cons_hover],
        # output_backend="webgl"
    )

    # Add aggregated rectangles (conservation bars)
    cons_rects = Rect(
        x="x", y=0, width=1, height="height", fill_color="fill_color", line_color=None
    )

    # cons_text = Text(
    #     x="x",
    #     y=0,
    #     text="cons_text",
    #     text_align="center",
    #     text_color="black",
    #     text_font_size=fontsize,
    # )

    p_cons.add_glyph(cons_source_aggregated, cons_rects)

    # Adding the legend items for each nucleotide from NUC_COLOR_DICT
    legend_items = []
    for nucleotide, color in dict(list(NUC_COLOR_DICT.items())[:5]).items():
        legend_source = ColumnDataSource(
            dict(x=[0], y=[0], width=[0], height=[0], fill_color=[color])
        )
        rect = p_cons.rect(
            x="x",
            y="y",
            width="width",
            height="height",
            fill_color="fill_color",
            line_color="#F5F5F5",
            source=legend_source,
            visible=False,
        )
        legend_items.append(LegendItem(label=nucleotide, renderers=[rect]))

    # Create and add the horizontal legend to the plot
    legend = Legend(items=legend_items, location="top_center", orientation="horizontal")
    p_cons.add_layout(legend, "above")  # Place legend above the plot

    # Remove the legend box and background
    p_cons.legend.border_line_color = None  # Removes the border of the legend
    p_cons.legend.background_fill_alpha = 0  # Makes the legend background transparent

    # Tighten the legend layout
    p_cons.legend.padding = 0
    p_cons.legend.spacing = 0

    # Hide the x-axis labels and keep the y-axis visible
    p_cons.xaxis.visible = False
    p_cons.yaxis.visible = True
    p_cons.yaxis.ticker = [1]
    p_cons.yaxis.axis_label = "Alignment conservation values"
    p_cons.yaxis.axis_label_orientation = "horizontal"

    p_cons.grid.visible = False
    p_cons.background_fill_color = "white"

    return gridplot(
        [[msa_title], [p_cons], [p_aln]],
        toolbar_location=None,
        sizing_mode=sizing_mode,
    )
