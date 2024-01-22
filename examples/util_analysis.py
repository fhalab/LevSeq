# Functions 

# Pandas releated

def rewrite_Plate_name(entry):
    entry = entry[-1]

    return int(entry)

def rewrite_Well_name(entry):
    num = int(entry[1:])
    new_entry = entry[0] + str(num)

    return new_entry

def single_plate_annotation(entry):
    row = ["A", "B", "C", "D", "E", "F", "G", "H"]
    new_well_name = row[int(entry["Plate"]) - 1] + entry["Well"][1:]
    entry["Well"] = new_well_name
    return entry

def variant_combo(entry):
    if entry in ["#PARENT#", "NA"]:
        return entry

    elif isinstance(entry, float) and math.isnan(entry):
        return entry

    else:
        return "_".join(entry)       
        

def transform(entry):
    # If the entry is either DEAD or PARENT, return empty lists
    if entry in ["#DEAD#", "#PARENT#"]:
        return [], [], []
    else:
        # Split the entry based on the underscore
        mutations = entry.split("_")
        parent_combo = []
        positions = []
        new_aa = []
        for mutation in mutations:
            # Append the first character to parent_combo, the middle to positions, and the last to new_aa

            #Skip if parent is equal to new_aa
            if mutation[0] == mutation[-1]:
                continue

            elif mutation.find("DEL") != -1:
                parent_combo.append(mutation[0])
                positions.append(mutation[1:-3])
                new_aa.append("DEL")
            
            else:
                parent_combo.append(mutation[0])
                positions.append(mutation[1:-1])
                new_aa.append(mutation[-1])

        return parent_combo, positions, new_aa

def transform_ref(entry):

    parent_combo = []
    positions = []
    new_aa = []

    try:

        for variant in entry:

            if variant in ["#DEAD#", "#PARENT#"]:
                return ["-"], ["-"], ["-"]

            elif variant.find("DEL") != -1:
                parent_combo.append(variant[0])
                positions.append(variant[1:-3])
                new_aa.append("DEL")
            
            else:
                parent_combo.append(variant[0])
                positions.append(variant[1:-1])
                new_aa.append(variant[-1])
        return parent_combo, positions, new_aa
    
    except:
        return ["NA"], ["NA"], ["NA"]
        
def substract_by_index(positions, index):
    new_positions = []
    for pos in positions:
        new_positions.append(str(int(pos) - index))
    
    return new_positions

def split_variant_combo(entry):
    
    if entry in ["#DEAD#", "#PARENT#"]:
        return [entry]
    
    
    elif entry == "NA":
        return ["#DEAD#"]

    
    elif isinstance(entry, str) and entry.startswith("wt"):
        return ["#PARENT#"]

    
    elif isinstance(entry, float) and math.isnan(entry):
        return ["#DEAD#"]
        
    
    else:
        
        return str(entry).split("_")
    
def combine_variant_combo(entry):

    if isinstance(entry, float) and math.isnan(entry):
        return "#DEAD#"

    elif "#PARENT#" in str(entry):
        return "#PARENT#"

    else:
        return "_".join(entry)

def compare_mutations(row, mutation_only = True):

    if isinstance(row['Variant'], float):
        return pd.Series([0, 0, [], []], index=['Correct', 'Missed', 'CorrectMutations', 'MissedMutations'])

    if mutation_only and row['Variant_Combo_Ref'] in [["#DEAD#"], ["#PARENT#"]]:
        return pd.Series([0, 0, [], []], index=['Correct', 'Missed', 'CorrectMutations', 'MissedMutations'])

    
    mutations1 = set(row['Variant_Combo_Ref'])
    mutations2 = set(row['Variant'])
    
    correct_mutations = mutations1.intersection(mutations2)
    missed_mutations = mutations1.difference(mutations2)
    
    return pd.Series([len(correct_mutations), len(missed_mutations), list(correct_mutations), list(missed_mutations)], index=['Correct', 'Missed', 'CorrectMutations', 'MissedMutations'])

def compare_only_Parent(row):

    if isinstance(row['Variant'], float):
        return pd.Series([0, 0, [], []], index=['Correct', 'Missed', 'CorrectMutations', 'MissedMutations'])
    
    if row['Variant_Combo_Ref'] in [["#PARENT#"]]:
        mutations1 = set(row['Variant_Combo_Ref'])
        mutations2 = set(row['Variant'])
        
        correct_mutations = mutations1.intersection(mutations2)
        missed_mutations = mutations1.difference(mutations2)
        return pd.Series([len(correct_mutations), len(missed_mutations), list(correct_mutations), list(missed_mutations)], index=['Correct', 'Missed', 'CorrectMutations', 'MissedMutations'])

    else:
         return pd.Series([0, 0, [], []], index=['Correct', 'Missed', 'CorrectMutations', 'MissedMutations'])

def count_mutations(row):
    if isinstance(row['Variant'], float):
        return 0
    
    if row['Variant_Combo_Ref'] in [["#DEAD#"], ["#PARENT#"]]:
        return 0
    
    else:
        return len(row['Variant_Combo_Ref'])
    
def count_parent(row):
    if isinstance(row['Variant'], float):
        return 0
    
    if row['Variant_Combo_Ref'] in [["#PARENT#"]]:
        return 1
    
    else:
        return 0

def count_all(row):

    if isinstance(row['Variant_Combo_Ref'], float):
        return 0
    
    elif row['Variant_Combo_Ref'] in [["#PARENT#"]]:
        return 1
    
    else:
        return len(row['Variant_Combo_Ref'])


# Function to create VariantCombo
def create_variant_combo(row):
    variant_combos = []
    for parent, pos, aa in zip(row['ParentCombo'], row['Positions'], row['NewAA']):
        if parent and pos and aa:
            variant_combos.append(f"{parent}{pos}{aa}")
    return variant_combos if variant_combos else ["#PARENT#"] if not row['ParentCombo'] else ["#DEAD#"]
    

def process_gzip_files(filename):
    
    # Open the gzip-compressed FASTQ file for reading
    with gzip.open(filename, 'rt') as handle:
        for record in SeqIO.parse(handle, "fastq"):
            # Access various attributes of the FASTQ record
            print("ID:", record.id)
            print("Sequence:", record.seq)
            print("Quality Scores:", record.letter_annotations["phred_quality"])
            print("\n-----\n")


def filter_single_chunk(records, min_length, max_length):
    """Filter a chunk of FASTQ records based on length."""
    return [record for record in records if min_length <= len(record.seq) <= max_length]

def filter_sequences(filename, min_length=650, max_length=1000, output_file="filtered_sequences.fastq.gz", n_workers=4):
    """Filter basecalled sequences based on length in parallel.
    
    Args:
        filename (str): .fastq.gz file name
        min_length (int): minimum length of sequence
        max_length (int): maximum length of sequence
        output_file (str): output file name for filtered sequences
        n_workers (int): number of processes for parallel processing

    Returns:
        str: Fastq file name with filtered sequences
    """
    
    # Break records into chunks for parallel processing
    chunk_size = 1000  # You can adjust this value based on your requirements
    all_records = []
    with gzip.open(filename, 'rt') as handle:
        all_records = list(SeqIO.parse(handle, "fastq"))
    chunks = [all_records[i:i + chunk_size] for i in range(0, len(all_records), chunk_size)]
    
    filtered_records = []
    
    # Parallel processing of chunks
    with ProcessPoolExecutor(max_workers=n_workers) as executor:
        futures = [executor.submit(filter_single_chunk, chunk, min_length, max_length) for chunk in chunks]
        for future in futures:
            filtered_records.extend(future.result())

    # Write filtered sequences to the output file
    with open(output_file, 'wt') as output_handle:
        SeqIO.write(filtered_records, output_handle, "fastq")

    return output_file

def get_seq_dist(filename):
    """
    Get Sequence distribution of .fastq files of the experiment

    Args:
        filename (str): .fastq file name
    
    Returns:
        seq_dist (list): list of length of each sequence
    """
    
    seq_dist = []
    with gzip.open(filename, 'rt') as handle:
        for record in SeqIO.parse(handle, "fastq"):
            seq_dist.append(len(record.seq))
    return seq_dist

def adjust_variant(variant, padding_start):
    """
    Adjusts the variant position based on the padding length. After variant calling,
    the position of the variant is based on the reference sequence.
    """

    # Check if variant is NaN
    if isinstance(variant, float) and math.isnan(variant):
        return "NA"

    # Convert variant to string for string-specific checks
    variant_str = str(variant)

    # Check for specific string patterns
    if variant_str == "NA" or "#PARENT#" in variant_str:
        return variant_str

    # Process the variant if it's a valid string
    else:
        variants = variant_str.split('_')
        adjusted_variants = []

        for v in variants:
            # Find the position number using regular expression
            match = re.search(r'([A-Za-z]+)(\d+)([A-Za-z]+)', v)
            if match:
                refAA, pos, newAA = match.groups()
                adjusted_pos = max(int(pos) - padding_start, 1)  
                adjusted_variants.append(f"{refAA}{adjusted_pos}{newAA}")

        return '_'.join(adjusted_variants)


def get_variant_df_guppy(demultiplexer_path, template_fasta):
    # Load barcode dictionaries
    
    #barcode_dicts = IO_processor.get_barcode_dict(Path(demultiplexer_path), "NB", "RB") Skipped for bacode_simulator

    # Initialize variants dictionary
    variants = {
        "RBC": [], "FBC": [], "Position": [],
        "Variant": [], "Quality-Score": [], "Sequence": []
    }

    # Get template sequence
    template = analyser.get_template_sequence(template_fasta)
    template = template.upper()
    

    summary = analyser.read_summary_file(demultiplexer_path)
    n_counts = summary.groupby(["RBC","FBC"])["FBC"].value_counts().reset_index()

    # # Process barcode dictionaries to calculate n_counts
    # n_counts = pd.DataFrame()
    # for barcode_id, barcode_dict in barcode_dicts.items():
    #     n_counts_rbc = analyser.read_summary_file(barcode_id)["barcode_arrangement"].value_counts().reset_index()
    #     rbc_basename = os.path.basename(barcode_id)
    #     n_counts_rbc = n_counts_rbc[n_counts_rbc["barcode_arrangement"] != "unclassified"]
    #     n_counts_rbc.rename(columns={"barcode_arrangement" : "FBC", "count" : "n_counts"}, inplace=True)
    #     n_counts_rbc["RBC"] = rbc_basename
    #     n_counts = pd.concat([n_counts, n_counts_rbc], axis=0)

    # Process each barcode
    for barcode_id, barcode_dict in barcode_dicts.items():
        rbc = os.path.basename(barcode_id)
        for front_barcode in barcode_dict:
            fbc = os.path.basename(front_barcode)
            front_barcode_path = os.path.join(front_barcode)


            #fasta_file = os.path.join(front_barcode_path, "consensus", "consensus.fastq")
            fasta_file = os.path.join(front_barcode_path, "consensus.fastq")

            # Check if consensus file exists
            if not os.path.exists(fasta_file):
                print(f"Consensus file in {front_barcode_path} does not exist, skipping {fbc}")
                continue

            try:
                consensus = analyser.get_consensus_sequence(fasta_file, True)
            except:
                print(f"Skipping {rbc}/{fbc}")
                continue

            variants["Sequence"].append(consensus["Sequence"][0])
            nn_variants = analyser.call_variant_nn(template, consensus["Sequence"][0], consensus["Quality-Score"][0])
            variants["RBC"].append(rbc)
            variants["FBC"].append(fbc)
            variants["Position"].append(nn_variants["Position"])
            variants["Variant"].append(nn_variants["Variant"])
            variants["Quality-Score"].append(nn_variants["Quality-Score"])

    # Create DataFrame from variants
    variant_df = analyser.rename_barcode_guppy(pd.DataFrame(variants))

    # Merge with template DataFrame
    variant_template_df = analyser.template_df(barcode_dicts, rowwise=False)
    variant_df = variant_df.merge(variant_template_df, on=["Plate", "Well"], how="right")

    return variant_df

def get_contingency(row):
    cont = {"TP": 0, "FP": 0, "TN": 0, "FN": 0}
    
    # Check if Variant is NaN or contains "#DEAD#"
    if isinstance(row['Variant'], float) or ("#DEAD#" in str(row['Variant'])) or (str(row['Variant']) == "['']"):
        if "PARENT" in str(row['Variant_Combo_Ref']):
            cont["FN"] += 1
        else:
            cont["FP"] += 1

    
    # Check if Variant contains "#PARENT#"
    elif "#PARENT#" in str(row['Variant']):
        if "#PARENT#" in str(row['Variant_Combo_Ref']):
            cont["TN"] += 1
        else:
            #num_var = len(str(row['Variant_Combo_Ref']).split("_"))
            cont["FP"] += 1
    
    else:
        called = False
        for var in row['Variant']:
            if var in str(row['Variant_Combo_Ref']):
                called = True
            else:
                called = False
        
        if called:
            cont["TP"] += 1
        else:
            cont["FP"] += 1
    
    return pd.Series(cont)

def calculate_metrics(variant_df):

    variant_df = variant_df.join(variant_df.apply(get_contingency, axis=1))

    sensitivity = variant_df["TP"].sum() / (variant_df["TP"].sum() + variant_df["FN"].sum())
    specificity = variant_df["TN"].sum() / (variant_df["TN"].sum() + variant_df["FP"].sum())
    precision = variant_df["TP"].sum() / (variant_df["TP"].sum() + variant_df["FP"].sum())
    accuracy = (variant_df["TP"].sum() + variant_df["TN"].sum()) / (variant_df["TP"].sum() + variant_df["TN"].sum() + variant_df["FP"].sum() + variant_df["FN"].sum())
    F1 = 2 * (precision * sensitivity) / (precision + sensitivity)

    # Printing the results
    print(f"Sensitivity: {sensitivity}")
    print(f"Specificity: {specificity}")
    print(f"Precision: {precision}")
    print(f"Accuracy: {accuracy}")
    print(f"F1: {F1}")

    # Creating a dictionary of results
    results = {"Sensitivity": [sensitivity], "Specificity": [specificity], "Precision": [precision], "Accuracy": [accuracy], "F1": [F1]}

    return results


def analyze_fastq_data(fastq_files_path):
    """
    Analyze basecalled reads in fastq files
    """
    fastq_files = glob.glob(fastq_files_path)
    names = []
    lengths = []
    qualities_mean = []

    read_data = {"Name" : [], "Length": [], "Mean_Quality": []}

    for fastq_file in fastq_files:
        with gzip.open(fastq_file, "rt") as handle:
            for record in SeqIO.parse(handle, "fastq"):
                if 750 <= len(record.seq) <= 1000:
                    read_data["Name"].append(record.id)
                    read_data["Length"].append(len(record.seq))
                    read_data["Mean_Quality"].append(np.mean(record.letter_annotations["phred_quality"]))

    if lengths:
        stats = {
            'n_reads': len(lengths),
            'mean_length': np.mean(lengths),
            'median_length': np.median(lengths),
            'max_length': np.max(lengths),
            'min_length': np.min(lengths),
            'cv': np.std(lengths) / np.mean(lengths),
            'upper_outlier_length': np.quantile(lengths, 0.95),
            'lower_outlier_length': np.quantile(lengths, 0.05)
        }
    else:
        stats = {'error': 'No valid reads found.'}

    return read_data, stats

def get_base_quality_scores(fastq_files_path, min_length, max_length, names = True):
    """
    Get base quality scores from FASTQ files.
    """
    fastq_files = glob.glob(fastq_files_path)
    if names:
        read_id = []
    
    qualities = []

    for fastq_file in fastq_files:
        with gzip.open(fastq_file, "rt") as handle:
            for record in SeqIO.parse(handle, "fastq"):
                if min_length <= len(record.seq) <= max_length:
                    if names: 
                        read_id.append(record.id)
                        qualities.append(record.letter_annotations["phred_quality"])

    return {"Name": read_id, "Qualities": qualities} if names else qualities


def plot_barcode_score(summary, bin_size, threshold, save_path = None, column = "RB", barcode_score_name = "barcode_score"):
    """
    Plot the barcode scores for each bin
    """
    fig, ax = plt.subplots(figsize=(12, 7))
    bins = np.arange(0, summary['barcode_score'].max() + bin_size, bin_size)
    sns.set_style("whitegrid")

    sns.histplot(data=summary[summary[column] >= threshold], x= barcode_score_name, bins=bins, color='lightgreen', edgecolor = "black" , ax=ax)

    sns.histplot(data=summary[summary[column] < threshold], x= barcode_score_name, bins=bins, color='lightcoral',  edgecolor = "black" , ax=ax)

    plt.axvline(x=threshold, color='r', linestyle='--', linewidth=2)
    plt.tick_params(size=14, labelsize=14)
    plt.xlabel("Barcode Score", size=16)
    plt.ylabel("# of Reads", size=16)
    plt.savefig(save_path, dpi = 600, bbox_inches='tight')
    plt.show()

def plot_quality_scores(qualities, save_path):
    min_length = min(len(q) for q in qualities)
    truncated_qualities = [q[:min_length] for q in qualities]

    quality_array = np.array(truncated_qualities)

    mean_qualities = np.mean(quality_array, axis=0)
    std_dev = np.std(quality_array, axis=0)

    # Create the plot
    plt.figure(figsize=(12, 7))
    plt.plot(mean_qualities, label='Mean Quality Score', linewidth = 2.5, color = "palevioletred")
    plt.fill_between(range(min_length), mean_qualities - std_dev, mean_qualities + std_dev, color='grey', alpha=0.5, label='Standard Deviation')
    plt.xlabel('Read Position')
    plt.ylabel('Phred Quality Score')
    #plt.title('Mean Quality Score at Each Position with Standard Deviation')
    plt.legend()
    plt.savefig(save_path, dpi = 600, bbox_inches='tight')
    plt.show()


def plot_mean_distribution(df, save_path, num_bins):


    bin_ranges = pd.cut(df['Read Length'], bins=num_bins)
    grouped = df.groupby(bin_ranges)['Phred Quality Score']
    mean_qualities = grouped.mean()
    std_qualities = grouped.std()
    #cv_qualities = grouped.std() / grouped.mean()
    bin_centers = mean_qualities.index.map(lambda x: x.mid)

    fig, ax1 = plt.subplots(figsize=(12, 7))

    # Plotting the histogram
    color_hist = 'tab:blue'
    ax1.set_xlabel('Read Length', size=16)
    ax1.set_ylabel('Count', color="black", size=16)
    ax1.hist(df['Read Length'], bins=num_bins, alpha=0.8, color="lightgrey", edgecolor = "black")
    ax1.tick_params(axis='y', labelcolor="black", labelsize=14)
    ax1.tick_params(axis='x', labelsize=14)

    # Creating the second axes for the mean quality score
    ax2 = ax1.twinx()

    # Adding standard deviation with fill_between
    ax2.fill_between(bin_centers, 
                    mean_qualities - std_qualities, 
                    mean_qualities + std_qualities, 
                    color='palevioletred', alpha=0.3)

    ax2.scatter(bin_centers, mean_qualities, color="palevioletred", s=50, zorder=3)
    ax2.plot(bin_centers, mean_qualities, color="palevioletred", zorder=3)

    color_line = 'tab:red'
    ax2.set_ylabel('Mean Phred Quality Score', color="palevioletred", size=16)
    ax2.tick_params(axis='y', labelcolor="palevioletred", labelsize=14)

    fig.tight_layout()
    plt.savefig(save_path, dpi = 300, bbox_inches='tight')
    plt.show()
