from minION.util.IO_processor import get_barcode_dict
from minION import analyser
from minION import consensus
import subprocess
import pysam
from collections import Counter



def call_and_filter_vcf(input_path, reference, allele_frequency, alignment_name = "alignment.bam", depth : int = 4000):

    """Call variants after getting alignments from mini_align. For now, bcftools is being used to filter variants with a given allele frequence. 
    This will be replaced in future custom scripts.
    
    Args:
        - input_path: Path to the folder containing alignment.bam
        - reference: Path to the reference genome
        - allele_frequency: Allele frequency threshold for filtering variants
    Return:
        - None """

    prompt_pileup = f"bcftools mpileup -d {depth} -Ou -f {reference}  {input_path}/{alignment_name} > {input_path}/pileup.bcf"

    prompt_call = f"bcftools call -mv -Ob -o {input_path}/raw_variants.bcf {input_path}/pileup.bcf"

    prompt_view = f"bcftools view -i 'INFO/AF>{allele_frequency}' -Ob -o {input_path}/filtered_variants.vcf {input_path}/raw_variants.bcf"


    subprocess.run(prompt_pileup, shell=True)

    subprocess.run(prompt_call, shell=True)

    subprocess.run(prompt_view, shell=True)

    print(f"Variant calling and filtering completed. Output saved to {input_path}/raw_variants_python.bcf")


def extract_positions_from_vcf(input_path):
    """Extract positions from a VCF file."""

    vcf_file = f"{input_path}/filtered_variants.vcf"

    positions = []

    # Open the VCF file using pysam 
    vcf = pysam.VariantFile(vcf_file)

    for record in vcf:
        # Get the position of the variant
        position = record.pos
        positions.append(position)
    
    # Close the VCF file
    vcf.close()

    return positions



def generate_heatmap_data(bam_file, chrom, start, end):
    data = []
    all_bases = set(['A', 'T', 'C', 'G', '-'])
    for position in range(start, end + 1):
        base_counts = get_base_counts_at_position(bam_file, chrom, position)
        col = [base_counts.get(base, 0) for base in all_bases]
        data.append(col)
    df = pd.DataFrame(data, index=range(start, end + 1), columns=list(all_bases))
    # Order as in all_bases
    df = df[['A', 'T', 'C', 'G', '-']]
    return df