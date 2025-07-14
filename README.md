# Variant Sequencing with Nanopore (LevSeq)

LevSeq provides a streamlined pipeline for sequencing and analyzing genetic variants using Oxford Nanopore technology. In directed evolution experiments, LevSeq enables sequencing of every variant, enhancing data insight and creating datasets suitable for AI/ML methods. Sequence variants can be generated within a day at an extremely low cost.

![Figure 1: LevSeq Workflow](manuscript/figures/LevSeq_Figure-1.jpeg)
Figure 1: Overview of the LevSeq variant sequencing workflow using Nanopore technology. This diagram illustrates the key steps in the process, from sample preparation to data analysis and visualization.

## Quick Start

### Docker Installation (Recommended)

1. Install Docker: [https://docs.docker.com/engine/install/](https://docs.docker.com/engine/install/)
2. Pull the appropriate image:
   ```bash
   # For Linux/Windows x86 systems:
   docker pull yueminglong/levseq:levseq-1.4-x86
   
   # For Mac M-series chips (M1, M2, M3, M4):
   docker pull yueminglong/levseq:levseq-1.4-arm64
   ```
3. Run LevSeq:
   ```bash
   docker run --rm -v "/full/path/to/data:/levseq_results" yueminglong/levseq:levseq-1.4-arm64 my_experiment levseq_results/ levseq_results/ref.csv
   ```
4. Connect function data to your sequence data
   ```bash
   docker run --rm -v "/full/path/to/data:/levseq_results" yueminglong/levseq:levseq-1.4-arm64 my_experiment levseq_results/ levseq_results/ref.csv --fitness_files "levseq_results/20250712_epPCR_Q06714_37.csv,levseq_results/20250712_epPCR_Q06714_39.csv,levseq_results/20250712_epPCR_Q06714_40.csv" --smiles 'O=P(OC1=CC=CC=C1)(OC2=CC=CC=C2)OC3=CC=CC=C3>>O=P(O)(OC4=CC=CC=C4)OC5=CC=CC=C5' --compound dPPi --variant_df "levseq_results/visualization_partial.csv"
   ```
### Pip Installation (Mac/Linux only)

**IMPORTANT**: On Mac M-series chips (M1-M4), gcc 13 and 14 are **REQUIRED**:
```bash
brew install gcc@13 gcc@14
```

1. Create and activate conda environment:
   ```bash
   conda create --name levseq python=3.12 -y
   conda activate levseq
   ```

2. Install dependencies:
   ```bash
   conda install -c bioconda -c conda-forge samtools minimap2
   ```

3. Install LevSeq:
   ```bash
   pip install levseq
   ```

4. Run LevSeq:
   ```bash
   levseq my_experiment /path/to/data/ /path/to/ref.csv
   ```

5. Combine function data:
   ```bash
   levseq my_experiment /path/to/data/ /path/to/ref.csv  "LCMS_file_{barcode1}.csv,LCMS_file_{barcode2}.csv," --smiles 'reaction_smiles_string' --compound "name_of_compound_in_LCMS_file" --variant_df "visualization_partial.csv"
   ```

Note for function data we currently expect a LCMS file e.g. with the columns: 
- `Sample Vial Number` (corresponding to the well that the sample was from). 
- `Area` (which becomes fitness value). 
- `Compound Name` which is the name of the compound we filter for that is passed as a parameter.
- The last `_X.csv` needs to be the barcode number to match that sample to your plate e.g. if you ran LevSeq with barcode 33 for plate 2 you need to have `_33.csv` for the fitness file for plate 2. e.g. `some_fitnes_for_plate_2_33.csv`.


## Data and Visualization

- **Test Data**: Sample data is available on Zenodo [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.13694463.svg)](https://doi.org/10.5281/zenodo.13694463)
- **Visualization Tool**: A web application is available at [https://levseqdb.streamlit.app/](https://levseqdb.streamlit.app/) - simply upload your LevSeq output and LCMS results
- **Self-hosted Solution**: You can deploy your own instance using our [LevSeq_db repository](https://github.com/fhalab/LevSeq_db)

## Reference File Format (ref.csv)

Your reference CSV file must contain the following columns:

| barcode_plate | name   | refseq    |
|---------------|--------|-----------|
| 33            | Q97A76 | ATGCGC... |

For oligopool experiments (multiple proteins per plate), use:

| barcode_plate | name   | refseq    |
|---------------|--------|-----------|
| 33            | Q97A76 | ATGCGCAAG |
| 33            | P96084 | ATGGATCA  |
| 34            | P46209 | ATGGGGCAA |
| 34            | Q60336 | ATGGGGCC  |

## Command Line Arguments

### Required Arguments
1. **name**: Name of the experiment (output folder)
2. **path**: Location of basecalled fastq files
3. **summary**: Path to reference CSV file

### Optional Arguments
- `--skip_demultiplexing`: Skip the demultiplexing step
- `--skip_variantcalling`: Skip the variant calling step
- `--output`: Custom save location (defaults to current directory)
- `--show_msa`: Show multiple sequence alignment for each well
- `--oligopool`: Process data as oligopool experiment

## Step-by-Step Tutorial

1. **Prepare your sequencing data**:
   - Your fastq files should be in a directory structure similar to Nanopore's output
   - Prepare a reference CSV file with barcode plates, sample names, and reference sequences

2. **Run LevSeq**:
   ```bash
   # Via Docker
   docker run --rm -v "/path/to/data:/levseq_results" yueminglong/levseq:levseq-1.4-arm64 my_experiment levseq_results/ levseq_results/ref.csv
   
   # Via pip
   levseq my_experiment /path/to/data/ /path/to/ref.csv
   ```

3. **Analyze results**:
   - Output includes variant data (CSV) and interactive visualizations (HTML)
   - Upload results to the LevSeq visualization tool for further analysis

## Experimental Setup

For the wet lab protocol:
- Refer to the [wiki](https://github.com/fhalab/LevSeq/wiki/Experimental-protocols)
- See the methods section of [our paper](https://pubs.acs.org/doi/10.1021/acssynbio.4c00625)
- Order forward and reverse primers compatible with your plasmid
- Install Oxford Nanopore's software for basecalling if needed

## Additional Resources

- **Example Notebook**: See `example/Example.ipynb` for a walkthrough
- **Advanced Usage**: See the [manuscript notebook](https://github.com/fhalab/LevSeq/blob/main/manuscript/notebooks/epPCR_10plates.ipynb)
- **Troubleshooting**: See our [computational protocols wiki](https://github.com/fhalab/LevSeq/wiki/Computational-protocols)

## Citing LevSeq

If you find LevSeq useful, please cite our paper:

```bibtex
@article{long2024levseq,
  title={LevSeq: Rapid Generation of Sequence-Function Data for Directed Evolution and Machine Learning},
  author={Long, Yueming and Mora, Ariane and Li, Francesca-Zhoufan and Gürsoy, Emre and Johnston, Kadina E and Arnold, Frances H},
  journal={ACS Synthetic Biology},
  year={2024},
  publisher={American Chemical Society}
}
```

## Contact

Leave a feature request in the issues or reach us via [email](mailto:levseqdb@gmail.com).