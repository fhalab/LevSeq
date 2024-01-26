# Global parameters for minIon 


# Parameter for Demultiplexing
SCORE_MATRIX =  {
        ('A', 'A'): 96,  ('A', 'C'): -316, ('A', 'G'): -192, ('A', 'T'): -369, ('A', 'N'): 0,
        ('C', 'A'): -316,('C', 'C'): 100,  ('C', 'G'): -352, ('C', 'T'): -295, ('C', 'N'): 0,
        ('G', 'A'): -192,('G', 'C'): -352, ('G', 'G'): 98,   ('G', 'T'): -329, ('G', 'N'): 0,
        ('T', 'A'): -369,('T', 'C'): -295, ('T', 'G'): -329, ('T', 'T'): 100,  ('T', 'N'): 0,
        ('N', 'A'): 0,   ('N', 'C'): 0,    ('N', 'G'): 0,    ('N', 'T'): 0,    ('N', 'N'): 0,
    } # Adapted from Guppy Barcoder

SW_ALIGN_PARAMS = {
    "start_gap1": 40,
    "end_gap1": 40,
    "open_gap1": 0,
    "extend_gap1": -40,
    "start_gap2": 40,
    "end_gap2": 40,
    "open_gap2": -160,
    "extend_gap2": -160,
    "min_score_barcode_front": 60.0,
    "front_window_size": 150,
    "rear_window_size": 150,
} # Adapted from Guppy Barcoder


# Path where the minION reads are stored 
MINKNOW_PATH = "/var/lib/minknow/data/"

# Defaul target folder names from ONT
DEFAULT_TARGETS = {"Not_basecalled": ["pod5"], "Basecalled": ["fastq_pass"]}

# Results Folder
RESULTS_FOLDER = "results"




# Dorado Models
DORADO_MODELS = {   "fast": "dna_r10.4.1_e8.2_400bps_fast@v4.2.0",  # Fast inference (least accurate)
                    "hac": "dna_r10.4.1_e8.2_400bps_hac@v4.2.0",    # Balanced inference (balanced accuracy and speed)
                    "sup": "dna_r10.4.1_e8.2_400bps_sup@v4.2.0"     # Most accurate inference (slowest)
                }

# Medaka Models

MEDAKA_MODELS = { "default" : "r1041_e82_400bps_sup_v4.2.0"}

# Barcodes
BARCODES = {    "Barcode-kit-FBC" : "FBC-sim",
                "Barcode-kit-RBC" : "RBC"}

# Codons
CODONS = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
    'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
    }


