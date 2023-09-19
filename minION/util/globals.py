# Global parameters for minIon 



# Path where the minION reads are stored 
MINKNOW_PATH = "/var/lib/minknow/data/"

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
BARCODES = {    "Barcode-kit" : "MINION-BARCODES-TRIMMED",
                "Barcode-kit-RBC" : "MINION-BARCODES-RBC-TRIMMED",
                "Barcode-kit-MASKED" : "MINION-BARCODES-MMM"
                }

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


