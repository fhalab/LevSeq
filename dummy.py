from minION import variantcaller
from pathlib import Path
demultiplex_path = Path("/home/emre/minION_results/20240112-RL-8Plates-FLongle-2_sup/demultiplexed/")

variantcaller.check_demultiplexing(demultiplex_path)