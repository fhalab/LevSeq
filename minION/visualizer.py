from minION.util.IO_processor import IO_processor
import numpy as np
import plotly.express as px


def barcode_score(summary, nbins = 100 , unclassified = True):
    """The Barcode score is a quality metric for each read. A read is only considered to be valid, if the barcode
    score is above a threshold. For lower quality, the reads are stored in unclassified folder
    Input: 
        - Summary file (.txt) by default. For each demultiplex run, a summary file is generated. For evSeq-minION, one summary file should exist
        for reverse barcode (# Plates). In addition, within each reverse barcode folder, a summary file should exist summarising the demultiplexing for 
        forward barcode (# Samples).
        - nbins: Integer. Number of bins for the histogram
        - unclassified: Boolean. If True, the unclassified reads are included in the histogram. If False, the unclassified reads are excluded.
    Output:
        - Histogram of Barcode Score Distribution """
    

    fig = px.histogram(summary, x="barcode_score", nbins=100, 
                       title="Barcode Score Distribution", 
                       labels={
                                'barcode_score': 'Barcode Score', 
                               'count': 'Frequency'})
    
    fig.update_traces(  marker_color='red',
                        marker_line_color='black', 
                        marker_line_width=1.5, 
                        opacity=0.5)
    
    fig.update_layout(  width=600, 
                        height=400)
    
    return fig

