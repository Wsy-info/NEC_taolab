import pandas as pd
import sys
import os

pd.set_option('display.max_columns', 100)
# Define our base directory for the analysis
os.chdir('../NEC/sc_RNA-seq/NEW/all_cell/cellphonedb')

print(sys.version)

out_path = '../NEC/sc_RNA-seq/NEW/all_cell/cellphonedb'
cpdb_file_path = '../reference/cellphonedb/cpdb_tutorial/db/database/v5.0.0/cellphonedb_10_14_2025_162716.zip'
meta_file_path = '../NEC/sc_RNA-seq/NEW/all_cell/cellphonedb/cellphonedb_meta.txt'
counts_file_path = '../NEC/round2_NEC.harmony_singlet.h5ad'


from cellphonedb.src.core.methods import cpdb_statistical_analysis_method

cpdb_results = cpdb_statistical_analysis_method.call(
    cpdb_file_path = cpdb_file_path,                 # mandatory: CellphoneDB database zip file.
    meta_file_path = meta_file_path,                 # mandatory: tsv file defining barcodes to cell label.
    counts_file_path = counts_file_path,             # mandatory: normalized count matrix - a path to the counts file, or an in-memory AnnData object
    counts_data = 'hgnc_symbol',                     # defines the gene annotation in counts matrix.
    # active_tfs_file_path = active_tf_path,         # optional: defines cell types and their active TFs.
    microenvs_file_path = None,                      # optional (default: None): defines cells per microenvironment.
    score_interactions = True,                       # optional: whether to score interactions or not.
    iterations = 1000,                               # denotes the number of shufflings performed in the analysis.
    threshold = 0.1,                                 # defines the min % of cells expressing a gene for this to be employed in the analysis.
    threads = 10,                                    # number of threads to use in the analysis.
    debug_seed = 0,                                  # debug randome seed. To disable >=0.
    result_precision = 3,                            # Sets the rounding for the mean values in significan_means.
    pvalue = 0.05,                                   # P-value threshold to employ for significance.
    subsampling = False,                             # To enable subsampling the data (geometri sketching).
    subsampling_log = False,                         # (mandatory) enable subsampling log1p for non log-transformed data inputs.
    subsampling_num_pc = 100,                        # Number of componets to subsample via geometric skectching (dafault: 100).
    subsampling_num_cells = 30000,                   # Number of cells to subsample (integer) (default: 1/3 of the dataset).
    separator = '|',                                 # Sets the string to employ to separate cells in the results dataframes "cellA|CellB".
    debug = False,                                   # Saves all intermediate tables employed during the analysis in pkl format.
    output_path = out_path,                          # Path to save results.
    output_suffix = None                             # Replaces the timestamp in the output files by a user defined string in the  (default: None).
    )
