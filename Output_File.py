"""Write expression, log-sum, and tissue-score tables."""

import os

def write_Expression_to_file(expression_tables, epi , id2, output_dir):
    epi_write = epi[:12]
    for df in expression_tables:

        df[1].iloc[:,0:].to_csv(f"{output_dir}/{epi_write}_EX_{str(os.path.basename(df[0]))}", mode='w', sep='\t', index=False)

def write_LogSums_to_file(Log_sums, epi , id2, output_dir):
    epi_write = epi[:12]
    for df in Log_sums:

        with open(f"{output_dir}/{epi_write}_LOG_{str(os.path.basename(df[0]))}", 'w') as f:

            df[1].to_csv(f, mode='a', header=f.tell() == 0, index=False, sep='\t')

def write_TissueSums_to_file(Tissue_sums, epi , id2, output_dir):
    epi_write = epi[:12]
    for df in Tissue_sums:

        with open(f"{output_dir}/{epi_write}_ICR_{str(os.path.basename(df[0]))}", 'w') as f:
            df[1].to_csv(f, mode='a', header=f.tell() == 0, index=False, sep='\t')
