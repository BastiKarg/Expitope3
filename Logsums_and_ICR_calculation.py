"""Tissue-weighted log sums and ICR calculations."""

import pandas as pd

import numpy as np

def calculate_log_sums_and_icr_v2(expression_data,qthresh,k,ouput_dir):
    tissue_sums = list()
    log_sums = list()
    normalized_mismatch_weight = [0.95, 0.0477, 0.00248, 0.000253, 0.000209, 0.000209] ## Fehler im Expitope2 Paper?
    normalized_mismatch_weight = [0.95, 0.0477, 0.002498, 0.00027353, 0.00023509] ## Fehler im Expitope2 Paper?
    normalized_mismatch_weight = [1, 0.0477, 0.002498, 0.00027353, 0.00023509] ## Fehler im Expitope2 Paper?

    tissue_weights = pd.read_csv('tissue_weights_exp4.unix_lowercase_unique', header=None, sep="\t")

    for df in expression_data:  ### BASTI was wird denn hier gemacht?
        tissue_profile = pd.DataFrame()
        log_sum_table = pd.DataFrame()

        if df[1].empty:
            continue

        for m in range(k+1):

            data = df[1][(df[1]["Mismatch"] == m)].iloc[:, 28:]  # setzen auf 20??   XXX zurücksetzten ### so ein Dreck, falsche Vorarbeit... evtl alle Tabellen aendern
            data = data.drop_duplicates()

            non_floats = []
            for col in data:

                if data[col].dtypes != "float64":
                    non_floats.append(col)

            sum = np.sum(data).to_frame()

            numeric_cols = [col for col in sum if sum[col].dtype.kind != 'O']
            numeric_cols
            sum += 1
            log_sum_table = pd.concat([log_sum_table, sum.transpose().applymap(np.log10)],ignore_index=True)

            sum = sum.applymap(np.log10)
            sum *= normalized_mismatch_weight[m]

            total_weight_of_tissues = 0  # sum (float(tissue_weights.loc[tissue_weights[0] == t, 1]) for t in sum.index )

            x = [(tissue_weights.loc[tissue_weights[0] == t, 1]) for t in sum.index]

            for i in range(0, len(x)):

                print(df[0])
                print(sum.index[i])
                total_weight_of_tissues += (float(x[i].drop(columns=1)))

            icr = 0
            if total_weight_of_tissues == 0:
                tissue_profile = pd.concat([tissue_profile, pd.concat([pd.DataFrame({'ICR': [icr]}).transpose(), sum]).transpose()], ignore_index=True)
            else:

                for index, row in sum.iterrows():
                    icr += float(tissue_weights.loc[tissue_weights[0] == index, 1]) *(row.values)

                icr = icr / total_weight_of_tissues
                if m > 0:
                    icr += tissue_profile.get('ICR').iloc[-1]
                tissue_profile = pd.concat([tissue_profile, pd.concat([pd.DataFrame({'ICR': [icr]}).transpose(), sum]).transpose()], ignore_index=True)

        tissue_sums.append((df[0], tissue_profile))
        log_sums.append((df[0], log_sum_table))

    return log_sums, tissue_sums
