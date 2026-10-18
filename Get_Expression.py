"""Join peptide scores with tissue expression databases."""

import pandas as pd
import pandas

import os

from tqdm import tqdm

import re

def search_expression_databases_new_v4_collapsed(epitopes, databases, output_dir, epi, f_string, hla_string):
    result_tables = list()

    def listdir_fullpath(d):
        return [os.path.join(d, f) for f in os.listdir(d)]

    def _grep(x, list_to_grep):
        """takes a string (x) and checks whether any string from a given
           list of strings (list_to_grep) exists in `x`"""
        for text in list_to_grep:
            if text != "":
                if text.lower() in x.lower():
                    return True
        return False

    epi_write = epi[:12]
    df = pandas.read_csv(output_dir + '/' + epi_write + f_string + hla_string + '_epitope_search_scores_mlp.tsv', sep='\t')
    df = df.fillna('')
    print("Expressionswerte werden gesucht")
    dbs = listdir_fullpath(databases)

    for db in tqdm(dbs):  # create Dataframes from .tsv Databases, [0] db name [1] filepath
        db_name = db.split(".")[0].split("\\")[-1]
        table = pandas.read_table(db)
        table.columns.values[0] = "id"

        p = pandas.DataFrame()
        grep_id = "Uniprot_gene_id"
        if _grep(table.iloc[1, 0], ["ENSG"]):  #changed: from _grep(table.iloc[1, 0], "ENSG")
            grep_id = "ENSG"
        elif _grep(table.iloc[1, 0], ["ENST"]):
            grep_id = "ENST"
        p = pandas.merge(df, table, left_on=grep_id, right_on="id")  # right_on="id")
        mask = df[grep_id].apply(_grep, list_to_grep=',')
        for index, row in df[mask].iterrows():
            s = re.sub(",+", ",", (str(row.ENST) + "," + str(row.ENSG) + "," + str(row.Uniprot_gene_id)))
            s = re.sub("^,", "", s)
            s = re.sub(",$", "|", s)
            s = re.sub(",", "|", s)
            g = s.split("|")
            mask2 = table.id.apply(_grep, list_to_grep=g)
            if any(mask2):
                t = table[mask2].iloc[:, 1:].mean()
                i = ",".join((table[mask2].iloc[:, 0].unique()).tolist())
                t = pandas.concat([pandas.Series(i), t])
                t = t.rename({0: "id"})
                toadd = row.to_frame().T.reset_index(drop=True).join(t.to_frame().T)  # pandas.concat([p, toadd])
                p = pandas.concat([p, toadd])
        p.reset_index(inplace=True)  # Resets the index, makes factor a column
        p = p.iloc[:, 1:]
        result_tables.append((db_name, p))

    return result_tables

def search_expression_databases_new_v5_collapsed(epitopes, databases, output_dir, epi, f_string, hla_string):
    result_tables = []

    epi_write = epi[:12]
    tsv_path = os.path.join(output_dir, f"{epi_write}{f_string}{hla_string}_epitope_search_scores_mlp.tsv")
    df = pd.read_csv(tsv_path, sep='\t').fillna('')
    print("Expressionswerte werden gesucht")

    dbs = [os.path.join(databases, f) for f in os.listdir(databases) if f.endswith(('.tsv', '.txt', '.csv'))]

    for db in tqdm(dbs):
        db_name = os.path.splitext(os.path.basename(db))[0]
        table = pd.read_table(db)
        table.rename(columns={table.columns[0]: "id"}, inplace=True)
        table["id"] = table["id"].astype(str)

        sample_id = str(table.iloc[0]["id"]) if len(table) > 0 else ""
        if "ENSG" in sample_id:
            grep_id = "ENSG"
        elif "ENST" in sample_id:
            grep_id = "ENST"
        else:
            grep_id = "Uniprot_gene_id"

        merged_direct = pd.merge(df, table, left_on=grep_id, right_on="id")

        numeric_cols = table.columns[1:]

        table_grouped = table.groupby("id")[numeric_cols].mean()

        has_comma = df[grep_id].str.contains(',', na=False)
        comma_df = df[has_comma]

        fallback_rows = []

        if not comma_df.empty:
            existing_db_ids = set(table_grouped.index)

            for idx, row in comma_df.iterrows():
                raw_ids = f"{row.get('ENST', '')},{row.get('ENSG', '')},{row.get('Uniprot_gene_id', '')}"
                candidate_ids = [i.strip() for i in re.split(r'[,|]+', raw_ids) if i.strip()]

                matched_ids = [cid for cid in candidate_ids if cid in existing_db_ids]

                if matched_ids:
                    mean_vals = table_grouped.loc[matched_ids].mean(axis=0)

                    row_dict = row.to_dict()
                    row_dict["id"] = ",".join(matched_ids)
                    for col in numeric_cols:
                        row_dict[col] = mean_vals[col]

                    fallback_rows.append(row_dict)

        if fallback_rows:
            fallback_df = pd.DataFrame(fallback_rows)
            p = pd.concat([merged_direct, fallback_df], ignore_index=True)
        else:
            p = merged_direct

        result_tables.append((db_name, p))

    return result_tables

def write_scores_mut(epitopes, output_dir, epi,hla_string):

    full_epitopes_scores = pd.DataFrame(columns =["epitope","Distance", "Start Position", "End Position", "Uniprot_gene_id", "external_gene_name", "ENST","ENSG", "refSeq_id",
                                                          "affinity","processing_score","presentation_score","presentation_percentile",
                      "tap", "tap10", "tap_new","pepsickle","processing_score_flanks","presentation_score_flanks",
                                                      "presentation_percentile_flanks","Q-Score","weblinks","uniprot_db", "expression", "flags"])
    full_epitopes_scores = pd.DataFrame(
        columns=["Epitope", "Mismatch", "HLA", "Start", "End", "Uniprot_gene_id",
                 "Gene_name",
                 "ENST", "ENSG", "RefSeq_id",
                 "Affinity", "Processing_score", "Presentation_score", "Presentation_percentile",
                 "Tap_score", "Tap10_score", "Tap04_score", "Cleavage_score", "Processing_score_flanks",
                 "Presentation_score_flanks",
                 "Presentation_percentile_flanks", "Q-Score", "Weblinks", "Uniprot_db", "Expression",
                 "Flags_annotation"])

    for epitope in epitopes:  # create DF with all attribute values for all epitopes
        for hit in epitope.protein_hits:
            full_epitopes_scores = pd.concat([full_epitopes_scores, pd.DataFrame(
                {"Epitope": [hit[0]], "HLA": hla_string[1:], "Affinity": [hit[7]],
                 "Processing_score": [hit[8]],
                 "Presentation_score": [hit[9]],
                 "Presentation_percentile": [hit[10]],
                 "Tap_score": [hit[15]], "Tap10_score": [hit[16]],
                 "Tap04_score": [hit[17]], "Cleavage_score": [hit[6]],
                 "Processing_score_flanks": [hit[12]], "Presentation_score_flanks": [hit[13]],
                 "Presentation_percentile_flanks": [hit[14]], "Q-Score": [hit[18]], "Mismatch": [hit[3]],
                 "Start": [hit[1]],
                 "End": [hit[2]], "Uniprot_gene_id": [epitope.uniprot_gene_id],
                 "Gene_name": [epitope.external_gene_name], "ENSG": [epitope.ensembl_gene_id],
                 "ENST": [epitope.ensembl_transcript_id], "RefSeq_id": [epitope.refSeq_id],
                 "Weblinks": [epitope.links], "Uniprot_db": [epitope.uniprot_db], "Expression": [epitope.expression],
                 "Flags_annotation": [epitope.flags]})])

    print ("Leaving out the print out for now")
