"""Checkpoint-based peptide classification and score exports."""

import os

from model.MLP_f1_strong_dropout import MLP

import csv

import pandas as pd

import torch

import numpy as np

import pandas

def mlp_prediction_loaded(epitopes,output_dir, epi, f_string, hla_string, model):
    print("Starting the binding prediction of " + str(len(epitopes)) + " Epitopes")
    d = {0: 'non', 1: 'weak', 2: 'strong'}
    d = {'0': 'non', '1': 'weak', '2': 'strong'}
    d = {'non':0, 'weak':1, 'strong' :2}
    full_epitopes_scores = pandas.DataFrame(
        columns=["Epitope", "Mismatch", "HLA", "Binding", "Start", "End", "Uniprot_gene_id",
                 "Gene_name",
                 "ENST", "ENSG", "RefSeq_id",
                 "Affinity", "Processing_score", "Presentation_score", "Presentation_percentile",
                 "Tap_score", "Tap10_score", "Tap04_score", "Cleavage_score", "Processing_score_flanks",
                 "Presentation_score_flanks",
                 "Presentation_percentile_flanks", "Q-Score", "Weblinks", "Uniprot_db", "Expression",
                 "Flags_annotation"])

    for epitope in epitopes:  # create DF with all attribute values for all epitopes
        for hit in epitope.protein_hits:
            test = hit[21]
            l = np.append(test, [hit[7], hit[8], hit[9], hit[10], hit[15], hit[16], hit[17], hit[6], hit[12], hit[13], hit[14]])

            pred = model(torch.FloatTensor([l]))
            pred_probabilities = torch.nn.functional.softmax(pred, dim=1)
            prediction = torch.argmax(pred_probabilities)
            y = str(prediction.item())
            hit.extend([k for k, v in d.items() if v == prediction.item()])
            full_epitopes_scores = pandas.concat([full_epitopes_scores, pandas.DataFrame(
                {"Epitope": [hit[0]], "HLA": hla_string[1:], "Affinity": [hit[7]], "Binding": [hit[22]],
                 "Processing_score": [hit[8]],
                 "Presentation_score": [hit[9]],
                 "Presentation_percentile": [hit[10]],
                 "Tap_score": [hit[15]], "Tap10_score": [hit[16]],
                 "Tap04_score": [hit[17]], "Cleavage_score": [hit[6]],
                 "Processing_score_flanks": [hit[12]], "Presentation_score_flanks": [hit[13]],
                 "Presentation_percentile_flanks": [hit[14]], "Q-Score": [hit[18]], "Mismatch": [hit[3]],
                 "Start": [hit[1]], "End": [hit[2]], "Uniprot_gene_id": [epitope.uniprot_gene_id],
                 "Gene_name": [epitope.external_gene_name], "ENSG": [epitope.ensembl_gene_id],
                 "ENST": [epitope.ensembl_transcript_id], "RefSeq_id": [epitope.refSeq_id],
                 "Weblinks": [epitope.links], "Uniprot_db": [epitope.uniprot_db], "Expression": [epitope.expression],
                 "Flags_annotation": [epitope.flags]})])

    epi_write = epi[:12]
    full_epitopes_scores.to_csv(path_or_buf=output_dir + '/' + epi_write + f_string + hla_string+'_epitope_search_scores_mlp.tsv', sep='\t', index=False, mode='w')
    df = pandas.read_csv(output_dir + '/' + epi_write + f_string + hla_string+'_epitope_search_scores_mlp.tsv', sep='\t')
    df = df.sort_values(['Presentation_score', 'Mismatch'], ascending=[False, True])
    df.to_csv(output_dir + '/' + epi_write + f_string + hla_string+'_epitope_search_scores_mlp.tsv', mode='w', sep='\t', index=False)

def mlp_prediction_reloaded(epitopes, output_dir, epi, f_string, hla_string, model):
    print(f"Starting the binding prediction of {len(epitopes)} Epitopes")

    label_map = {0: 'non', 1: 'weak', 2: 'strong'}

    feature_list = []
    hit_refs = []
    meta_records = []

    extra_indices = [7, 8, 9, 10, 15, 16, 17, 6, 12, 13, 14]

    for epitope in epitopes:
        for hit in epitope.protein_hits:
            embedding_vector = hit[21]
            extra_features = [hit[idx] for idx in extra_indices]

            combined = np.hstack([embedding_vector, extra_features])

            feature_list.append(combined)
            hit_refs.append(hit)

            meta_records.append({
                "Epitope": hit[0],
                "HLA": hla_string[1:] if len(hla_string) > 1 else hla_string,
                "Affinity": hit[7],
                "Processing_score": hit[8],
                "Presentation_score": hit[9],
                "Presentation_percentile": hit[10],
                "Tap_score": hit[15],
                "Tap10_score": hit[16],
                "Tap04_score": hit[17],
                "Cleavage_score": hit[6],
                "Processing_score_flanks": hit[12],
                "Presentation_score_flanks": hit[13],
                "Presentation_percentile_flanks": hit[14],
                "Q-Score": hit[18],
                "Mismatch": hit[3],
                "Start": hit[1],
                "End": hit[2],
                "Uniprot_gene_id": epitope.uniprot_gene_id,
                "Gene_name": epitope.external_gene_name,
                "ENSG": epitope.ensembl_gene_id,
                "ENST": epitope.ensembl_transcript_id,
                "RefSeq_id": epitope.refSeq_id,
                "Weblinks": epitope.links,
                "Uniprot_db": epitope.uniprot_db,
                "Expression": epitope.expression,
                "Flags_annotation": epitope.flags
            })

    if not feature_list:
        print("No epitope hits found to predict.")
        return

    features_array = np.asarray(feature_list, dtype=np.float32)
    device = next(model.parameters()).device
    inputs = torch.from_numpy(features_array).to(device)

    predictions = []
    batch_size = 4096

    with torch.no_grad():
        for i in range(0, len(inputs), batch_size):
            batch_inputs = inputs[i:i + batch_size]
            logits = model(batch_inputs)
            probs = torch.nn.functional.softmax(logits, dim=1)
            preds = torch.argmax(probs, dim=1).cpu().tolist()
            predictions.extend(preds)

    for hit, record, pred_idx in zip(hit_refs, meta_records, predictions):
        binding_label = label_map.get(pred_idx, 'non')
        hit.append(binding_label)  # Appends to hit[22]
        record["Binding"] = binding_label

    columns_order = [
        "Epitope", "Mismatch", "HLA", "Binding", "Start", "End", "Uniprot_gene_id",
        "Gene_name", "ENST", "ENSG", "RefSeq_id", "Affinity", "Processing_score",
        "Presentation_score", "Presentation_percentile", "Tap_score", "Tap10_score",
        "Tap04_score", "Cleavage_score", "Processing_score_flanks",
        "Presentation_score_flanks", "Presentation_percentile_flanks", "Q-Score",
        "Weblinks", "Uniprot_db", "Expression", "Flags_annotation"
    ]

    full_epitopes_scores = pd.DataFrame(meta_records)[columns_order]

    full_epitopes_scores.sort_values(
        by=['Presentation_score', 'Mismatch'],
        ascending=[False, True],
        inplace=True
    )

    epi_write = epi[:12]
    out_filepath = os.path.join(output_dir, f"{epi_write}{f_string}{hla_string}_epitope_search_scores_mlp.tsv")
    full_epitopes_scores.to_csv(out_filepath, sep='\t', index=False)

def mlp_prediction_multi_in_one_loaded_bench_create_trainingset(epitopes,output_dir, epi, f_string, hla_string, model):
    d = {0: 'non', 1: 'weak', 2: 'strong'}
    d = {'0': 'non', '1': 'weak', '2': 'strong'}
    d = {'non':0, 'weak':1, 'strong' :2}
    full_epitopes_scores = pandas.DataFrame(
        columns=["Epitope", "Mismatch", "HLA", "Binding", "Start", "End", "Uniprot_gene_id",
                 "Gene_name",
                 "ENST", "ENSG", "RefSeq_id",
                 "Affinity", "Processing_score", "Presentation_score", "Presentation_percentile",
                 "Tap_score", "Tap10_score", "Tap04_score", "Cleavage_score", "Processing_score_flanks",
                 "Presentation_score_flanks",
                 "Presentation_percentile_flanks", "Q-Score", "Weblinks", "Uniprot_db", "Expression",
                 "Flags_annotation", "Prediction"])

    for epitope in epitopes:  # create DF with all attribute values for all epitopes
        for hit in epitope.protein_hits:
            test = hit[21]
            l = np.append(test, [hit[7], hit[8], hit[9], hit[10], hit[15], hit[16], hit[17], hit[6], hit[12], hit[13], hit[14]])
            pred = model(torch.FloatTensor([l]))
            print(pred)
            pred_probabilities = torch.nn.functional.softmax(pred, dim=1)
            prediction = torch.argmax(pred_probabilities)
            print(pred_probabilities)
            y = str(prediction.item())
            print(y)
            hit.extend([k for k, v in d.items() if v == prediction.item()])
            print(hit[22])
            full_epitopes_scores = pandas.concat([full_epitopes_scores, pandas.DataFrame(
                {"Epitope": [hit[0]], "HLA": hla_string[1:], "Affinity": [hit[7]], "Binding": [hit[22]],
                 "Processing_score": [hit[8]],
                 "Presentation_score": [hit[9]],
                 "Presentation_percentile": [hit[10]],
                 "Tap_score": [hit[15]], "Tap10_score": [hit[16]],
                 "Tap04_score": [hit[17]], "Cleavage_score": [hit[6]],
                 "Processing_score_flanks": [hit[12]], "Presentation_score_flanks": [hit[13]],
                 "Presentation_percentile_flanks": [hit[14]], "Q-Score": [hit[18]], "Mismatch": [hit[3]], "Start": [hit[1]],
                 "End": [hit[2]], "Uniprot_gene_id": [epitope.uniprot_gene_id],
                 "Gene_name": [epitope.external_gene_name], "ENSG": [epitope.ensembl_gene_id],
                 "ENST": [epitope.ensembl_transcript_id], "RefSeq_id": [epitope.refSeq_id],
                 "Weblinks": [epitope.links], "Uniprot_db": [epitope.uniprot_db], "Expression": [epitope.expression],
                 "Flags_annotation": [epitope.flags], "Prediction": str(pred_probabilities)})])

    full_epitopes_scores.to_csv(path_or_buf=output_dir + '/Mutation_epitope_search_scores_mlp.tsv', sep='\t',
                                    index=False, mode='a')

    row = [hit[0], hla_string[1:]] + l.tolist()

    with open(output_dir + '/bench_scores_mlp.tsv', 'a', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(row)
