"""Precomputed and on-demand proteasomal cleavage scores."""

import pandas as pd

import numpy as np

from tqdm import tqdm

import pepsickle

def run_pepsickle_new_loaded(epitopes, pepsickle_file_loaded, thresh_pcleav=0.5):
    print("Lesen der Proteosomal Cleavage aktuell loaded")
    from pprint import pprint
    list_of_ids = []
    for e in epitopes:
        list_of_ids.append(e.database_id)
    print(list_of_ids)
    pep_res = pepsickle_file_loaded
    pep_res = pep_res[pep_res['protein_id'].isin(list_of_ids)]
    print(pep_res)
    for e in tqdm(epitopes):
        pep_epitope = pep_res[(pep_res['protein_id']).astype(int) == int(e.database_id)]

        for protein_hit in e.protein_hits:
            probs = list(pep_epitope.iloc[protein_hit[1]:protein_hit[2] + 1, 1])  ## BASTI TODO: evtl andere Spalte nehmen

            if probs and len(probs) > 0:
                 pC = probs.pop()
            else:
                pC = 0.0

            pCon = np.prod([1 - prob for prob in probs if prob > thresh_pcleav])
            protein_hit.append(pCon * pC)

def run_pepsickle_new_loaded2(epitopes, pepsickle_file_loaded, thresh_pcleav=0.5):
    print("Lesen der Proteosomal Cleavage aktuell loaded")

    target_ids = {int(e.database_id) for e in epitopes if hasattr(e, 'database_id') and e.database_id is not None}

    pep_res = pepsickle_file_loaded.copy()
    pep_res['protein_id'] = pep_res['protein_id'].astype(int)
    pep_res = pep_res[pep_res['protein_id'].isin(target_ids)]

    cleav_col = pep_res.columns[1]
    cleavage_dict = {
        pid: group[cleav_col].to_numpy()
        for pid, group in pep_res.groupby('protein_id')
    }

    for e in tqdm(epitopes):
        pid = int(e.database_id)

        if pid not in cleavage_dict:
            for protein_hit in e.protein_hits:
                protein_hit.append(0.0)
            continue

        prob_array = cleavage_dict[pid]

        for protein_hit in e.protein_hits:
            start_pos = protein_hit[1]
            end_pos = protein_hit[2] + 1

            hit_probs = prob_array[start_pos:end_pos]

            if len(hit_probs) > 0:
                pC = float(hit_probs[-1])
                internal_probs = hit_probs[:-1]
                filtered = internal_probs[internal_probs > thresh_pcleav]
                pCon = float(np.prod(1.0 - filtered)) if len(filtered) > 0 else 1.0
            else:
                pC = 0.0
                pCon = 1.0

            protein_hit.append(pCon * pC)

def run_pepsickle_adhoc_modelnumber(epitopes, thresh_pcleav, pep_model):

    for e in tqdm(epitopes):
        e.protein_sequence = e.protein_sequence.rstrip(".1")
        for protein_hit in e.protein_hits:
            p = pepsickle.predict_protein_cleavage_locations(e.protein_sequence, model=pepsickle.initialize_epitope_model())
            if pep_model == 2:
                try:
                    p = pepsickle.predict_protein_cleavage_locations(e.protein_sequence,pepsickle.initialize_digestion_gb_model(), mod_type="in-vitro", proteasome_type="C")
                except:
                    print("gbm model not initializable, using epitope model")
            elif pep_model == 3:
                try:
                    p = pepsickle.predict_protein_cleavage_locations(e.protein_sequence, pepsickle.initialize_digestion_gb_model(),mod_type="in-vitro", proteasome_type="I")
                except:
                    print("gbm model not initializable, using epitope model")
            elif pep_model == 4:
                p = pepsickle.predict_protein_cleavage_locations(e.protein_sequence,pepsickle.initialize_digestion_model(), mod_type="in-vitro-2", proteasome_type="C")
            elif pep_model == 5:
                p = pepsickle.predict_protein_cleavage_locations(e.protein_sequence, pepsickle.initialize_digestion_model(),mod_type="in-vitro-2", proteasome_type="I")

            pep_epitope = pd.DataFrame(p, columns=['position', 'residue', 'cleav_prob', 'cleaved', 'protein_id'])
            probs = list(pep_epitope.iloc[protein_hit[1]:protein_hit[2] + 1, 2])  ## BASTI TODO: evtl andere Spalte nehmen
            pC = probs.pop()  #### Not sure about this calculation? müsste man sich nicht auch die probabilities in den nächsten Positionen / ganzen Protein anschauen?
            pCon = np.prod([1 - prob for prob in probs if prob > thresh_pcleav])
            protein_hit.append(pCon * pC)
