"""Protein-hit container and Q-score calculations."""

import math

class Epitope:
    def __init__(self, protein_hits, protein_sequence, database_id, uniprot_gene_id, external_gene_name,
                 ensembl_transcript_id,ensembl_gene_id, refSeq_id, links,uniprot_db, expression, flags):
        self.protein_sequence = protein_sequence
        self.database_id = database_id
        self.protein_hits = protein_hits  # [protein_sequence, start_pos, end_pos, distance]
        self.uniprot_gene_id = uniprot_gene_id
        self.external_gene_name = external_gene_name
        self.ensembl_gene_id = ensembl_gene_id
        self.ensembl_transcript_id = ensembl_transcript_id
        self.refSeq_id = refSeq_id
        self.links = links
        self.pcleav_score = None
        self.tap_score = None
        self.tap_score10 = None
        self.mhc_score = None
        self.mhc_processing_score = None
        self.q_score = None
        self.tap_new = None
        self.mhc_affinity = None
        self.q_score_corrected = None
        self.q_score_new = None

        self.uniprot_db = uniprot_db
        self.expression = expression
        self.flags = flags

        self.hla = None

def calculate_q_score_new(epitopes):
    print("start of calculate q score")
    for epitope in epitopes:
        for hit in epitope.protein_hits:
            print(hit)
            print(str(len(hit)))
            q_score = hit[6]/(hit[16] * hit[11])
            hit.append(q_score)

            tap_log = math.log10(hit[15] + 8.75)
            q_score = hit[6]/(tap_log * hit[11])
            hit.append(q_score)

            q_score = hit[6] * hit[17] /  hit[11]
            hit.append(q_score)

def calculate_q_score_new2(epitopes):
    print(f"start of calculate q score for {len(epitopes)} epitopes")

    for epitope in epitopes:
        for hit in epitope.protein_hits:

            mhc_val = hit[11] if hit[11] != 0 else 1e-6
            tap10_val = hit[16] if hit[16] != 0 else 1e-6

            q_score_std = hit[6] / (tap10_val * mhc_val)
            hit.append(q_score_std)

            tap_log = math.log10(hit[15] + 8.75) if (hit[15] + 8.75) > 0 else 1e-6
            q_score_log = hit[6] / (tap_log * mhc_val)
            hit.append(q_score_log)

            q_score_2004 = (hit[6] * hit[17]) / mhc_val
            hit.append(q_score_2004)

    print("finished calculate q score")
