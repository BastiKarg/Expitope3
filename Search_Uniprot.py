"""Protein sequence matching, fixed positions, and indel searches."""

from Epitope import Epitope

from tqdm import tqdm

from operator import itemgetter

from Bio import pairwise2

import re

def missmatches_fix(a, b, fix): # missmatches() returns list with position a[i]!=b[i]. len(list)=mismatches
    for f in fix:
        if a[f-1] != b[f-1]:
            return [1,2,3,4,5,6,7,8,9]

    return [i for i in range(len(a)) if a[i] != b[i]]

def search_in_mut(epi, fullseq):
    epitopes = set()
    protein_hits = []
    for i in range(len(fullseq) - len(epi) + 1):
        comparison = fullseq[i:i + len(epi)]
        distance = len(missmatches_fix(epi, comparison,[]))
        if distance <= 0:
            max_flanking_residues = 9
            n_flanks = fullseq[max(i - max_flanking_residues, 0):i]
            c_flanks = fullseq[i + len(epi):min(len(fullseq), i + len(epi) + max_flanking_residues)]
            protein_hits.append([comparison, i, i + len(epi), distance, n_flanks, c_flanks])
    if protein_hits:
        protein_hits = sorted(protein_hits, key=itemgetter(3))

        uniprot_gene_id, external_gene_name, ensembl_transcript_id, ensembl_gene_id, refSeq_id, links, uniprot_db, expression, flags = "","","","","","","","",""
        database_id = 1
        epitopes.add(Epitope(protein_hits, fullseq, database_id, uniprot_gene_id, external_gene_name, ensembl_transcript_id, ensembl_gene_id, refSeq_id, links, uniprot_db, expression, flags))
        return epitopes

    else:
        return "No hits found!"

def search_for_potential_Epitopes_new_loaded_multicore_indel(epi,k,ouput_dir, fix, uniprot_loaded, database_id,indel):
    proteins = uniprot_loaded
    epitopes = set()
    epitopesX = set()
    print("Epitope werden identifiziert")
    g=4
    score_thresh_small = (len(epi)-1) * 2 - (k-1) * 3 - g
    score_thresh_big = (len(epi) + 1) * 2 - (k-1) * 3 - g

    g_small =  (len(epi) - 1) * 2
    g_big =  (len(epi)) * 2
    score_thresh_small_singlegap = (len(epi) - 1) * 2 - (k - 1) * 3 - g_small
    score_thresh_big_singlegap = (len(epi) + 1) * 2 - (k - 1) * 3 - g_big
    score_thresh_big_singlegap = (len(epi) ) * 2 - (k - 1) * 3 - g_big

    g_small_multi =  (len(epi) - 1)
    g_big_multi =  (len(epi))
    g_small_multi =  2
    g_big_multi =  2

    score_thresh_small_multigap = (len(epi) - 1) * 2 - (k - 1) * 3 - 2 * g_small_multi
    score_thresh_big_multigap = (len(epi) ) * 2 - (k - 1) * 3 - 2 * g_big_multi

    max_res_flank = 5
    for protein in tqdm(proteins): #itarate through protein sequences
        protein_seq = str(protein.seq)
        database_id += 1
        protein_hits = []
        protein_X_hits = []
        for i in range(len(protein_seq)  - len(epi) + 1):
            comparison = protein_seq[i:i + len(epi)]

            distance = len(missmatches_fix(epi,comparison, fix))

            if distance <= k: #distance < min:

                n_flanks = protein_seq[max(i-max_res_flank,0):i]
                c_flanks = protein_seq[i + len(epi):min(len(protein_seq) , i + len(epi) + max_res_flank)]

                if 'X' in comparison:
                    protein_X_hits.append([comparison, i, i + len(epi), distance, n_flanks, c_flanks])
                else:
                    protein_hits.append([comparison, i, i + len(epi), distance, n_flanks, c_flanks])
        if indel == "indel":
            for i in range(len(protein_seq) - len(epi)):  ### sequences that are +1 longer than the epitope
                comparison = protein_seq[i:i + len(epi) + 1]  ## added +1 weil eins länger
                try:
                    for a in pairwise2.align.globalms(comparison, epi, 2, -1, - (g_big), - (g_big)):
                        if a[2] >= score_thresh_big_singlegap:

                            n_flanks = protein_seq[max(i - max_res_flank, 0):i]
                            c_flanks = protein_seq[i + len(epi) + 1:min(len(protein_seq), i + len(
                                epi) + max_res_flank + 1)]  # added +1 weil eins länger

                            if 'X' in comparison:
                                protein_X_hits.append([re.sub("-", "", a[0]), i, i + len(epi) + 1, a[2] / 3, n_flanks, c_flanks])
                            else:
                                protein_hits.append([re.sub("-", "", a[0]), i, i + len(epi) + 1, a[2] / 3, n_flanks, c_flanks])
                except:
                    print("uncommon amino acid in sequence")

            for i in range(len(protein_seq) - len(epi) + 2):  ### sequences that are +1 longer than the epitope
                comparison = protein_seq[i:i + len(epi) - 1]  ## added +1 weil eins länger
                try:
                    for a in pairwise2.align.globalms(comparison, epi, 2, -1, - (g_small), - (g_small)):
                        if a[2] >= score_thresh_small_singlegap:

                            n_flanks = protein_seq[max(i - max_res_flank, 0):i]
                            c_flanks = protein_seq[i + len(epi) - 1:min(len(protein_seq), i + len(
                                epi) + max_res_flank - 1)]  # added +1 weil eins länger

                            if 'X' in comparison:
                                protein_X_hits.append([re.sub("-", "", a[0]), i, i + len(epi) + 1, a[2] / 3, n_flanks, c_flanks])
                            else:
                                protein_hits.append([re.sub("-", "", a[0]), i, i + len(epi) + 1, a[2] / 3 , n_flanks, c_flanks])
                except:
                    print("uncommon amino acid in sequence")

        elif indel == "multigap": # this is now changed ot the previous indel search version, that
            for i in range(len(protein_seq) - len(epi)):  ### sequences that are +1 longer than the epitope
                comparison = protein_seq[i:i + len(epi) + 1]  ## added +1 weil eins länger
                try:
                    for a in pairwise2.align.globalms(comparison, epi, 2, -1, -2, -1):
                        if a[2] >= score_thresh_big:

                            n_flanks = protein_seq[max(i - max_res_flank, 0):i]
                            c_flanks = protein_seq[i + len(epi) + 1:min(len(protein_seq), i + len(
                                epi) + max_res_flank + 1)]  # added +1 weil eins länger

                            if 'X' in comparison:
                                protein_X_hits.append([re.sub("-", "", a[0]), i, i + len(epi) + 1, a[2] / 3, n_flanks, c_flanks])
                            else:
                                protein_hits.append([re.sub("-", "", a[0]), i, i + len(epi) + 1, a[2] / 3, n_flanks, c_flanks])
                except:
                    print("uncommon amino acid in sequence")

            for i in range(len(protein_seq) - len(epi) + 2):  ### sequences that are +1 longer than the epitope
                comparison = protein_seq[i:i + len(epi) - 1]  ## added +1 weil eins länger
                try:
                    for a in pairwise2.align.globalms(comparison, epi, 2, -1, -2, -1):
                        if a[2] >= score_thresh_small:

                            n_flanks = protein_seq[max(i - max_res_flank, 0):i]
                            c_flanks = protein_seq[i + len(epi) - 1:min(len(protein_seq), i + len(
                                epi) + max_res_flank - 1)]  # added +1 weil eins länger

                            if 'X' in comparison:
                                protein_X_hits.append([re.sub("-", "", a[0]), i, i + len(epi) + 1, a[2] / 3, n_flanks, c_flanks])
                            else:
                                protein_hits.append([re.sub("-", "", a[0]), i, i + len(epi) + 1, a[2] / 3 , n_flanks, c_flanks])
                except:
                    print("uncommon amino acid in sequence")

        elif indel == "multigap_not_working":
            for i in range(len(protein_seq) - len(epi)):  ### sequences that are +1 longer than the epitope
                comparison = protein_seq[i:i + len(epi) + 1]  ## added +1 weil eins länger
                try:
                    for a in pairwise2.align.globalms(comparison, epi, 2, -1, - (g_big_multi), 0):
                        if a[2] >= score_thresh_big_multigap:

                            n_flanks = protein_seq[max(i - max_res_flank, 0):i]
                            c_flanks = protein_seq[i + len(epi) + 1:min(len(protein_seq), i + len(
                                epi) + max_res_flank + 1)]  # added +1 weil eins länger

                            if 'X' in comparison:
                                protein_X_hits.append([re.sub("-", "", a[0]), i, i + len(epi) + 1, a[2] / 3, n_flanks, c_flanks])
                            else:
                                protein_hits.append([re.sub("-", "", a[0]), i, i + len(epi) + 1, a[2] / 3, n_flanks, c_flanks])
                except:
                    print("uncommon amino acid in sequence")

            for i in range(len(protein_seq) - len(epi) + 2):  ### sequences that are +1 longer than the epitope
                comparison = protein_seq[i:i + len(epi) - 1]  ## added +1 weil eins länger
                try:
                    for a in pairwise2.align.globalms(comparison, epi, 2, -1, - (g_small_multi), 0):
                        if a[2] >= score_thresh_small_multigap:

                            n_flanks = protein_seq[max(i - max_res_flank, 0):i]
                            c_flanks = protein_seq[i + len(epi) - 1:min(len(protein_seq), i + len(
                                epi) + max_res_flank - 1)]  # added +1 weil eins länger

                            if 'X' in comparison:
                                protein_X_hits.append([re.sub("-", "", a[0]), i, i + len(epi) + 1, a[2] / 3, n_flanks, c_flanks])
                            else:
                                protein_hits.append([re.sub("-", "", a[0]), i, i + len(epi) + 1, a[2] / 3 , n_flanks, c_flanks])
                except:
                    print("uncommon amino acid in sequence")

        if protein_hits:
            protein_hits = sorted(protein_hits, key=itemgetter(3))

            uniprot_gene_id, external_gene_name, ensembl_transcript_id, ensembl_gene_id, refSeq_id, links, uniprot_db, expression, flags = str(protein.description).split(sep='|')[0:9]

            print(str(database_id))

            epitopes.add(Epitope(protein_hits, protein_seq, database_id, uniprot_gene_id, external_gene_name,ensembl_transcript_id, ensembl_gene_id, refSeq_id, links,uniprot_db, expression, flags))
            for i in protein_hits, protein_seq, database_id, uniprot_gene_id, external_gene_name, ensembl_transcript_id, ensembl_gene_id, refSeq_id, links,uniprot_db, expression, flags:
                print(str(i))

        if protein_X_hits:
            protein_X_hits = sorted(protein_X_hits, key=itemgetter(3))

            uniprot_gene_id, external_gene_name, ensembl_transcript_id, ensembl_gene_id, refSeq_id, links = str(protein.name).split(sep='|')[0:6]

            print(str(database_id))
            epitopesX.add(Epitope(protein_X_hits, protein_seq, database_id, uniprot_gene_id, external_gene_name,ensembl_transcript_id, ensembl_gene_id, refSeq_id, links,uniprot_db, expression, flags))
            for i in protein_X_hits, protein_seq, database_id, uniprot_gene_id, external_gene_name, ensembl_transcript_id, ensembl_gene_id, refSeq_id, links,uniprot_db, expression, flags:
                print(str(i))

    from pprint import pprint
    for e in epitopes:
        pprint(vars(e))
    f_string = ""
    for f in fix:
        f_string = f_string + "_" + str(f)

    f = open(ouput_dir + '/' + epi + f_string + '_epitope_search.tsv', 'w')

    f.write("epitope_seq\tdistance\tstart\tend\tuniprot_gene_id\texternal_gene_name\tensembl_transcript_id\tensembl_gene_id\trefSeq_id\tlinks\tuniprot_db\texpression\tflags\n") #

    for e in epitopes:
        for hit in e.protein_hits:

            f.write(str(hit[0]) + "\t" + str(hit[3]) + "\t" + str(hit[1]) + "\t" + str(hit[2]) + "\t" + str(
                e.uniprot_gene_id) + "\t" + str(e.external_gene_name) + "\t" +
                    str(e.ensembl_transcript_id) + "\t" + str(e.ensembl_gene_id) + "\t" + str(e.refSeq_id) + "\t" + str(
                e.links) + "\t" + str(e.uniprot_db) + "\t" + str(e.expression) + "\t" + str(e.flags)+ "\n")
    f.close()
    import pandas
    df = pandas.read_csv(ouput_dir + '/'+ epi + f_string +'_epitope_search.tsv', sep='\t')
    df = df.sort_values(['distance'], ascending = [True])
    df.to_csv(ouput_dir + '/'+ epi + f_string +'_epitope_search.tsv', mode='w', sep='\t',index=False)

    if not epitopesX == set():
        f = open(ouput_dir + '/'+ epi + f_string + '_epitope_search_lowComplexity.tsv', 'w')
        f.write("epitope_seq\tdistance\tstart\tend\tuniprot_gene_id\texternal_gene_name\tensembl_transcript_id\tensembl_gene_id\trefSeq_id\tlinks\n") # braucht man die newline?
        for e in epitopesX:
            for hit in e.protein_hits:
                f.write(str(hit[0]) + "\t" + str(hit[3]) + "\t" + str(hit[1]) + "\t" + str(hit[2]) + "\t" + str(
                    e.uniprot_gene_id) + "\t" + str(e.external_gene_name) + "\t" +
                        str(e.ensembl_transcript_id) + "\t" + str(e.ensembl_gene_id) + "\t" + str(e.refSeq_id) + "\t" + str(
                    e.links) + "\n")
        f.close()

    return epitopes#, epitopesX
