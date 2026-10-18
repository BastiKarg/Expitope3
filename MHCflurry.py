"""MHCflurry affinity and presentation prediction adapters."""

def predict_new_loaded(epitopes,allele, predictor):
    if allele == []:
        allele = ['HLA-A0201']

    for epitope in epitopes:
        for hit in epitope.protein_hits:
            peptide = [hit[0]]
            results = predictor.predict(peptides=peptide, alleles=allele)
            print(results.to_string())

            score = results.iloc[0,3]
            processing_score = results.iloc[0, 5]
            presentation_score = results.iloc[0, 6]
            percentile = results.iloc[0, 7]

            hit.extend((score, processing_score,presentation_score, percentile))
            print(results.to_string())

            if str(hit[4]) == 'nan':
                results = predictor.predict(peptides=[hit[0]], alleles=allele,
                                            n_flanks=[''], c_flanks=[hit[5]],
                                            verbose=0)  # HLA)
            elif str(hit[5]) == 'nan':
                results = predictor.predict(peptides=[hit[0]], alleles=allele,
                                            n_flanks=[hit[4]], c_flanks=[''],
                                            verbose=0)  # HLA)
            else:
                results = predictor.predict(peptides=[hit[0]], alleles=allele,
                                            n_flanks=[hit[4]], c_flanks=[hit[5]],
                                            verbose=0)  # HLA)
            print(results.to_string())

            score = results.iloc[0,5]
            processing_score = results.iloc[0, 7]
            presentation_score = results.iloc[0, 8]
            percentile = results.iloc[0, 9]
            hit.extend((score, processing_score,presentation_score, percentile))

def get_all_kmers(protein_seq, k):
    res = [protein_seq[i: j] for i in range(len(protein_seq)) for j in range(i + 1, len(protein_seq) + 1) if
           len(protein_seq[i:j]) == k]
    return res

def predict_new_loaded_once(epitopes,allele, predictor):
    if allele == []:
        allele = ['HLA-A0201']

    peptides = []
    n_flanks = []
    c_flanks = []

    for epitope in epitopes:
        for hit in epitope.protein_hits:
            peptides.append(hit[0])
            n_flanks.append(hit[4])
            c_flanks.append(hit[5])

    c_flanks = [s.replace('nan', '') for s in c_flanks]
    n_flanks = [s.replace('nan', '') for s in n_flanks]
    print(peptides)
    results = predictor.predict(peptides=peptides, alleles=allele)
    results_flank = predictor.predict(peptides=peptides, alleles=allele, c_flanks=c_flanks, n_flanks=n_flanks)
    print(results.to_string())
    print(results_flank.to_string())
    counter = 0
    for epitope in epitopes:
        for hit in epitope.protein_hits:
            score = results.iloc[counter, 3]
            processing_score = results.iloc[counter, 5]
            presentation_score = results.iloc[counter, 6]
            percentile = results.iloc[counter, 7]
            hit.extend((score, processing_score, presentation_score, percentile))

            score = results_flank.iloc[counter,5]
            processing_score = results_flank.iloc[counter, 7]
            presentation_score = results_flank.iloc[counter, 8]
            percentile = results_flank.iloc[counter, 9]
            hit.extend((score, processing_score,presentation_score, percentile))
            counter += 1

def predict_new_loaded_once2(epitopes, allele, predictor):
    if not allele:
        allele = ['HLA-A0201']

    peptides = []
    n_flanks = []
    c_flanks = []
    all_hits = []

    for epitope in epitopes:
        for hit in epitope.protein_hits:
            peptides.append(hit[0])

            n_val = str(hit[4]) if hit[4] is not None else ""
            c_val = str(hit[5]) if hit[5] is not None else ""

            n_flanks.append("" if n_val in ("nan", "None") else n_val)
            c_flanks.append("" if c_val in ("nan", "None") else c_val)

            all_hits.append(hit)

    if not peptides:
        return

    results = predictor.predict(peptides=peptides, alleles=allele)
    results_flank = predictor.predict(
        peptides=peptides,
        alleles=allele,
        c_flanks=c_flanks,
        n_flanks=n_flanks
    )

    res_matrix = results.iloc[:, [3, 5, 6, 7]].to_numpy()

    flank_matrix = results_flank.iloc[:, [5, 7, 8, 9]].to_numpy()

    for hit, res_row, flank_row in zip(all_hits, res_matrix, flank_matrix):
        hit.extend(res_row.tolist())
        hit.extend(flank_row.tolist())
