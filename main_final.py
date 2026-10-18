"""Expitope3 file-queue worker and prediction pipelines."""

import pandas

import numpy as np

import os

# Must be set before MHCflurry imports TensorFlow.
os.environ.setdefault("TF_USE_LEGACY_KERAS", "1")

import Epitope

import Search_Uniprot

import Get_Expression

import MHCflurry

import TAP_calculation as tab

import Run_Pepsickle

import Logsums_and_ICR_calculation

import Output_File

import embedding

import mlp_predict

from multiprocessing import SimpleQueue

from multiprocessing.pool import Pool

import time

from mhcflurry import class1_presentation_predictor as mhc

def init_worker(shared_queue):
    global queue
    queue = shared_queue

def task_process(identifier):
    global queue, fix, epi, k, output_dir, proteins_split, pep_model, HLA_allel, predictor, search_index, indel_bool
    alpha = 0.2
    es = Search_Uniprot.search_for_potential_Epitopes_new_loaded_multicore_indel(epi, k, output_dir, fix, proteins_split[identifier], search_index[identifier],indel_bool)

    queue.put(es)

def list_full_paths(directory):
    return [os.path.join(directory, file) for file in os.listdir(directory)]

def get_fixed_pos(epi):
    fix=[]
    if epi.isupper():
        return epi, fix
    elif epi.islower():
        return epi.upper(), fix
    else:
        fix = [idx + 1 for idx in range(len(epi)) if epi[idx].isupper()]
        return epi.upper(), fix

def run_mut_loaded(epi, fullseq, id, output_dir, HLA_allele, pep_model, predictor, model, tokenizer, mlp_model):
    alpha = 0.2
    thresh_pcleav=0.5
    pcleav_predictor="normal"

    epitopes = Search_Uniprot.search_in_mut(epi, fullseq)
    Run_Pepsickle.run_pepsickle_adhoc_modelnumber(epitopes, thresh_pcleav, pep_model)  ## TODO mit Args verknüpfen

    MHCflurry.predict_new_loaded(epitopes, HLA_allele,predictor)
    tab.calculate(epitopes, alpha)  # calculates and sets tap_score for epitopes

    Epitope.calculate_q_score_new(epitopes)
    f_string = ""
    hla_string = '_' + HLA_allele[0].replace('*', '_').replace(':', '_').split('-')[1]

    Get_Expression.write_scores_mut(epitopes, output_dir, epi,hla_string)

    print("doing some casual embedding")
    embedding.do_embedding_loaded(epitopes, model, tokenizer)

    mlp_predict.mlp_prediction_multi_in_one_loaded_bench_create_trainingset(epitopes, output_dir, epi, f_string, hla_string, mlp_model)

def run_waiting_loaded_multicore():  # epi, k, id, output_dir, HLA_allele=[], pcleav_predictor="e", thresh_pcleav=0.5, alpha=0.2, threshQ=0.02,
    global fix, epi, k, output_dir, pep_model, HLA_allel, predictor, indel_bool

    from Bio import SwissProt
    from Bio import SeqIO
    indel_bool = False

    pepsickle_filese = ['database5/pepsickle_default.tsv', 'database5/pepsickle_gbm_C.tsv',
                        'database5/pepsickle_gbm_I.tsv', 'database5/pepsickle_nn_C.tsv',
                        'database5/pepsickle_nn_I.tsv', ]
    pep_modles = [pandas.read_table(pepsickle_filese[0]),pandas.read_table(pepsickle_filese[1]),pandas.read_table(pepsickle_filese[2]),pandas.read_table(pepsickle_filese[3]),pandas.read_table(pepsickle_filese[4]),]

    from mhcflurry import class1_presentation_predictor as mhc
    predictor = mhc.Class1PresentationPredictor.load()

    import torch
    device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')
    if device.type != 'cpu':
        torch.cuda.set_per_process_memory_fraction(0.45, 0)

    print("Using {}".format(device))
    print("we are in the function run_waiting_loaded_multicore")
    model, tokenizer = embedding.get_T5_model()
    proteinDB="database5/database_october.fa"
    proteins = list(SeqIO.parse(open(proteinDB),'fasta'))

    mlp_models_loaded = [mlp_predict.MLP.load_from_checkpoint("./best_models/best_val_f1.2_9547_113.ckpt", weights_only=True), mlp_predict.MLP.load_from_checkpoint("./best_models/atlas_best_val_f1.2_08846_47.ckpt", weights_only=True), mlp_predict.MLP.load_from_checkpoint("./best_models/best_val_f1.2_9547_113.ckpt", weights_only=True)]
    mlp_models_loaded[0].eval()

    mlp_models_loaded[1].eval()
    print("loaded second model")

    mlp_models_loaded[2].eval()
    print("loaded all models")

    id = "test"
    import time
    file_path = "running"
    os.makedirs(file_path, exist_ok=True)
    os.makedirs("recent_finish", exist_ok=True)
    print("Waiting for jobs")
    while True:
        while len(list_full_paths(file_path)) == 0:
            time.sleep(1)
        else:
            job_todo = max(list_full_paths(file_path), key=os.path.getctime)  ### this was changed to work from latest to oldest
            job = list(pandas.read_csv(job_todo, sep=' ').columns)
            if (len(job) == 7):
                epi, hla, mismatches, outdir, pep_modle, mlp_num, kmer = job
                HLA_allel = hla
                mismatches = int(mismatches)
                mlp_num = int(float(mlp_num))

                pep_modle = int(float(pep_modle))
                k = mismatches
                pep_model = pep_modles[pep_modle - 1]
                outdir = "./results/" + outdir
                if kmer == "indel":
                    indel_bool = "indel"
                    kmer = "normal_seq"
                elif kmer == "multigap":
                    indel_bool = "multigap"
                    kmer = "normal_seq"
                else:
                    indel_bool = False
                if not os.path.exists(outdir):
                    os.makedirs(outdir)
                if kmer == "normal_seq":
                    run_new_app_v3_loaded_multicore(epi, mismatches, id, outdir, [hla], pep_modles[pep_modle - 1], predictor,
                                          model, tokenizer, mlp_models_loaded[mlp_num - 1], proteins)
                else:
                    run_new_app_v3_loaded_multicore_kmer(epi, mismatches, id, outdir, [hla], pep_modles[pep_modle - 1], predictor,
                                          model, tokenizer, mlp_models_loaded[mlp_num - 1], proteins, int(kmer))
                os.replace(job_todo, "./recent_finish/last_finish")

            elif (len(job) == 8):
                epi, hla, mismatches, outdir, pep_modle, mlp_num, kmer, custom_seq = job # custom_seq, kmer
                HLA_allel = hla
                mismatches = int(mismatches)
                mlp_num = int(float(mlp_num))

                pep_modle = int(float(pep_modle))
                k = mismatches
                pep_model = pep_modles[pep_modle - 1]
                outdir = "./results/" + outdir

                if kmer == "indel":
                    indel_bool = "indel"
                    kmer = "normal_seq"
                elif kmer == "multigap":
                    indel_bool = "multigap"
                    kmer = "normal_seq"

                custom_seq = custom_seq.upper()
                if not os.path.exists(outdir):
                    os.makedirs(outdir)
                if kmer == "normal_seq":
                    run_mut_loaded(epi, custom_seq, id, outdir, [hla], pep_modle, predictor, model, tokenizer,
                                   mlp_models_loaded[mlp_num - 1])
                else:
                    for seq in MHCflurry.get_all_kmers(epi, int(kmer)):
                        print(seq)
                        run_mut_loaded(seq, custom_seq, id, outdir, [hla], pep_modle, predictor, model, tokenizer,
                                       mlp_models_loaded[mlp_num - 1])
                os.replace(job_todo, "./recent_finish/last_finish")

def run_new_app_v3_loaded_multicore(epi_x, k_x, id, output_dir_x, HLA_allele, pep_model_x, predictor_x, model, tokenizer, mlp_model, proteins):

    global queue, fix, epi, k, output_dir, proteins_split, pep_model, HLA_allel, predictor
    global search_index
    global indel_bool

    search_index = [0, 39432, 78864, 118296, 157728, 197160, 236592, 276024, 315456, 354888, 394320, 433752, 473182, 512612, 552042, 591472]

    epi = epi_x
    k = k_x
    output_dir = output_dir_x
    pep_model = pep_model_x
    predictor = predictor_x

    expression_databases = "Datenbanken_Ordner_v5"

    threshQ = 0.02
    alpha = 0.2
    thresh_pcleav = 0
    pcleav_predictor = "normal"
    epi, fix = get_fixed_pos(epi)

    number_of_task = 15
    print("Starting " + str(number_of_task) + " tasks")
    proteins_split = np.array_split(np.asarray(proteins,dtype="object"), number_of_task)
    epitopes = set()
    shared_queue = SimpleQueue()
    print("Indel is " + str(indel_bool))
    print("starting with pool")

    with Pool(initializer=init_worker, initargs=(shared_queue,)) as pool:
        _ = pool.map_async(task_process, range(number_of_task))
        for _ in range(number_of_task):
            epitopes.update(shared_queue.get())

    print(len(epitopes))

    if (len(epitopes) == 0):
        print(f"[SKIP] Keine Epitope gefunden für {epi} (HLA: ...). Überspringe Eintrag...")
        return

    Run_Pepsickle.run_pepsickle_new_loaded2(epitopes, pep_model, thresh_pcleav=0.5)

    MHCflurry.predict_new_loaded_once2(epitopes, HLA_allele, predictor)

    tab.calculate(epitopes, alpha)  # calculates and sets tap_score for epitopes

    Epitope.calculate_q_score_new2(epitopes)
    f_string = ""
    for f in fix:
        f_string = f_string + "_" + str(f)
    if indel_bool:
        f_string = f_string + "_INDEL"
    hla_string = '_' + HLA_allele[0].replace('*', '_').replace(':', '_').split('-')[1]

    print("doing some casual embedding")
    embedding.do_embedding_loaded_once(epitopes, model, tokenizer)

    mlp_predict.mlp_prediction_reloaded(epitopes, output_dir, epi, f_string, hla_string, mlp_model)

    expression_data = Get_Expression.search_expression_databases_new_v5_collapsed(epitopes, expression_databases, output_dir, epi, f_string, hla_string)

    print("expression fertig")

    log_sums, tissue_sums = Logsums_and_ICR_calculation.calculate_log_sums_and_icr_v2(expression_data, threshQ, k, output_dir)
    Output_File.write_Expression_to_file(expression_data, epi=epi + f_string, id2=id, output_dir=output_dir)  ### added + f_string here to print the positions of fixed amino acids. Quick and dirty
    Output_File.write_LogSums_to_file(log_sums, epi + f_string, id, output_dir)
    Output_File.write_TissueSums_to_file(tissue_sums, epi + f_string, id, output_dir)
    print("Wrote tables")

def task_process_kmer(identifier):
    global queue, fix_k, epi_k, k, output_dir, proteins_split, pep_model, HLA_allel, predictor, search_index, indel_bool#, thread_unifier
    es = Search_Uniprot.search_for_potential_Epitopes_new_loaded_multicore_indel(epi_k, k, output_dir, fix_k, proteins_split[identifier], search_index[identifier],indel_bool)
    queue.put(es)

def run_new_app_v3_loaded_multicore_kmer(epi_x, k_x, id, output_dir_x, HLA_allele, pep_model_x, predictor_x, model, tokenizer, mlp_model, proteins, kmer):
    global queue, fix, epi, k, output_dir, proteins_split, pep_model, HLA_allel, predictor
    global search_index
    global epi_k, fix_k
    global thread_unifier

    epi = epi_x
    k = k_x
    output_dir = output_dir_x
    pep_model = pep_model_x
    predictor = predictor_x

    expression_databases = "Datenbanken_Ordner_v5"

    threshQ = 0.02
    alpha = 0.2
    thresh_pcleav = 0
    pcleav_predictor = "normal"
    epi, fix = get_fixed_pos(epi)

    number_of_task = 24
    thread_unifier = []

    print("Starting " + str(number_of_task) + " tasks")
    proteins_split = np.array_split(np.asarray(proteins, dtype="object"), number_of_task)

    search_index = [0]
    for i in range(len(proteins_split)):
        search_index.append(search_index[-1] + len(proteins_split[i]))
    print(search_index)
    epitopes = set()
    shared_queue = SimpleQueue()
    print("starting with pool")

    max_retries = 5
    print("starting loop with chat")

    from multiprocessing import Queue
    shared_queue = Queue()

    for seq in MHCflurry.get_all_kmers(epi, int(kmer)):
        epi_k, fix_k = get_fixed_pos(seq)
        with Pool(initializer=init_worker, initargs=(shared_queue,)) as pool:
            success = False

            for attempt in range(max_retries):
                result = pool.map_async(task_process_kmer, range(number_of_task))

                result.wait()

                task_results = []
                try:
                    for _ in range(number_of_task):
                        task_results.append(shared_queue.get(timeout=100))  # Timeout to avoid indefinite waiting

                    success = True
                    break
                except Exception as e:
                    print(f"Attempt {attempt + 1}/{max_retries} failed: {e}", flush=True)
                    time.sleep(10) # bei Problemen zu 10 wechseln

                    while not shared_queue.empty():
                        shared_queue.get()

            if success:
                for result in task_results:
                    for r in result:
                        epitopes.add(r)
            else:
                print("Failed to process sequence after maximum retries.", flush=True)
                exit(0)

            print(f"Current number of epitopes: {len(epitopes)}")
            pool.close()  # Prevent new tasks
            pool.join()   # Wait for all worker processes to finish
            pool.terminate()

    print(f"Final number of epitopes: {len(epitopes)}")

    Run_Pepsickle.run_pepsickle_new_loaded(epitopes, pep_model, thresh_pcleav=0.5)

    MHCflurry.predict_new_loaded_once(epitopes, HLA_allele, predictor)

    tab.calculate(epitopes, alpha)  # calculates and sets tap_score for epitopes

    Epitope.calculate_q_score_new(epitopes)
    f_string = ""
    for f in fix:
        f_string = f_string + "_" + str(f)
    hla_string = '_' + HLA_allele[0].replace('*', '_').replace(':', '_').split('-')[1]

    print("doing some casual embedding")
    epi = epi + str(kmer)
    embedding.do_embedding_loaded_once(epitopes, model, tokenizer)
    mlp_predict.mlp_prediction_loaded(epitopes, output_dir, epi, f_string, hla_string, mlp_model)

    expression_data = Get_Expression.search_expression_databases_new_v4_collapsed(epitopes, expression_databases,
                                                                                  output_dir, epi, f_string, hla_string)
    print("expression fertig")

    log_sums, tissue_sums = Logsums_and_ICR_calculation.calculate_log_sums_and_icr_v2(expression_data, threshQ, k,
                                                                                      output_dir)
    Output_File.write_Expression_to_file(expression_data, epi=epi + f_string, id2=id,
                                         output_dir=output_dir)  ### added + f_string here to print the positions of fixed amino acids. Quick and dirty
    Output_File.write_LogSums_to_file(log_sums, epi + f_string, id, output_dir)
    Output_File.write_TissueSums_to_file(tissue_sums, epi + f_string, id, output_dir)
    print("Wrote tables")

if __name__ == "__main__":
    run_waiting_loaded_multicore()
