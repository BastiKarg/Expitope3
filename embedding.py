"""ProtT5 embeddings for peptide-level classification."""

from transformers import T5EncoderModel, T5Tokenizer

import torch

import time

torch.cuda.empty_cache()

device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')

print("Using {}".format(device))

def get_T5_model():
    model = T5EncoderModel.from_pretrained("Rostlab/prot_t5_xl_half_uniref50-enc")
    model = model.to(device) # move model to GPU
    model = model.eval() # set model to evaluation model
    tokenizer = T5Tokenizer.from_pretrained('Rostlab/prot_t5_xl_half_uniref50-enc', do_lower_case=False)

    return model, tokenizer

def get_embeddings_local( model, tokenizer, seqs, per_residue, per_protein, sec_struct,
                   max_residues=4000, max_seq_len=1000, max_batch=100 ):
    results = {"residue_embs" : dict(),
               "protein_embs" : dict(),
               "sec_structs" : dict()
               }

    start = time.time()
    batch = list()

    seq = seqs
    seq_len = len(seq)
    pdb_id = "id"
    seq = ' '.join(list(seq))
    batch.append((pdb_id,seq,seq_len))

    n_res_batch = sum([ s_len for  _, _, s_len in batch ]) + seq_len
    if True:
        pdb_ids, seqs, seq_lens = zip(*batch)
        batch = list()

        token_encoding = tokenizer.batch_encode_plus(seqs, add_special_tokens=True, padding="longest")
        input_ids      = torch.tensor(token_encoding['input_ids']).to(device)
        attention_mask = torch.tensor(token_encoding['attention_mask']).to(device)

        try:
            with torch.no_grad():
                embedding_repr = model(input_ids, attention_mask=attention_mask)
        except RuntimeError:
            print("RuntimeError during embedding for {} (L={})".format(pdb_id, seq_len))

        for batch_idx, identifier in enumerate(pdb_ids): # for each protein in the current mini-batch
            s_len = seq_lens[batch_idx]
            emb = embedding_repr.last_hidden_state[batch_idx,:s_len]

            if per_protein: # apply average-pooling to derive per-protein embeddings (1024-d)

                protein_emb = emb.mean(dim=0)
                results["protein_embs"][identifier] = protein_emb.detach().cpu().numpy().squeeze()
    return results

def do_embedding_loaded(epitopes, model, tokenizer):

    per_protein = True
    per_residue = False
    sec_struct = False
    d = {}
    for epitope in epitopes:
        for hit in epitope.protein_hits:
            peptide = hit[0]
            results = get_embeddings_local( model, tokenizer, peptide, per_residue, per_protein, sec_struct)
            hit.extend((results.get('protein_embs').values()))

def do_embedding_loaded_once(epitopes, model, tokenizer):
    per_protein = True
    per_residue = False
    sec_struct = False

    peptides = []
    for epitope in epitopes:
        for hit in epitope.protein_hits:
            peptides.append(hit[0])

    peptides_unique = []
    for x in peptides:
        if x not in peptides_unique:
            peptides_unique.append(x)

    for peptide in peptides_unique:
        results = get_embeddings_local(model, tokenizer, peptide, per_residue, per_protein, sec_struct)
        for epitope in epitopes:
            for hit in epitope.protein_hits:
                if peptide == hit[0]:
                    hit.extend((results.get('protein_embs').values()))
