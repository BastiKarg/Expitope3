"""Compare actual ProtT5 peptide embeddings between the validated CPU and a GPU environment."""
import argparse
import json
from pathlib import Path
import sys
import time
import numpy as np
import torch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import embedding


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--device', choices=['cpu', 'cuda'], required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--reference', type=Path)
    args = parser.parse_args()
    torch.set_num_threads(4)
    if args.device == 'cuda':
        if not torch.cuda.is_available():
            raise RuntimeError('CUDA is unavailable')
        torch.backends.cuda.matmul.allow_tf32 = False
        torch.backends.cudnn.allow_tf32 = False
    embedding.device = torch.device(args.device)
    model, tokenizer = embedding.get_T5_model()
    peptides = ['ALSYIFCLV', 'SYFPEITHI', 'LLGPLLVLQA', 'SLLMWITQC', 'GILGFVFTL', 'MKTIIALSYIFCLVF']

    def calculate():
        return np.stack([embedding.get_embeddings_local(model, tokenizer, peptide, False, True, False)
                         ['protein_embs']['id'] for peptide in peptides])

    calculate()  # Warm-up is excluded from timings.
    timings = []
    for _ in range(3):
        if args.device == 'cuda':
            torch.cuda.synchronize()
        start = time.perf_counter()
        values = calculate()
        if args.device == 'cuda':
            torch.cuda.synchronize()
        timings.append(time.perf_counter() - start)
    assert values.shape == (len(peptides), 1024) and np.isfinite(values).all()
    report = dict(torch=torch.__version__, device=args.device, peptides=len(peptides),
                  median_seconds=float(np.median(timings)), timings=timings)
    if args.device == 'cuda':
        report.update(gpu=torch.cuda.get_device_name(), peak_gpu_bytes=torch.cuda.max_memory_allocated())
    if args.reference:
        reference = np.load(args.reference)
        assert np.array_equal(reference['peptides'], peptides)
        report['max_absolute_difference'] = float(np.max(np.abs(values-reference['values'])))
        np.testing.assert_allclose(values, reference['values'], rtol=1e-3, atol=1e-4)
        report['comparison_passed'] = True
    np.savez(args.output, peptides=peptides, values=values)
    args.output.with_suffix('.json').write_text(json.dumps(report, indent=2))
    print(json.dumps(report, indent=2), flush=True)


if __name__ == '__main__':
    main()
