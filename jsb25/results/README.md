# Paper Gate-Stream Reproduction Results

`n163_full_stream.csv` was generated with:

```bash
python3 point_addition_resources.py \
    --n 163 --mode both --format csv \
    --output results/n163_full_stream.csv
```

The project retains only the paper scheduling order used by `recursive_karatsuba/inv`:

| mode | width | Toffoli | CNOT | full depth | current depth | TD |
|---|---:|---:|---:|---:|---:|---:|
| FLT-in | 53,134 | 193,354 | 2,417,882 | 5,888 | 5,892 | 46 |
| FLT-out | 57,521 | 48,583 | 650,277 | 1,494 | 1,494 | 12 |

After exact global gate cancellation, FLT-in contains 2,417,882 CNOT gates. No cancellable
gate pair is found in FLT-out under the current non-interference criterion. Every resource
field in `n163_full_stream.csv` is derived from this paper gate stream.

`n163_paper_shor.csv` contains the complete-Shor composition of the same stream:

| mode | Toffoli | CNOT | width | Toffoli-depth | NCT-depth |
|---|---:|---:|---:|---:|---:|---:|
| FLT-in | 63,420,112 | 793,065,296 | 53,134 | 15,088 | 1,931,264 |
| FLT-out | 31,870,448 | 426,582,038 | 15,847,322 | 7,872 | 980,065 |

Both structures have completed full-stream verification for `n=163,233,283`, and FLT-out
has also completed end-to-end verification for `n=571`. All runs passed `q=0/1`, checkpoint
state, clearing, and compute-copy-uncompute checks. The default `n=163` test includes
equivalence of the complete physical state before and after global cancellation. The
`n=571` FLT-in stream was also fully built and scanned, reproducing width 500,450, CNOT
count 30,657,812, and Toffoli count 1,871,402; its raw-stream build reached approximately
305 MiB peak RSS.
