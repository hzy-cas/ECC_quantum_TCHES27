# Rebuttal reproduction bundle

This directory separates the reproducible code and circuit implementations from the
generated tables and PDFs.

```text
rebuttal/
├── codes_and_circuits/
│   ├── ac_balanced_dw_search/             # C++ AC-arithmetic Balanced estimator
│   ├── algorithm1_inplace_point_add/      # Python in-place AC comparison 
│   ├── binary_ecc_point_add_full_stream/  # Python JSB25 FLT-in/FLT-out reproduction
│   └── windowed_qrom_nct/                 # Coherent Window-QROM artifact
└── results_and_tables/                    # Generated comparison PDFs
```

The full-stream project is self-contained. The three AC-backed projects contain all source
code and locate the repository's canonical, read-only arithmetic matrices in
`mul_line/C++/data` by searching their ancestor directories. 

## Quick checks

Run from the repository root unless a command changes directory explicitly.

```bash
cd rebuttal/codes_and_circuits/binary_ecc_point_add_full_stream
python3 -m unittest discover -s tests -v

cd ../ac_balanced_dw_search
make test
build/ac_balanced_estimator primitive --n 163 --kind multiplication

cd ../algorithm1_inplace_point_add
python3 tests/test_point_addition.py
python3 verification/verify_point_addition.py \
  --sizes 163 233 283 571 --algebra-only

cd ..
python3 -m unittest discover -s windowed_qrom_nct/tests -v
python3 -m windowed_qrom_nct.cli analytic-scan \
  --output-dir /tmp/windowed_qrom_scan
```
