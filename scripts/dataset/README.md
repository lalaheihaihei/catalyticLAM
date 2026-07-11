# Dataset conversion

`convert_vasp_dataset.py` converts VASP AIMD `OUTCAR` files to DeepMD mixed NPY
and extxyz datasets. It retains the original CLAM choices:

- sample every 20 frames by default;
- let dpdata discard electronically unconverged selected frames;
- keep frames with `energy < 0`;
- allow usable frames from an AIMD job that did not reach its requested final step;
- split randomly at frame level.

The random-frame split intentionally allows frames from one trajectory to occur
in train, validation, and test. It measures random-frame interpolation accuracy,
not cross-trajectory generalization.

## Dependencies

```bash
python -m pip install dpdata ase numpy tqdm
```

## Example

```bash
python scripts/dataset/convert_vasp_dataset.py \
  /home/gengzi/workflow-M-N-C/POSCAR \
  --output-dir /data/clam/M-N-C-v2 \
  --stride 20 \
  --seed 2026 \
  --workers 8
```

The validation split uses `seed + 1`. Input paths are sorted before parsing and
parallel results are consumed in that order, so a fixed seed produces a stable
split for unchanged inputs.

## Audit outputs

- `failed_outcars.csv`: OUTCAR files that failed parsing or merging.
- `excluded_frames.csv`: frames removed by the explicit `energy < 0` filter.
  Frames internally rejected by dpdata for electronic nonconvergence are not
  logged separately.
- `conversion_summary.json`: counts and the complete split configuration.

The generated NPY/extxyz datasets are data products and should normally remain
outside Git. Commit the script and configuration as needed.
