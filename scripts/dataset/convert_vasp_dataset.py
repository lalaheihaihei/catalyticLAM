#!/usr/bin/env python3
"""Convert VASP AIMD OUTCAR files to reproducible CLAM datasets.

The split intentionally operates at the sampled-frame level. Frames from one
trajectory may therefore occur in train, validation, and test. This measures
random-frame interpolation accuracy rather than cross-trajectory generalization.
"""

from __future__ import annotations

import argparse
import csv
import json
import shutil
import warnings
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
from tqdm import tqdm


@dataclass
class ParseResult:
    outcar_path: str
    system: Any | None
    excluded_frames: list[dict[str, Any]]
    error_type: str = ""
    error_message: str = ""


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Convert VASP AIMD OUTCAR files to DeepMD NPY and extxyz datasets."
    )
    parser.add_argument(
        "input_roots",
        nargs="+",
        help="Directories containing one-level calculation folders with MD/OUTCAR.",
    )
    parser.add_argument("--output-dir", required=True, help="New output directory.")
    parser.add_argument(
        "--pattern",
        default="*/MD/OUTCAR",
        help="Glob below every input root (default: */MD/OUTCAR).",
    )
    parser.add_argument("--stride", type=int, default=20, help="OUTCAR frame stride.")
    parser.add_argument("--heldout-fraction", type=float, default=0.10)
    parser.add_argument(
        "--validation-fraction-of-heldout",
        type=float,
        default=0.50,
        help="Fraction of held-out frames assigned to validation; the rest become test.",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=2026,
        help="Random seed. The validation split uses seed + 1.",
    )
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument(
        "--formats",
        nargs="+",
        choices=("deepmd-npy", "extxyz"),
        default=("deepmd-npy", "extxyz"),
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Replace an existing output directory.",
    )
    args = parser.parse_args()

    if args.stride < 1:
        parser.error("--stride must be at least 1")
    if not 0.0 <= args.heldout_fraction <= 1.0:
        parser.error("--heldout-fraction must be between 0 and 1")
    if not 0.0 <= args.validation_fraction_of_heldout <= 1.0:
        parser.error("--validation-fraction-of-heldout must be between 0 and 1")
    if args.workers < 1:
        parser.error("--workers must be at least 1")
    return args


def discover_outcars(input_roots: list[str], pattern: str) -> list[str]:
    discovered: set[str] = set()
    for raw_root in input_roots:
        root = Path(raw_root).expanduser().resolve()
        if not root.is_dir():
            raise FileNotFoundError(f"Input root does not exist: {root}")
        for path in root.glob(pattern):
            if path.is_file():
                resolved = str(path.resolve())
                discovered.add(resolved)
    return sorted(discovered)


def process_outcar(task: tuple[str, int]) -> ParseResult:
    outcar_path, stride = task
    try:
        import dpdata as dp

        # dpdata already discards electronically unconverged selected frames.
        # The project intentionally does not maintain a separate report for them.
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            data = dp.LabeledSystem(outcar_path, fmt="vasp/outcar", step=stride)

        energies = np.asarray(data["energies"], dtype=float)
        keep = energies < 0
        excluded = []
        for index in np.flatnonzero(~keep):
            energy = float(energies[index])
            reason = "nonfinite_energy" if not np.isfinite(energy) else "energy_nonnegative"
            excluded.append(
                {
                    "outcar_path": outcar_path,
                    "parsed_frame_index": int(index),
                    "energy": energy,
                    "reason": reason,
                }
            )

        return ParseResult(
            outcar_path=outcar_path,
            system=data[keep],
            excluded_frames=excluded,
        )
    except Exception as error:
        return ParseResult(
            outcar_path=outcar_path,
            system=None,
            excluded_frames=[],
            error_type=type(error).__name__,
            error_message=str(error).replace("\n", " "),
        )


def write_csv(path: Path, fieldnames: list[str], rows: list[dict[str, Any]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def write_extxyz(dataset: Any, path: Path) -> None:
    from ase.io import write

    first_write = True
    for system in dataset:
        structures = system.to_ase_structure()
        if not structures:
            continue
        write(
            path,
            structures,
            format="extxyz",
            append=not first_write,
        )
        first_write = False


def prepare_output_directory(path: Path, overwrite: bool) -> None:
    if path.exists():
        if not overwrite:
            raise FileExistsError(
                f"Output directory already exists: {path}. Use --overwrite to replace it."
            )
        shutil.rmtree(path)
    path.mkdir(parents=True)


def main() -> None:
    args = parse_arguments()

    import dpdata as dp

    output_dir = Path(args.output_dir).expanduser().resolve()
    prepare_output_directory(output_dir, args.overwrite)

    outcars = discover_outcars(args.input_roots, args.pattern)
    if not outcars:
        raise FileNotFoundError(
            f"No OUTCAR files matched pattern '{args.pattern}' in {args.input_roots}"
        )

    tasks = [(path, args.stride) for path in outcars]
    failures: list[dict[str, Any]] = []
    excluded_frames: list[dict[str, Any]] = []
    successful_trajectories = 0
    systems = dp.MultiSystems()

    # executor.map preserves the sorted input order. This is required for the
    # fixed seed to produce the same split across repeated runs.
    with ProcessPoolExecutor(max_workers=args.workers) as executor:
        results = executor.map(process_outcar, tasks, chunksize=1)
        for result in tqdm(results, total=len(tasks), desc="Parsing OUTCAR files"):
            excluded_frames.extend(result.excluded_frames)
            if result.system is None:
                failures.append(
                    {
                        "outcar_path": result.outcar_path,
                        "stage": "parse",
                        "error_type": result.error_type,
                        "error_message": result.error_message,
                    }
                )
                continue
            if len(result.system) == 0:
                failures.append(
                    {
                        "outcar_path": result.outcar_path,
                        "stage": "filter",
                        "error_type": "EmptyAfterFilter",
                        "error_message": "No frames remained after energy < 0 filtering.",
                    }
                )
                continue

            try:
                systems.append(result.system)
                successful_trajectories += 1
            except Exception as error:
                failures.append(
                    {
                        "outcar_path": result.outcar_path,
                        "stage": "merge",
                        "error_type": type(error).__name__,
                        "error_message": str(error).replace("\n", " "),
                    }
                )

    if successful_trajectories == 0:
        write_csv(
            output_dir / "failed_outcars.csv",
            ["outcar_path", "stage", "error_type", "error_message"],
            failures,
        )
        raise RuntimeError("No usable frames were parsed.")

    train_data, heldout_data, _ = systems.train_test_split(
        test_size=args.heldout_fraction,
        seed=args.seed,
    )
    test_data, val_data, _ = heldout_data.train_test_split(
        test_size=args.validation_fraction_of_heldout,
        seed=args.seed + 1,
    )

    write_csv(
        output_dir / "failed_outcars.csv",
        ["outcar_path", "stage", "error_type", "error_message"],
        failures,
    )
    write_csv(
        output_dir / "excluded_frames.csv",
        ["outcar_path", "parsed_frame_index", "energy", "reason"],
        excluded_frames,
    )
    split_data = {
        "train": train_data,
        "test": test_data,
        "val": val_data,
    }
    if "deepmd-npy" in args.formats:
        for name, dataset in split_data.items():
            dataset.to_deepmd_npy_mixed(str(output_dir / name))
    if "extxyz" in args.formats:
        for name, dataset in split_data.items():
            write_extxyz(dataset, output_dir / f"{name}.extxyz")

    summary = {
        "input_roots": [str(Path(root).expanduser().resolve()) for root in args.input_roots],
        "outcar_count": len(outcars),
        "successful_trajectories": successful_trajectories,
        "failed_trajectories": len(failures),
        "energy_filtered_frames": len(excluded_frames),
        "frames": {name: data.get_nframes() for name, data in split_data.items()},
        "stride": args.stride,
        "seed": args.seed,
        "split_strategy": "random_frame",
        "heldout_fraction": args.heldout_fraction,
        "validation_fraction_of_heldout": args.validation_fraction_of_heldout,
        "formats": list(args.formats),
    }
    with (output_dir / "conversion_summary.json").open("w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2, ensure_ascii=False)
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
