#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import math
import os
import shutil
import sys
from dataclasses import dataclass
from pathlib import Path


SUPPORTED_SUFFIXES = {".sdf", ".csv", ".smiles", ".smi"}


@dataclass
class SdfRecord:
    source_file: Path
    record_index: int
    title: str
    content: str


@dataclass
class SmilesRecord:
    source_file: Path
    record_index: int
    ligand_name: str
    smiles_line: str
    suffix: str


def sanitize_name(value: str, default: str = "ligand") -> str:
    safe = "".join(ch if ch.isalnum() or ch in ("-", "_", ".") else "_" for ch in (value or ""))
    safe = safe.strip("._")
    return safe or default


def _ignore_entry(name: str) -> bool:
    return (
        not name
        or name in {"__MACOSX", ".DS_Store"}
        or name.startswith(".")
        or name.startswith("._")
    )


def _ignore_dir(name: str) -> bool:
    return (
        _ignore_entry(name)
        or name.startswith("Ligands_TMP_SDF_")
        or (name.startswith("Ligands_") and "_PDBQT_" in name)
    )


def discover_input_files(input_path: Path) -> list[Path]:
    input_path = Path(input_path)
    if input_path.is_file():
        if input_path.suffix.lower() not in SUPPORTED_SUFFIXES:
            raise ValueError(f"Unsupported input file: {input_path.name}")
        return [input_path]

    found: list[Path] = []
    for root, dirnames, filenames in os.walk(input_path):
        dirnames[:] = [name for name in dirnames if not _ignore_dir(name)]
        root_path = Path(root)
        for filename in sorted(filenames):
            if _ignore_entry(filename):
                continue
            candidate = root_path / filename
            if candidate.suffix.lower() in SUPPORTED_SUFFIXES:
                found.append(candidate)
    return sorted(found)


def detect_dataset_mode(files: list[Path]) -> str:
    if not files:
        raise ValueError("No supported ligand files were found.")

    suffixes = {path.suffix.lower() for path in files}
    if ".csv" in suffixes:
        if suffixes != {".csv"}:
            raise ValueError("CSV inputs cannot be mixed with SDF or SMILES inputs in one split run.")
        if len(files) != 1:
            raise ValueError("Please provide a single CSV file per split run.")
        return "csv"

    if suffixes <= {".smiles", ".smi"}:
        return "smiles"

    if suffixes == {".sdf"}:
        return "sdf"

    raise ValueError(
        "Mixed ligand input types were detected. Use one run for SDF inputs and a separate run for SMILES inputs."
    )


def count_nonempty_csv_rows(csv_path: Path) -> int:
    count = 0
    with csv_path.open("r", newline="", encoding="utf-8-sig") as handle:
        reader = csv.reader(handle)
        next(reader, None)
        for row in reader:
            if any((cell or "").strip() for cell in row):
                count += 1
    return count


def load_csv_rows(csv_path: Path) -> tuple[list[str], list[list[str]]]:
    with csv_path.open("r", newline="", encoding="utf-8-sig") as handle:
        reader = csv.reader(handle)
        header = next(reader, [])
        rows = [row for row in reader if any((cell or "").strip() for cell in row)]
    return header, rows


def iter_sdf_records(sdf_path: Path):
    buffer: list[str] = []
    record_index = 0
    with sdf_path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            buffer.append(line)
            if line.strip() == "$$$$":
                record_index += 1
                content = "".join(buffer)
                title = buffer[0].strip() if buffer else ""
                yield SdfRecord(
                    source_file=sdf_path,
                    record_index=record_index,
                    title=title,
                    content=content,
                )
                buffer = []

    if buffer and any(part.strip() for part in buffer):
        record_index += 1
        if not buffer[-1].endswith("\n"):
            buffer[-1] = buffer[-1] + "\n"
        buffer.append("$$$$\n")
        title = buffer[0].strip() if buffer else ""
        yield SdfRecord(
            source_file=sdf_path,
            record_index=record_index,
            title=title,
            content="".join(buffer),
        )


def count_sdf_records(sdf_path: Path) -> int:
    return sum(1 for _ in iter_sdf_records(sdf_path))


def parse_smiles_line(line: str, source_file: Path, record_index: int) -> SmilesRecord | None:
    stripped = line.strip()
    if not stripped:
        return None

    parts = stripped.split()
    smiles = parts[0]
    raw_name = " ".join(parts[1:]).strip()
    ligand_name = sanitize_name(raw_name or f"{source_file.stem}_rec{record_index:04d}")
    return SmilesRecord(
        source_file=source_file,
        record_index=record_index,
        ligand_name=ligand_name,
        smiles_line=stripped,
        suffix=source_file.suffix.lower(),
    )


def iter_smiles_records(smiles_path: Path):
    with smiles_path.open("r", encoding="utf-8", errors="replace") as handle:
        record_index = 0
        for line in handle:
            record = parse_smiles_line(line, smiles_path, record_index + 1)
            if record is None:
                continue
            record_index += 1
            yield SmilesRecord(
                source_file=record.source_file,
                record_index=record_index,
                ligand_name=record.ligand_name,
                smiles_line=record.smiles_line,
                suffix=record.suffix,
            )


def count_smiles_records(smiles_path: Path) -> int:
    return sum(1 for _ in iter_smiles_records(smiles_path))


def count_items_for_mode(mode: str, files: list[Path]) -> int:
    if mode == "csv":
        return count_nonempty_csv_rows(files[0])
    if mode == "sdf":
        return sum(count_sdf_records(path) for path in files)
    if mode == "smiles":
        return sum(count_smiles_records(path) for path in files)
    raise ValueError(f"Unknown dataset mode: {mode}")


def plan_groups(total_items: int, group_size: int | None = None, num_groups: int | None = None) -> list[int]:
    if total_items < 1:
        raise ValueError("No ligands were found to split.")
    if group_size is None and num_groups is None:
        raise ValueError("Either group_size or num_groups must be provided.")
    if group_size is not None and num_groups is not None:
        raise ValueError("Provide only one of group_size or num_groups.")

    if group_size is not None:
        if group_size < 1:
            raise ValueError("Group size must be at least 1.")
        num_groups = math.ceil(total_items / group_size)
    else:
        if num_groups is None or num_groups < 1:
            raise ValueError("Number of groups must be at least 1.")
        group_size = math.ceil(total_items / num_groups)

    sizes: list[int] = []
    remaining = total_items
    while remaining > 0:
        current = min(group_size, remaining)
        sizes.append(current)
        remaining -= current
    return sizes


def make_group_dirs(output_root: Path, prefix: str, count: int) -> list[Path]:
    width = max(3, len(str(count)))
    dirs: list[Path] = []
    for idx in range(1, count + 1):
        folder = output_root / f"{prefix}_{idx:0{width}d}"
        folder.mkdir(parents=True, exist_ok=True)
        dirs.append(folder)
    return dirs


def assign_groups(items: list, group_sizes: list[int]) -> list[list]:
    groups: list[list] = []
    cursor = 0
    for size in group_sizes:
        groups.append(items[cursor:cursor + size])
        cursor += size
    return groups


def unique_path(dest_dir: Path, base_name: str, suffix: str, seen: dict[str, int]) -> Path:
    stem = sanitize_name(base_name, default="ligand")
    key = f"{stem}{suffix}".lower()
    count = seen.get(key, 0) + 1
    seen[key] = count
    filename = f"{stem}{suffix}" if count == 1 else f"{stem}_{count:03d}{suffix}"
    return dest_dir / filename


def write_csv_groups(csv_path: Path, output_root: Path, prefix: str, group_sizes: list[int]) -> list[tuple[str, int]]:
    header, rows = load_csv_rows(csv_path)
    groups = assign_groups(rows, group_sizes)
    group_dirs = make_group_dirs(output_root, prefix, len(groups))
    summary: list[tuple[str, int]] = []
    width = max(3, len(str(len(groups))))

    for idx, (group_dir, group_rows) in enumerate(zip(group_dirs, groups), start=1):
        chunk_path = group_dir / f"{sanitize_name(csv_path.stem, 'ligands')}_part_{idx:0{width}d}.csv"
        with chunk_path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle)
            if header:
                writer.writerow(header)
            writer.writerows(group_rows)
        summary.append((group_dir.name, len(group_rows)))

    return summary


def write_sdf_groups(files: list[Path], output_root: Path, prefix: str, group_sizes: list[int]) -> list[tuple[str, int]]:
    records: list[SdfRecord] = []
    for path in files:
        records.extend(iter_sdf_records(path))

    groups = assign_groups(records, group_sizes)
    group_dirs = make_group_dirs(output_root, prefix, len(groups))
    summary: list[tuple[str, int]] = []

    for group_dir, group_records in zip(group_dirs, groups):
        seen: dict[str, int] = {}
        for record in group_records:
            base_name = record.title or f"{record.source_file.stem}_rec{record.record_index:04d}"
            out_path = unique_path(group_dir, base_name, ".sdf", seen)
            out_path.write_text(record.content, encoding="utf-8")
        summary.append((group_dir.name, len(group_records)))

    return summary


def write_smiles_groups(files: list[Path], output_root: Path, prefix: str, group_sizes: list[int]) -> list[tuple[str, int]]:
    records: list[SmilesRecord] = []
    for path in files:
        records.extend(iter_smiles_records(path))

    groups = assign_groups(records, group_sizes)
    group_dirs = make_group_dirs(output_root, prefix, len(groups))
    summary: list[tuple[str, int]] = []

    for group_dir, group_records in zip(group_dirs, groups):
        seen: dict[str, int] = {}
        for record in group_records:
            suffix = ".smiles" if record.suffix == ".smiles" else ".smi"
            out_path = unique_path(group_dir, record.ligand_name, suffix, seen)
            out_path.write_text(record.smiles_line.rstrip() + "\n", encoding="utf-8")
        summary.append((group_dir.name, len(group_records)))

    return summary


def write_summary_file(output_root: Path, dataset_mode: str, input_path: Path, summary: list[tuple[str, int]]):
    summary_path = output_root / "split_summary.csv"
    with summary_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["group_folder", "ligand_count", "dataset_mode", "source_input"])
        for folder_name, count in summary:
            writer.writerow([folder_name, count, dataset_mode, str(input_path)])


def split_ligands(
    input_path: Path,
    output_root: Path,
    group_size: int | None = None,
    num_groups: int | None = None,
    prefix: str = "Ligands_split",
    overwrite: bool = False,
) -> dict:
    files = discover_input_files(input_path)
    dataset_mode = detect_dataset_mode(files)
    total_items = count_items_for_mode(dataset_mode, files)
    group_sizes = plan_groups(total_items, group_size=group_size, num_groups=num_groups)

    output_root = Path(output_root)
    if output_root.exists():
        if not overwrite:
            raise FileExistsError(f"Output folder already exists: {output_root}")
        shutil.rmtree(output_root)
    output_root.mkdir(parents=True, exist_ok=True)

    if dataset_mode == "csv":
        summary = write_csv_groups(files[0], output_root, prefix, group_sizes)
    elif dataset_mode == "sdf":
        summary = write_sdf_groups(files, output_root, prefix, group_sizes)
    else:
        summary = write_smiles_groups(files, output_root, prefix, group_sizes)

    write_summary_file(output_root, dataset_mode, Path(input_path), summary)
    return {
        "mode": dataset_mode,
        "total_items": total_items,
        "group_sizes": group_sizes,
        "output_root": str(output_root),
        "groups": summary,
    }


def prompt_path(prompt: str, default: str | None = None) -> Path:
    suffix = f" [{default}]" if default else ""
    while True:
        raw = input(f"{prompt}{suffix}: ").strip()
        chosen = raw or (default or "")
        if not chosen:
            print("Please enter a path.")
            continue
        path = Path(chosen).expanduser()
        if path.exists():
            return path
        print(f"Path not found: {path}")


def prompt_positive_int(prompt: str, default: int | None = None) -> int:
    suffix = f" [{default}]" if default is not None else ""
    while True:
        raw = input(f"{prompt}{suffix}: ").strip()
        if not raw and default is not None:
            return default
        if raw.isdigit() and int(raw) > 0:
            return int(raw)
        print("Enter a positive whole number.")


def prompt_yes_no(prompt: str, default: bool = True) -> bool:
    default_text = "Y/n" if default else "y/N"
    while True:
        raw = input(f"{prompt} [{default_text}]: ").strip().lower()
        if not raw:
            return default
        if raw in {"y", "yes"}:
            return True
        if raw in {"n", "no"}:
            return False
        print("Please answer y or n.")


def interactive_args() -> argparse.Namespace:
    default_input = "Ligands" if Path("Ligands").exists() else None
    input_path = prompt_path("Ligand input path", default_input)
    files = discover_input_files(input_path)
    dataset_mode = detect_dataset_mode(files)
    total_items = count_items_for_mode(dataset_mode, files)

    print(f"\nDetected input mode: {dataset_mode}")
    print(f"Total ligands found: {total_items}")

    print("\nHow would you like to split them?")
    print(" [1] Set ligands per folder")
    print(" [2] Set total number of folders")
    choice = ""
    while choice not in {"1", "2"}:
        choice = input("Choose 1 or 2: ").strip()

    group_size = None
    num_groups = None
    if choice == "1":
        group_size = prompt_positive_int("Ligands per folder", default=min(100, total_items))
    else:
        num_groups = prompt_positive_int("Number of folders", default=min(10, total_items))

    group_sizes = plan_groups(total_items, group_size=group_size, num_groups=num_groups)
    default_output = f"{sanitize_name(input_path.stem if input_path.is_file() else input_path.name, 'Ligands')}_split"
    output_root = Path(input(f"Output folder [{default_output}]: ").strip() or default_output).expanduser()
    prefix = input("Group folder prefix [Ligands_split]: ").strip() or "Ligands_split"

    print("\nPlanned split:")
    for idx, size in enumerate(group_sizes, start=1):
        print(f" - group {idx}: {size} ligand(s)")

    overwrite = output_root.exists() and prompt_yes_no(f"{output_root} already exists. Overwrite it?", default=False)
    if output_root.exists() and not overwrite:
        print("Aborted before writing files.")
        sys.exit(1)

    if not prompt_yes_no("Proceed with this split?", default=True):
        print("Aborted.")
        sys.exit(1)

    return argparse.Namespace(
        input=input_path,
        output=output_root,
        group_size=group_size,
        num_groups=num_groups,
        prefix=prefix,
        overwrite=overwrite,
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Split ligand inputs into batch-friendly folders. Supports a single CSV, "
            "one or more SDF files, or one or more SMILES/SMI files."
        )
    )
    parser.add_argument("--input", help="Ligand file or folder to split")
    parser.add_argument("--output", help="Destination root folder for split groups")
    parser.add_argument("--group-size", type=int, help="Maximum ligands per output folder")
    parser.add_argument("--num-groups", type=int, help="Number of output folders to create")
    parser.add_argument("--prefix", default="Ligands_split", help="Folder name prefix for each split group")
    parser.add_argument("--overwrite", action="store_true", help="Overwrite the output folder if it already exists")
    parser.add_argument("--yes", action="store_true", help="Run non-interactively using the provided flags")
    args = parser.parse_args()

    if not args.input and not args.yes:
        return interactive_args()

    if not args.input:
        parser.error("--input is required in non-interactive mode")
    if not args.output:
        parser.error("--output is required in non-interactive mode")
    if bool(args.group_size) == bool(args.num_groups):
        parser.error("Provide exactly one of --group-size or --num-groups")
    return args


def print_completion(result: dict):
    print("\nSplit complete.")
    print(f"Mode: {result['mode']}")
    print(f"Total ligands: {result['total_items']}")
    print(f"Output: {result['output_root']}")
    print("Groups:")
    for folder_name, count in result["groups"]:
        print(f" - {folder_name}: {count}")

    print("\nDownstream note:")
    if result["mode"] == "csv":
        print("Use 1_ConformerGeneration.py or 1B_confgen_batch.py in CSV-folder mode for each split folder.")
    elif result["mode"] == "sdf":
        print("Use 1_ConformerGeneration.py or 1B_confgen_batch.py in folder/SDF mode for each split folder.")
    else:
        print("Use 1_ConformerGeneration.py or 1B_confgen_batch.py in folder/SMILES mode for each split folder.")


def main():
    args = parse_args()
    result = split_ligands(
        input_path=Path(args.input).expanduser(),
        output_root=Path(args.output).expanduser(),
        group_size=args.group_size,
        num_groups=args.num_groups,
        prefix=args.prefix,
        overwrite=bool(args.overwrite),
    )
    print_completion(result)


if __name__ == "__main__":
    main()
