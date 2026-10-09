#!/usr/bin/env python3
"""Portable Window-QROM: retained/min-TD or Balanced/min-TDW estimation."""
from __future__ import annotations

import argparse
import concurrent.futures as futures
import csv
import hashlib
import json
import os
from pathlib import Path
import platform
import shutil
import signal
import subprocess
import sys
import threading
import tempfile
import time

import model
from support.basis import poly_to_l
from support.curves import CURVES
from support.fixtures import make_fixture

ROOT = Path(__file__).resolve().parent
EXECUTABLE = ROOT / "build" / "retained_estimator"
CHILDREN = set()
CHILDREN_LOCK = threading.Lock()
STOP = threading.Event()
# Retained gate-generation sources/data are unchanged. The added Balanced
# entry point/header/build target changes the package-wide hash, not retained
# arithmetic. Only these audited releases can reuse old RETAINED records;
# Balanced records always require an exact fingerprint match. Build stamps
# are never covered by this exception (the old binary lacks the new option).
LEGACY_LAYER_IDENTITIES = {
    "12617d405bb508009e3d386491d294ee1a44ee1a11d3c5bbb69082598d206b9d":
    "3f1b57558749c2c40466cfb013c6df2938859f11007eaa17c5a17b802fddeace",
    "79ad7909db7c967357e1bb8069fb96198863f891cf2256daddbf332d9229c0c9":
    "3f1b57558749c2c40466cfb013c6df2938859f11007eaa17c5a17b802fddeace"
}


def atomic_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + f".{os.getpid()}.{threading.get_ident()}.tmp")
    try:
        with temporary.open("w", encoding="utf-8") as handle:
            json.dump(value, handle, indent=2, sort_keys=True)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def layer_fingerprint():
    paths = [ROOT / "CMakeLists.txt", ROOT / "data_manifest.json"]
    paths += sorted((ROOT / "backend").rglob("*.cpp"))
    paths += sorted((ROOT / "backend").rglob("*.h"))
    paths += sorted((ROOT / "support").glob("*.py"))
    paths += sorted((ROOT.parent / "circuit_backend" / "cpp").glob("GateManager.*"))
    digest = hashlib.sha256()
    for path in paths:
        digest.update(os.path.relpath(path, ROOT).encode())
        digest.update(sha256(path).encode())
    return digest.hexdigest()


def fingerprint(objective="td"):
    """Network identity is separate from the unchanged point-add circuits."""
    digest = hashlib.sha256(layer_fingerprint().encode())
    for name in ("model.py", "run.py"):
        digest.update(sha256(ROOT / name).encode())
    digest.update(model.architecture(objective).encode())
    return digest.hexdigest()


def compatible_layer_identity(stored, current):
    return stored == current or LEGACY_LAYER_IDENTITIES.get(stored) == current


def manifest(create=False):
    destination = ROOT / "data_manifest.json"
    files = sorted(path for path in (ROOT / "data").rglob("*") if path.is_file())
    if any(path.is_symlink() for path in files):
        raise ValueError("data must be real bundled files, not symlinks")
    actual = {str(path.relative_to(ROOT)): sha256(path) for path in files}
    if not actual:
        raise ValueError("no bundled data")
    if create:
        if destination.exists():
            raise ValueError("manifest already exists; use verify-data to audit it")
        atomic_json(destination, actual)
    else:
        expected = json.loads(destination.read_text())
        if actual != expected:
            bad = sorted(key for key in actual.keys() | expected.keys() if actual.get(key) != expected.get(key))
            raise ValueError(f"bundled data checksum mismatch: {bad[:10]}")
    print(f"Verified {len(actual)} bundled data files.", flush=True)


def check_build():
    info = json.loads((ROOT / "build" / "build_stamp.json").read_text())
    if info.get("fingerprint") != layer_fingerprint() or info["executable_sha256"] != sha256(EXECUTABLE):
        raise ValueError("source/data manifest/build changed; run python3 run.py build again")
    if "platform" in info and info["platform"] != [platform.system(), platform.machine()]:
        raise ValueError("the executable was built for another OS/CPU architecture")
    # Old build stamps did not record the platform. Probe once, before workers,
    # so a copied Mach-O/ARM binary never produces one failure per layer.
    subprocess.run([str(EXECUTABLE), "--help"], check=True, timeout=15,
                   stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)


def build_project(jobs):
    if not shutil.which("cmake") or not shutil.which("ctest"):
        raise RuntimeError("CMake/CTest missing. On Ubuntu: apt-get install cmake build-essential")
    directory = ROOT / "build"
    foreign = False
    cache = directory / "CMakeCache.txt"
    if cache.is_file():
        for line in cache.read_text().splitlines():
            for key, expected in (("CMAKE_HOME_DIRECTORY:INTERNAL=", ROOT),
                                  ("CMAKE_CACHEFILE_DIR:INTERNAL=", directory)):
                if line.startswith(key) and Path(line[len(key):]).resolve() != expected.resolve():
                    foreign = True
    if EXECUTABLE.is_file():
        try:
            subprocess.run([str(EXECUTABLE), "--help"], check=True, timeout=15,
                           stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
        except (OSError, subprocess.SubprocessError):
            foreign = True
    if foreign and directory.exists():
        backup = Path(tempfile.mkdtemp(prefix="build.backup.", dir=ROOT)) / "build"
        directory.rename(backup)
        print(f"Preserved incompatible build at {backup}; compiling for this machine.", flush=True)
    subprocess.run(["cmake", "-S", str(ROOT), "-B", str(directory), "-DCMAKE_BUILD_TYPE=Release"], check=True)
    subprocess.run(["cmake", "--build", str(directory), "-j", str(jobs)], check=True)
    # CTest --test-dir is newer than Ubuntu 20.04's CMake 3.16.
    subprocess.run(["ctest", "--output-on-failure"], cwd=directory, check=True)
    atomic_json(directory / "build_stamp.json", dict(fingerprint=layer_fingerprint(),
        executable_sha256=sha256(EXECUTABLE), platform=[platform.system(), platform.machine()]))


def ensure_build(jobs):
    try:
        check_build()
    except (OSError, ValueError, KeyError, subprocess.SubprocessError) as error:
        print(f"Native build required ({error}); rebuilding automatically.", flush=True)
        build_project(min(jobs, 8))
        check_build()


def cache_path(output, task):
    n, kind, k = task
    return output / "layers" / f"n{n}_{kind}_k{k}.json"


def read_record(path, task, identity, objective="td"):
    record = json.loads(path.read_text())
    if not (record.get("fingerprint") == identity or
            (objective == "td" and compatible_layer_identity(record.get("fingerprint"), identity))):
        raise ValueError(f"stale/foreign cache: {path}; use a new output directory")
    model.validate_record(record, *task, objective)
    return record


def estimated_gb(task, objective="td"):
    n, _, k = task
    # Deliberately coarse planning estimates, NOT measured upper bounds.
    # Balanced emits the arithmetic inverse as well. Keep a separate,
    # conservative estimate rather than reusing retained RSS measurements.
    factor = 1 if objective == "td" else 2
    return 1.0 + factor * k * {163: 0.08, 233: 0.20, 283: 0.35, 571: 1.0}[n]


def task_plan(args):
    tasks = set()
    candidates = []
    for n in args.n:
        selected = model.formula_minimizers(n, args.s_min, args.s_max, getattr(args, "objective", "td"))
        candidates.extend(selected)
        if args.scope == "smoke":
            tasks.update((n, kind, k) for kind in ("tail", "reduction") for k in (1, 2, 3))
        elif args.scope == "winners":
            for row in selected:
                tasks.update(model.required_layers(n, row["s"], row["p"]))
        else:
            for s, p in model.grid(n, args.s_min, args.s_max):
                tasks.update(model.required_layers(n, s, p))
    return sorted(tasks), candidates


def check_output_architecture(args):
    """Fail before any writes if this directory belongs to the other mode."""
    summary = args.output / "summary.json"
    if summary.exists():
        previous = json.loads(summary.read_text())
        if previous.get("objective", "td") != args.objective:
            raise ValueError("output directory belongs to another objective; use separate --output directories")
    for path in (args.output / "layers").glob("*.json"):
        if json.loads(path.read_text()).get("model_version") != model.layer_version(args.objective):
            raise ValueError(f"output directory contains another arithmetic model: {path}; use separate --output directories")


def terminate_children():
    STOP.set()
    with CHILDREN_LOCK:
        for process in tuple(CHILDREN):
            if process.poll() is None:
                try:
                    os.killpg(process.pid, signal.SIGTERM)
                except ProcessLookupError:
                    pass


def run_layer(task, args, identity):
    n, kind, k = task
    destination = cache_path(args.output, task)
    lock = destination.with_suffix(".lock")
    lock.parent.mkdir(parents=True, exist_ok=True)
    try:
        descriptor = os.open(lock, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o600)
    except FileExistsError:
        raise RuntimeError(f"layer is locked: {lock}; see README for interrupted-run recovery")
    try:
        with os.fdopen(descriptor, "w") as handle:
            handle.write(f"pid={os.getpid()} host={os.uname().nodename}\n")
        command = [str(EXECUTABLE), "--n", str(n), "--kind", kind, "--lanes", str(k),
                   "--data-root", str(ROOT / "data"),
                   "--architecture", model.architecture(args.objective),
                   "--curve-a-hex", hex(poly_to_l(CURVES[n].a, n))]
        if args.verify:
            fixture = args.output / "fixtures" / f"n{n}_k{k}.txt"
            fixture.parent.mkdir(parents=True, exist_ok=True)
            # Separate kind-specific filenames avoid concurrent fixture writes.
            fixture = fixture.with_name(f"n{n}_{kind}_k{k}.txt")
            fixture.write_text(make_fixture(n, k), encoding="utf-8")
            command += ["--verify-fixture", str(fixture)]
        if args.job_memory_gb:
            command += ["--memory-limit-bytes", str(int(args.job_memory_gb * 1024**3))]
        if STOP.is_set():
            raise RuntimeError("run interrupted")
        log = args.output / "logs" / f"n{n}_{kind}_k{k}.log"
        log.parent.mkdir(parents=True, exist_ok=True)
        with log.open("w", encoding="utf-8") as stderr:
            with CHILDREN_LOCK:
                if STOP.is_set():
                    raise RuntimeError("run interrupted")
                process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=stderr,
                                           text=True, start_new_session=True)
                CHILDREN.add(process)
            try:
                stdout, _ = process.communicate(timeout=args.timeout or None)
            except subprocess.TimeoutExpired:
                os.killpg(process.pid, signal.SIGKILL)
                process.communicate()
                raise RuntimeError(f"layer timed out: {task}; see {log}")
            finally:
                with CHILDREN_LOCK:
                    CHILDREN.discard(process)
        if process.returncode:
            raise RuntimeError(f"layer {task} failed with exit {process.returncode}; see {log}")
        record = json.loads(stdout)
        model.validate_record(record, *task, args.objective)
        if args.verify and not record.get("basis_verified"):
            raise RuntimeError("backend did not verify requested fixture")
        record["fingerprint"] = identity
        record["qrom_model"] = "not_applicable_arithmetic_layer"
        atomic_json(destination, record)
        return record
    finally:
        lock.unlink(missing_ok=True)


def generate(args):
    check_output_architecture(args)
    manifest()
    identity = layer_fingerprint()
    print("Windowing: lane_first; p = manuscript w; scanning 1 <= p <= 2n+2.", flush=True)
    print(f"Objective: {args.objective}; point addition: {model.architecture(args.objective)}.", flush=True)
    tasks, candidates = task_plan(args)
    pending = []
    for task in tasks:
        path = cache_path(args.output, task)
        if path.exists():
            record = read_record(path, task, identity, args.objective)
            if not args.verify or record.get("basis_verified"):
                continue
        pending.append(task)
    print(f"Layers: {len(tasks)}, cached: {len(tasks)-len(pending)}, pending: {len(pending)}", flush=True)
    if args.memory_budget_gb and any(estimated_gb(task, args.objective) > args.memory_budget_gb for task in pending):
        raise ValueError("a task exceeds the estimated memory budget; increase --memory-budget-gb or use --scope winners")
    if pending:
        ensure_build(args.jobs)
    atomic_json(args.output / "last_plan.json", dict(model_version=model.model_version(args.objective),
        objective=args.objective, point_addition=model.architecture(args.objective),
        fingerprint=fingerprint(args.objective), layer_fingerprint=identity,
        field_sizes=args.n, scope=args.scope, candidates=candidates,
        tasks=[list(task) for task in tasks], pending=[list(task) for task in pending]))
    pool = futures.ThreadPoolExecutor(max_workers=args.jobs)
    active = {}
    used = 0.0
    failed = []
    completed = 0
    try:
        while pending or active:
            for task in list(pending):
                estimate = estimated_gb(task, args.objective)
                if len(active) >= args.jobs:
                    break
                if args.memory_budget_gb and used + estimate > args.memory_budget_gb:
                    continue
                pending.remove(task)
                job = pool.submit(run_layer, task, args, identity)
                active[job] = (task, estimate)
                used += estimate
                print(f"start {task}, estimated {estimate:.1f} GiB", flush=True)
            done, _ = futures.wait(active, return_when=futures.FIRST_COMPLETED, timeout=1)
            for job in done:
                task, estimate = active.pop(job)
                used -= estimate
                try:
                    record = job.result()
                    completed += 1
                    print(f"done {task}: TD={record['toffoli_depth']}, "
                          f"{record['seconds']:.1f}s, peak RSS={record['peak_rss_bytes']/1024**3:.2f} GiB", flush=True)
                except Exception as error:
                    failed.append(dict(task=list(task), error=str(error)))
                    print(f"FAILED {task}: {error}", file=sys.stderr, flush=True)
            atomic_json(args.output / "run_progress.json", dict(completed_this_run=completed,
                failed=failed, pending=len(pending), running=len(active)))
    except KeyboardInterrupt:
        terminate_children()
        raise
    finally:
        pool.shutdown(wait=True)
    scan(args)
    if failed:
        raise RuntimeError(f"{len(failed)} layers failed; completed layers retained; rerun the same command")


def scan(args):
    check_output_architecture(args)
    identity = layer_fingerprint()
    records = {}
    for path in sorted((args.output / "layers").glob("*.json")):
        raw = json.loads(path.read_text())
        task = raw["n"], raw["kind"], raw["lanes"]
        records[task] = read_record(path, task, identity, args.objective)
    rows = []
    summaries = {}
    for n in args.n:
        local = [model.evaluate(n, s, p, records, args.objective) for s, p in model.grid(n, args.s_min, args.s_max)]
        rows.extend(local)
        available = [row for row in local if row["available"]]
        candidates = model.formula_minimizers(n, args.s_min, args.s_max, args.objective)
        td_key = lambda row: (row["toffoli_depth"], row["width"], row["toffoli"], row["s"], row["p"])
        tdw_key = lambda row: (row["tdw"], row["toffoli_depth"], row["width"], row["toffoli"], row["s"], row["p"])
        best_td = min(available, key=td_key) if available else None
        best_tdw = min(available, key=tdw_key) if available else None
        best = best_td if args.objective == "td" else best_tdw
        by_key = {(row["s"], row["p"]): row for row in local}
        summaries[str(n)] = dict(grid_points=len(local), available_points=len(available),
            complete_full_resource_grid=len(local)==len(available), formula_minimizers=candidates,
            formula_minimizers_all_available=all(by_key[row["s"],row["p"]]["available"] for row in candidates),
            minimum_toffoli_depth_over_available=best_td,
            minimum_tdw_over_available=best_tdw, selected_optimum=best)
    columns = ["n", "s", "p", "w", "original_addends", "addends", "accumulation_layers",
               "tree_layers", "windowing_model", "objective", "point_addition", "available", "missing_layers",
               "predicted_toffoli_depth", "predicted_width", "predicted_tdw",
               "toffoli", "cnot", "width", "toffoli_depth", "nct_depth_barrier", "dw", "tdw",
               "forward_retained_qubits", "qrom_model", "depth_model", "width_model"]
    args.output.mkdir(parents=True, exist_ok=True)
    temporary = args.output / f"scan.{os.getpid()}.csv.tmp"
    with temporary.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader(); writer.writerows(rows)
    os.replace(temporary, args.output / "scan.csv")
    summary = dict(model_version=model.model_version(args.objective), fingerprint=fingerprint(args.objective),
                   objective=args.objective, point_addition=model.architecture(args.objective),
                   arithmetic_model_version=model.layer_version(args.objective), layer_fingerprint=identity,
                   windowing_model="lane_first", parameter_p="manuscript w",
                   search_range=dict(s_min=args.s_min, s_max=args.s_max, p="1..2n+2"),
                   scope="module-composed estimates; not monolithic gate-stream synthesis",
                   qrom="table-independent h_u=2n", fields=summaries)
    atomic_json(args.output / "summary.json", summary)
    for n, info in summaries.items():
        best = info["selected_optimum"]
        print(f"n={n}: {info['available_points']}/{info['grid_points']} complete resource rows; "
              + (f"best available by {args.objective.upper()} (s,p)=({best['s']},{best['p']}), "
                 f"TD={best['toffoli_depth']}, width={best['width']}, TDW={best['tdw']}"
                 if best else "no complete resource rows")
              + f"; all {args.objective.upper()} formula minimizers available: {info['formula_minimizers_all_available']}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    build = sub.add_parser("build")
    build.add_argument("--jobs", type=int, default=2)
    sub.add_parser("verify-data")
    sub.add_parser("create-manifest", help=argparse.SUPPRESS)
    for name in ("plan", "run", "scan"):
        child = sub.add_parser(name)
        child.add_argument("--n", type=int, nargs="+", choices=model.FIELDS, default=[163])
        child.add_argument("--s-min", type=int, default=2)
        child.add_argument("--s-max", type=int, default=18)
        child.add_argument("--scope", choices=("smoke", "winners", "full"), default="winners")
        child.add_argument("--objective", choices=("td", "tdw"), default="td",
                           help="td: retained point addition; tdw: Balanced compute-copy-uncompute")
        child.add_argument("--output", type=Path, default=None,
                           help="default: results for td, results_tdw for tdw; never mix architectures")
        if name == "run":
            child.add_argument("--jobs", type=int, default=1)
            child.add_argument("--memory-budget-gb", type=float, default=0)
            child.add_argument("--job-memory-gb", type=float, default=0,
                               help="Linux RLIMIT_AS per process; 0 disables the hard cap")
            child.add_argument("--timeout", type=float, default=0, help="seconds per layer; 0 disables timeout")
            child.add_argument("--verify", action="store_true", help="basis-state and clean-scratch checks")
    args = parser.parse_args()
    if hasattr(args, "jobs") and args.jobs < 1:
        parser.error("--jobs must be positive")
    if hasattr(args, "n"):
        args.n = sorted(set(args.n))
        args.output = (args.output or ROOT / ("results" if args.objective == "td" else "results_tdw")).resolve()
        if not 2 <= args.s_min <= args.s_max <= 18:
            parser.error("require 2 <= s-min <= s-max <= 18")
    for name in ("memory_budget_gb", "job_memory_gb", "timeout"):
        if getattr(args, name, 0) < 0:
            parser.error(f"{name} must be nonnegative")
    if getattr(args, "job_memory_gb", 0) and sys.platform != "linux":
        parser.error("--job-memory-gb is supported on Linux only")
    if args.command in ("verify-data", "create-manifest"):
        manifest(create=args.command == "create-manifest")
    elif args.command == "build":
        manifest()
        build_project(args.jobs)
    elif args.command == "plan":
        tasks, candidates = task_plan(args)
        print(json.dumps(dict(model_version=model.model_version(args.objective), windowing_model="lane_first",
            objective=args.objective, point_addition=model.architecture(args.objective),
            parameter_p="manuscript w", scope=args.scope, task_count=len(tasks), candidates=candidates,
            maximum_estimated_job_gb=max((estimated_gb(task,args.objective) for task in tasks), default=0),
            tasks=[dict(n=n, kind=kind, lanes=k, estimated_gb=estimated_gb((n,kind,k),args.objective)) for n,kind,k in tasks]), indent=2))
    elif args.command == "run":
        generate(args)
    else:
        scan(args)


if __name__ == "__main__":
    def on_termination(signum, frame):
        raise KeyboardInterrupt
    signal.signal(signal.SIGTERM, on_termination)
    try:
        main()
    except KeyboardInterrupt:
        terminate_children()
        print("Interrupted; completed layer JSON files are reusable.", file=sys.stderr)
        raise SystemExit(130)
    except Exception as error:
        print(f"error: {error}", file=sys.stderr)
        raise SystemExit(1)
