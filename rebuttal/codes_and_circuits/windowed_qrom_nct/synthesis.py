"""Gate-exact synthesis and Table-6-style aggregation for one (n,s,p) point."""

from __future__ import annotations

import hashlib
import json
import csv
import subprocess
from dataclasses import asdict, dataclass
from pathlib import Path

from .basis import poly_to_l
from .curves import CURVES
from .paths import (
    AC_DATA_ROOT,
    AC_EXECUTABLE,
    AC_RESULTS,
    AC_ROOT,
    PACKAGE_ROOT,
)
from .tables import WindowDescriptor, descriptor, write_basis_table


DEFAULT_CACHE = PACKAGE_ROOT / "artifacts" / "cache"
AC_EXE = AC_EXECUTABLE
AC_DATA = AC_DATA_ROOT
GATE_FLOW_VERSION = "window-qrom-early-unlookup-fanout-v1"


@dataclass(frozen=True)
class ExactLayer:
    qubits: int
    toffoli: int
    cnot: int
    full_depth: int
    current_depth: int
    toffoli_depth: int


@dataclass(frozen=True)
class ExactDesign:
    model_version: str
    arithmetic: str
    n: int
    scalar_bits: int
    s: int
    p: int
    addends: int
    phase1_layers: int
    tree_layers: int
    toffoli: int
    cnot: int
    width: int
    toffoli_depth: int
    nct_depth: int
    current_depth: int
    dw: int
    tdw: int
    table_seed: str
    data_basis: str
    qrom_data_write: str
    scheduler: str
    layer_records: tuple[dict[str, object], ...]


def scheduled_descriptors(n: int, s: int) -> list[WindowDescriptor]:
    """Interleave P_j,Q_j and put both short top windows last."""
    curve = CURVES[n]
    windows = (curve.n + 1 + s - 1) // s
    result: list[WindowDescriptor] = []
    for j in range(windows):
        result.append(descriptor(curve, s, "P", j))
        result.append(descriptor(curve, s, "Q", j))
    return result


def _table_path(desc: WindowDescriptor, basis: str, cache: Path) -> Path:
    return cache / "tables" / basis / f"{desc.table_id}.bin"


def ensure_table(desc: WindowDescriptor, basis: str, cache: Path = DEFAULT_CACHE) -> Path:
    if basis not in {"polynomial", "l"}:
        raise ValueError("basis must be polynomial or l")
    path = _table_path(desc, basis, cache)
    expected = desc.entries * 2 * CURVES[desc.n].coordinate_bytes
    if path.is_file() and path.stat().st_size == expected:
        return path
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("wb") as output:
        write_basis_table(desc, basis, output)
    temporary.replace(path)
    return path


def _table_spec(descriptors: list[WindowDescriptor], basis: str, cache: Path) -> str:
    return ",".join(
        f"{item.window_bits}:{ensure_table(item, basis, cache)}" for item in descriptors
    )


def _run_json(command: list[str], *, cwd: Path) -> ExactLayer:
    executable = Path(command[0])
    if not executable.is_file():
        raise FileNotFoundError(
            f"estimator executable not found: {executable}; "
            "build the bundled AC backend as documented in README.md"
        )
    process = subprocess.run(command, cwd=cwd, text=True, capture_output=True, check=False)
    if process.returncode:
        raise RuntimeError(
            f"gate synthesis failed ({process.returncode}): {' '.join(command)}\n{process.stderr}"
        )
    payload = json.loads(process.stdout.strip().splitlines()[-1])
    return ExactLayer(**{field: int(payload[field]) for field in ExactLayer.__dataclass_fields__})


def _cache_key(kind: str, arithmetic: str, n: int, ids: list[str]) -> str:
    payload = json.dumps(
        {"version": GATE_FLOW_VERSION, "kind": kind, "arithmetic": arithmetic, "n": n, "ids": ids},
        sort_keys=True,
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def synthesize_qrom_layer(
    n: int,
    arithmetic: str,
    tables: list[WindowDescriptor],
    *,
    initialization: bool,
    cache: Path = DEFAULT_CACHE,
) -> ExactLayer:
    if arithmetic != "ac":
        raise ValueError("this submission artifact supports only AC arithmetic")
    kind = "initialization" if initialization else "early-window-add"
    key = _cache_key(kind, arithmetic, n, [table.table_id for table in tables])
    record = cache / "layers" / arithmetic / f"n{n}-{key}.json"
    if record.is_file():
        return ExactLayer(**json.loads(record.read_text(encoding="utf-8"))["stats"])

    # Initialization QROMs act on pairwise-disjoint address, target, selector,
    # and fanout registers.  Per-line ASAP therefore composes them exactly by
    # summing width/counts and taking the maximum depth.  Synthesizing each
    # table separately avoids retaining hundreds of millions of unrelated
    # Clifford gates at once for p=L.
    if initialization and len(tables) > 1:
        parts = [
            synthesize_qrom_layer(n, arithmetic, [table], initialization=True, cache=cache)
            for table in tables
        ]
        stats = ExactLayer(
            qubits=sum(part.qubits for part in parts),
            toffoli=sum(part.toffoli for part in parts),
            cnot=sum(part.cnot for part in parts),
            full_depth=max(part.full_depth for part in parts),
            current_depth=max(part.current_depth for part in parts),
            toffoli_depth=max(part.toffoli_depth for part in parts),
        )
        record.parent.mkdir(parents=True, exist_ok=True)
        record.write_text(
            json.dumps(
                {"model_version": GATE_FLOW_VERSION, "kind": kind, "arithmetic": arithmetic,
                 "composition": "exact-disjoint-register-sum-max", "tables": [x.table_id for x in tables],
                 "stats": asdict(stats)},
                indent=2, sort_keys=True,
            ) + "\n",
            encoding="utf-8",
        )
        return stats

    basis = "l"
    specification = _table_spec(tables, basis, cache)
    command = [
        str(AC_EXE),
        "window-init" if initialization else "window-layer",
        "--n",
        str(n),
        "--tables",
        specification,
        "--data-root",
        str(AC_DATA),
    ]
    if not initialization:
        command.extend(("--curve-a-hex", hex(poly_to_l(CURVES[n].a, n))))
    cwd = AC_ROOT
    stats = _run_json(command, cwd=cwd)
    record.parent.mkdir(parents=True, exist_ok=True)
    record.write_text(
        json.dumps(
            {
                "model_version": GATE_FLOW_VERSION,
                "kind": kind,
                "arithmetic": arithmetic,
                "tables": [table.table_id for table in tables],
                "stats": asdict(stats),
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )
    return stats


def synthesize_tree_layer(
    n: int, arithmetic: str, lanes: int, *, cache: Path = DEFAULT_CACHE
) -> ExactLayer:
    if arithmetic != "ac":
        raise ValueError("this submission artifact supports only AC arithmetic")
    key = _cache_key("tree", arithmetic, n, [str(lanes)])
    record = cache / "layers" / arithmetic / f"n{n}-{key}.json"
    if record.is_file():
        return ExactLayer(**json.loads(record.read_text(encoding="utf-8"))["stats"])
    candidates = sorted(AC_RESULTS.glob(f"layers_n{n}*.csv"))
    # Prefer the most complete cache when multiple cache fragments coexist.
    for path in candidates:
        if not path.is_file():
            continue
        with path.open(newline="", encoding="utf-8") as handle:
            for row in csv.DictReader(handle):
                if row.get("mode") != "reduction" or int(row["lanes"]) != lanes:
                    continue
                stats = ExactLayer(
                    qubits=int(row["qubits"]), toffoli=int(row["toffoli"]),
                    cnot=int(row["cnot"]), full_depth=int(row["full_depth"]),
                    current_depth=int(row["current_depth"]),
                    toffoli_depth=int(row["toffoli_depth"]),
                )
                record.parent.mkdir(parents=True, exist_ok=True)
                record.write_text(
                    json.dumps(
                        {"model_version": GATE_FLOW_VERSION, "kind": "tree",
                         "arithmetic": arithmetic, "lanes": lanes,
                         "source_cache": str(path.relative_to(PACKAGE_ROOT)),
                         "stats": asdict(stats)},
                        indent=2, sort_keys=True,
                    ) + "\n",
                    encoding="utf-8",
                )
                return stats
    command = [
        str(AC_EXE), "layer", "--n", str(n), "--lanes", str(lanes),
        "--mode", "reduction", "--data-root", str(AC_DATA),
    ]
    cwd = AC_ROOT
    stats = _run_json(command, cwd=cwd)
    record.parent.mkdir(parents=True, exist_ok=True)
    record.write_text(
        json.dumps(
            {"model_version": GATE_FLOW_VERSION, "kind": "tree", "arithmetic": arithmetic,
             "lanes": lanes, "stats": asdict(stats)},
            indent=2, sort_keys=True,
        ) + "\n",
        encoding="utf-8",
    )
    return stats


def synthesize_design(
    n: int, s: int, p: int, arithmetic: str, *, cache: Path = DEFAULT_CACHE
) -> ExactDesign:
    if arithmetic != "ac":
        raise ValueError("this submission artifact supports only AC arithmetic")
    windows = scheduled_descriptors(n, s)
    lane_count = min(p, len(windows))
    if p < 1 or p > len(windows):
        raise ValueError(f"p must be in 1..{len(windows)}")
    m = n + 1
    records: list[dict[str, object]] = []

    init_tables = windows[:lane_count]
    init = synthesize_qrom_layer(n, arithmetic, init_tables, initialization=True, cache=cache)
    active_address = sum(table.window_bits for table in init_tables)
    # All scalar bits exist globally; the standalone layer includes only the
    # active address slices, so replace those by the complete 2m-bit input.
    peak = 2 * m + init.qubits - active_address
    clear = init.qubits - active_address - 2 * n * lane_count
    forward_t = init.toffoli
    forward_c = init.cnot
    forward_d = init.full_depth
    forward_cd = init.current_depth
    forward_td = init.toffoli_depth
    records.append({"phase": "initialization", "tables": [x.table_id for x in init_tables],
                    "stats": asdict(init)})

    phase1_layers = 0
    for start in range(lane_count, len(windows), lane_count):
        batch = windows[start : start + lane_count]
        k = len(batch)
        layer = synthesize_qrom_layer(n, arithmetic, batch, initialization=False, cache=cache)
        address_bits = sum(table.window_bits for table in batch)
        demand = layer.qubits - address_bits - 2 * n * k
        if demand > clear:
            peak += demand - clear
            clear = layer.qubits - address_bits - 7 * n * k
        else:
            clear = layer.qubits - address_bits - 7 * n * k + (clear - demand)
        forward_t += layer.toffoli
        forward_c += layer.cnot
        forward_d += layer.full_depth
        forward_cd += layer.current_depth
        forward_td += layer.toffoli_depth
        phase1_layers += 1
        records.append({"phase": "phase1", "tables": [x.table_id for x in batch],
                        "stats": asdict(layer)})

    tree_layers = 0
    items = lane_count
    while items > 1:
        k = items // 2
        layer = synthesize_tree_layer(n, arithmetic, k, cache=cache)
        demand = layer.qubits - 4 * n * k
        if demand > clear:
            peak += demand - clear
            clear = layer.qubits - 8 * n * k
        else:
            clear = layer.qubits - 8 * n * k + (clear - demand)
        forward_t += layer.toffoli
        forward_c += layer.cnot
        forward_d += layer.full_depth
        forward_cd += layer.current_depth
        forward_td += layer.toffoli_depth
        tree_layers += 1
        records.append({"phase": "tree", "lanes": k, "stats": asdict(layer)})
        items = k + items % 2

    width = peak + 2 * n
    toffoli = 2 * forward_t
    cnot = 2 * forward_c + 2 * n
    nct_depth = 2 * forward_d + 2
    current_depth = 2 * forward_cd + 2
    toffoli_depth = 2 * forward_td
    return ExactDesign(
        model_version=GATE_FLOW_VERSION,
        arithmetic=arithmetic,
        n=n,
        scalar_bits=m,
        s=s,
        p=p,
        addends=len(windows),
        phase1_layers=phase1_layers,
        tree_layers=tree_layers,
        toffoli=toffoli,
        cnot=cnot,
        width=width,
        toffoli_depth=toffoli_depth,
        nct_depth=nct_depth,
        current_depth=current_depth,
        dw=nct_depth * width,
        tdw=toffoli_depth * width,
        table_seed="de0faddf37eaf09c24fe74b8dc130a5f75e73234f62f1ca55dc66b7607a67340",
        data_basis="AC L-basis",
        qrom_data_write="balanced CNOT fanout; 2n-1 reusable fanout lines per active lane",
        scheduler="repository GateManager per-line ASAP after self-inverse gate cancellation",
        layer_records=tuple(records),
    )


def write_design(path: Path, design: ExactDesign) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(asdict(design), indent=2, sort_keys=True) + "\n", encoding="utf-8")
