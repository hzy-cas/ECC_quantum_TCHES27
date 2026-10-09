"""Frozen shifted-window tables and canonical manifest records."""

from __future__ import annotations

import hashlib
import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import BinaryIO, Iterator

from .basis import DEFAULT_DATA_ROOT, poly_to_l
from .curves import CURVES, BinaryCurve
from .openssl_ec import PointGenerator


SEED_LABEL = b"Binary_ECC/windowed-QROM/shifted-tables/v1"
PUBLIC_SEED = hashlib.sha256(SEED_LABEL).digest()
DOMAIN_TAG = "binary-ecc-window-shift-v1"
GENERATOR_VERSION = "shifted-table-generator-v1"


@dataclass(frozen=True)
class WindowDescriptor:
    n: int
    scalar_bits: int
    window_size: int
    stream: str
    window_index: int
    window_bits: int
    entries: int
    base_scalar: int | None
    step_scalar: int
    delta_scalar: int
    scalar_reference: str
    base_x: int
    base_y: int
    public_q_tag: str
    derivation_attempt: int

    @property
    def table_id(self) -> str:
        base = (
            f"n{self.n}-s{self.window_size}-{self.stream}-"
            f"j{self.window_index}-b{self.window_bits}"
        )
        if self.stream == "Q" and self.public_q_tag != "q7p":
            return f"{base}-q{self.public_q_tag}"
        return base


@dataclass(frozen=True)
class TableStats:
    table_id: str
    entries: int
    coordinate_bytes: int
    polynomial_sha256: str
    l_basis_sha256: str
    polynomial_hamming_weight: int
    l_basis_hamming_weight: int
    delta_x: str
    delta_y: str


def _derive_candidate(curve: BinaryCurve, s: int, stream: str, j: int, attempt: int) -> int:
    payload = b"\x00".join(
        (
            PUBLIC_SEED,
            DOMAIN_TAG.encode("ascii"),
            str(curve.n).encode("ascii"),
            str(s).encode("ascii"),
            stream.encode("ascii"),
            str(j).encode("ascii"),
            str(attempt).encode("ascii"),
        )
    )
    width = (curve.order.bit_length() + 7) // 8 + 16
    return int.from_bytes(hashlib.shake_256(payload).digest(width), "big") % curve.order


def _default_q(curve: BinaryCurve) -> tuple[int, int]:
    with PointGenerator(curve) as generator:
        [(qx, qy)] = list(generator.progression(7, 1, 1))
    return qx, qy


def _public_q_tag(curve: BinaryCurve, q_point: tuple[int, int] | None) -> str:
    if q_point is None:
        return "q7p"
    width = curve.coordinate_bytes
    payload = q_point[0].to_bytes(width, "big") + q_point[1].to_bytes(width, "big")
    return hashlib.sha256(payload).hexdigest()[:16]


def validate_public_q(curve: BinaryCurve, q_point: tuple[int, int]) -> None:
    with PointGenerator(curve) as generator:
        generator.validate_prime_subgroup_point(*q_point)


def descriptor(
    curve: BinaryCurve,
    s: int,
    stream: str,
    j: int,
    q_point: tuple[int, int] | None = None,
) -> WindowDescriptor:
    if not 2 <= s <= 18:
        raise ValueError("window size must be in 2..18")
    if stream not in {"P", "Q"}:
        raise ValueError("stream must be P or Q")
    scalar_bits = curve.n + 1
    windows = (scalar_bits + s - 1) // s
    if not 0 <= j < windows:
        raise ValueError("window index out of range")
    top_bits = scalar_bits - s * (windows - 1)
    bits = top_bits if j == windows - 1 else s
    entries = 1 << bits
    if q_point is not None:
        validate_public_q(curve, q_point)
    if stream == "P":
        base_scalar: int | None = 1
        scalar_reference = "P"
        base_x, base_y = curve.gx, curve.gy
        step = pow(2, s * j, curve.order)
    elif q_point is None:
        # Preserve the frozen Q=[7]P test instance and its original tables.
        base_scalar = 7
        scalar_reference = "P"
        base_x, base_y = curve.gx, curve.gy
        step = (7 * pow(2, s * j, curve.order)) % curve.order
    else:
        # No discrete logarithm of Q relative to P is needed: both the shift
        # and the window step are expressed as known multiples of public Q.
        base_scalar = None
        scalar_reference = "Q"
        base_x, base_y = q_point
        step = pow(2, s * j, curve.order)
    inverse_step = pow(step, -1, curve.order)
    attempt = 0
    while True:
        delta = _derive_candidate(curve, s, stream, j, attempt)
        # Exactly one u modulo r could make delta+u*step zero.  Reject if it is
        # one of the represented table addresses; also reject delta=0 directly.
        bad_u = (-delta * inverse_step) % curve.order
        if delta and bad_u >= entries:
            break
        attempt += 1
    return WindowDescriptor(
        n=curve.n,
        scalar_bits=scalar_bits,
        window_size=s,
        stream=stream,
        window_index=j,
        window_bits=bits,
        entries=entries,
        base_scalar=base_scalar,
        step_scalar=step,
        delta_scalar=delta,
        scalar_reference=scalar_reference,
        base_x=base_x,
        base_y=base_y,
        public_q_tag=_public_q_tag(curve, q_point),
        derivation_attempt=attempt,
    )


def descriptors(
    n: int, s: int, q_point: tuple[int, int] | None = None
) -> list[WindowDescriptor]:
    curve = CURVES[n]
    if q_point is not None:
        validate_public_q(curve, q_point)
    count = (curve.n + 1 + s - 1) // s
    # Keep P and Q windows separate and in a fixed, documented order.
    return [
        descriptor(curve, s, stream, j, q_point)
        for stream in ("P", "Q")
        for j in range(count)
    ]


def iter_table(desc: WindowDescriptor) -> Iterator[tuple[int, int]]:
    curve = CURVES[desc.n]
    with PointGenerator(curve) as generator:
        yield from generator.affine_progression(
            desc.base_x,
            desc.base_y,
            desc.delta_scalar,
            desc.step_scalar,
            desc.entries,
        )


def _encoded_pair(x: int, y: int, width: int) -> bytes:
    return x.to_bytes(width, "big") + y.to_bytes(width, "big")


def materialize_table(
    desc: WindowDescriptor,
    *,
    polynomial_output: BinaryIO | None = None,
    l_basis_output: BinaryIO | None = None,
    data_root: Path = DEFAULT_DATA_ROOT,
) -> TableStats:
    curve = CURVES[desc.n]
    width = curve.coordinate_bytes
    poly_hash = hashlib.sha256()
    l_hash = hashlib.sha256()
    poly_weight = 0
    l_weight = 0
    delta_xy: tuple[int, int] | None = None
    for index, (x, y) in enumerate(iter_table(desc)):
        if index == 0:
            delta_xy = (x, y)
        encoded = _encoded_pair(x, y, width)
        poly_hash.update(encoded)
        poly_weight += x.bit_count() + y.bit_count()
        if polynomial_output is not None:
            polynomial_output.write(encoded)

        lx = poly_to_l(x, curve.n, data_root)
        ly = poly_to_l(y, curve.n, data_root)
        l_encoded = _encoded_pair(lx, ly, width)
        l_hash.update(l_encoded)
        l_weight += lx.bit_count() + ly.bit_count()
        if l_basis_output is not None:
            l_basis_output.write(l_encoded)
    if delta_xy is None:
        raise AssertionError("a table must contain at least one entry")
    return TableStats(
        table_id=desc.table_id,
        entries=desc.entries,
        coordinate_bytes=width,
        polynomial_sha256=poly_hash.hexdigest(),
        l_basis_sha256=l_hash.hexdigest(),
        polynomial_hamming_weight=poly_weight,
        l_basis_hamming_weight=l_weight,
        delta_x=f"0x{delta_xy[0]:0{2 * width}x}",
        delta_y=f"0x{delta_xy[1]:0{2 * width}x}",
    )


def write_basis_table(
    desc: WindowDescriptor,
    basis: str,
    output: BinaryIO,
    *,
    data_root: Path = DEFAULT_DATA_ROOT,
) -> None:
    """Write one basis without computing the unused basis' manifest statistics."""
    if basis not in {"polynomial", "l"}:
        raise ValueError("basis must be polynomial or l")
    curve = CURVES[desc.n]
    width = curve.coordinate_bytes
    for x, y in iter_table(desc):
        if basis == "l":
            x = poly_to_l(x, curve.n, data_root)
            y = poly_to_l(y, curve.n, data_root)
        output.write(_encoded_pair(x, y, width))


def curve_manifest(
    curve: BinaryCurve, q_point: tuple[int, int] | None = None
) -> dict[str, object]:
    width = 2 * curve.coordinate_bytes
    if q_point is None:
        qx, qy = _default_q(curve)
        q_scalar: int | None = 7
        q_source = "frozen reproducible test point Q=[7]P"
    else:
        validate_public_q(curve, q_point)
        qx, qy = q_point
        q_scalar = None
        q_source = "user-supplied public affine coordinates"
    return {
        "name": curve.openssl_name,
        "n": curve.n,
        "field_modulus": hex(curve.modulus),
        "a": hex(curve.a),
        "b": f"0x{curve.b:0{width}x}",
        "P": {"x": f"0x{curve.gx:0{width}x}", "y": f"0x{curve.gy:0{width}x}"},
        "Q": {"x": f"0x{qx:0{width}x}", "y": f"0x{qy:0{width}x}"},
        "Q_scalar_relative_to_P": q_scalar,
        "Q_source": q_source,
        "order": hex(curve.order),
        "cofactor": curve.cofactor,
    }


def write_manifest(
    path: Path,
    *,
    n: int,
    s: int,
    include_table_stats: bool = True,
    data_root: Path = DEFAULT_DATA_ROOT,
    q_point: tuple[int, int] | None = None,
) -> dict[str, object]:
    curve = CURVES[n]
    records = []
    table_descriptors = descriptors(n, s, q_point)
    for desc in table_descriptors:
        record: dict[str, object] = asdict(desc)
        record["table_id"] = desc.table_id
        if include_table_stats:
            record["table_stats"] = asdict(materialize_table(desc, data_root=data_root))
        records.append(record)
    with PointGenerator(curve) as generator:
        delta_points = []
        for desc in table_descriptors:
            point = generator.multiple_coordinates(
                desc.base_x, desc.base_y, desc.delta_scalar
            )
            if point is None:
                raise AssertionError("a shifted-table offset must be finite")
            delta_points.append(point)
        total_point = generator.sum_coordinates(iter(delta_points))
    coordinate_width = 2 * curve.coordinate_bytes
    if total_point is None:
        delta_total: dict[str, object] = {"point": "infinity"}
    else:
        delta_total = {
            "x": f"0x{total_point[0]:0{coordinate_width}x}",
            "y": f"0x{total_point[1]:0{coordinate_width}x}",
        }
    if q_point is None:
        delta_total["scalar_relative_to_P"] = (
            sum(desc.delta_scalar for desc in table_descriptors) % curve.order
        )
    else:
        delta_total["scalar_relative_to_P"] = None
    manifest: dict[str, object] = {
        "schema": "binary-ecc-window-tables-manifest-v1",
        "generator_version": GENERATOR_VERSION,
        "public_seed_hex": PUBLIC_SEED.hex(),
        "public_seed_label": SEED_LABEL.decode("ascii"),
        "domain_tag": DOMAIN_TAG,
        "bit_order": "field bit i is qubit i (little endian); canonical hashes use fixed-width big-endian bytes",
        "curve": curve_manifest(curve, q_point),
        "window_size": s,
        "descriptor_order": "all P windows by ascending j, then all Q windows by ascending j",
        "delta_total": delta_total,
        "tables": records,
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return manifest


def verify_manifest(path: Path, *, verify_table_stats: bool = False) -> dict[str, object]:
    """Fail closed if a manifest no longer matches the frozen derivation."""
    payload = json.loads(path.read_text(encoding="utf-8"))
    if payload.get("schema") != "binary-ecc-window-tables-manifest-v1":
        raise ValueError(f"unknown manifest schema: {path}")
    if payload.get("public_seed_hex") != PUBLIC_SEED.hex():
        raise ValueError(f"public seed mismatch: {path}")
    curve_record = payload.get("curve")
    if not isinstance(curve_record, dict):
        raise ValueError(f"missing curve record: {path}")
    n = int(curve_record["n"])
    s = int(payload["window_size"])
    q_record = curve_record.get("Q")
    if not isinstance(q_record, dict):
        raise ValueError(f"missing public Q record: {path}")
    q_coordinates = (int(q_record["x"], 16), int(q_record["y"], 16))
    q_point = None if curve_record.get("Q_scalar_relative_to_P") == 7 else q_coordinates
    if curve_record != curve_manifest(CURVES[n], q_point):
        raise ValueError(f"curve constants mismatch: {path}")
    actual = payload.get("tables")
    if not isinstance(actual, list):
        raise ValueError(f"missing table list: {path}")
    expected = descriptors(n, s, q_point)
    if len(actual) != len(expected):
        raise ValueError(f"descriptor count mismatch: {path}")
    for record, desc in zip(actual, expected):
        if not isinstance(record, dict):
            raise ValueError(f"malformed descriptor: {path}")
        expected_descriptor = asdict(desc)
        expected_descriptor["table_id"] = desc.table_id
        for key, value in expected_descriptor.items():
            if record.get(key) != value:
                raise ValueError(f"descriptor {desc.table_id} field {key} mismatch")
        if verify_table_stats:
            if "table_stats" not in record:
                raise ValueError(f"{desc.table_id} has no table statistics")
            regenerated = asdict(materialize_table(desc))
            if record["table_stats"] != regenerated:
                raise ValueError(f"table hash/statistics mismatch: {desc.table_id}")
    expected_total = payload.get("delta_total")
    if not isinstance(expected_total, dict):
        raise ValueError(f"delta_total mismatch: {path}")
    with PointGenerator(CURVES[n]) as generator:
        delta_points = []
        for desc in expected:
            point = generator.multiple_coordinates(
                desc.base_x, desc.base_y, desc.delta_scalar
            )
            if point is None:
                raise ValueError(f"infinite table offset: {desc.table_id}")
            delta_points.append(point)
        total_point = generator.sum_coordinates(iter(delta_points))
    width = 2 * CURVES[n].coordinate_bytes
    if total_point is None:
        if expected_total.get("point") != "infinity":
            raise ValueError(f"zero delta_total is not marked as infinity: {path}")
    elif (
        expected_total.get("x") != f"0x{total_point[0]:0{width}x}"
        or expected_total.get("y") != f"0x{total_point[1]:0{width}x}"
    ):
        raise ValueError(f"delta_total point mismatch: {path}")
    if q_point is None:
        delta_scalar = sum(desc.delta_scalar for desc in expected) % CURVES[n].order
        if expected_total.get("scalar_relative_to_P") != delta_scalar:
            raise ValueError(f"delta_total scalar mismatch: {path}")
    elif expected_total.get("scalar_relative_to_P") is not None:
        raise ValueError(f"unknown Q discrete logarithm must not be recorded: {path}")
    return {
        "path": str(path),
        "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
        "n": n,
        "s": s,
        "tables": len(expected),
        "short_window_bits": expected[-1].window_bits,
        "table_stats_verified": verify_table_stats,
    }


def write_manifest_index(
    directory: Path,
    output: Path,
    *,
    field_sizes: set[int] | None = None,
) -> dict[str, object]:
    records = []
    for path in sorted(directory.glob("n*_s*.json")):
        record = verify_manifest(path)
        if field_sizes is not None and int(record["n"]) not in field_sizes:
            continue
        payload = json.loads(path.read_text(encoding="utf-8"))
        record["contains_table_hashes"] = all("table_stats" in row for row in payload["tables"])
        records.append(record)
    result: dict[str, object] = {
        "schema": "binary-ecc-window-manifest-index-v1",
        "public_seed_hex": PUBLIC_SEED.hex(),
        "field_sizes": sorted({int(record["n"]) for record in records}),
        "manifests": records,
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return result
