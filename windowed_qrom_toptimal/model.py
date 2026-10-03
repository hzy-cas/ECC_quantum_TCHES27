"""Explicit-barrier Window QROM with retained/min-TD and Balanced/min-TDW.

Complete arithmetic costs come only from gate-generated records. QROM costs
use the table-independent h_u=2n model.  Full network depths are module sums,
not a monolithic gate-stream optimum. Retained conservatively keeps square
workspace; Balanced cleans arithmetic targets/squares but retains copied data.
"""
from __future__ import annotations

MODEL_VERSION = "window-retained-lane-first-v2"
# Point-add layers do not depend on how scalar bits are grouped into windows.
LAYER_VERSION = "window-retained-clean-zv-v1"
BALANCED_MODEL_VERSION = "window-balanced-lane-first-tdw-v1"
BALANCED_LAYER_VERSION = "window-balanced-clean-zv-v1"
FIELDS = (163, 233, 283, 571)
PARAMETERS = {163: (906, 9), 233: (1341, 10), 283: (1668, 11), 571: (3569, 13)}


def architecture(objective="td"):
    if objective not in ("td", "tdw"):
        raise ValueError("objective must be td or tdw")
    return "retained" if objective == "td" else "balanced"


def model_version(objective="td"):
    return MODEL_VERSION if architecture(objective) == "retained" else BALANCED_MODEL_VERSION


def layer_version(objective="td"):
    return LAYER_VERSION if architecture(objective) == "retained" else BALANCED_LAYER_VERSION


def ceil_log2(x):
    if x < 1:
        raise ValueError("positive value required")
    return (x - 1).bit_length()


def lane_lengths(n, p):
    """Paper lane j contains S_j, S_(p+j), ...; p is the paper's w."""
    if n not in FIELDS or not 1 <= p <= 2*n+2:
        raise ValueError("unsupported n or p outside 1..2n+2")
    return [(2*n+1-j)//p+1 for j in range(p)]


def window_indices(n, s, p, lane, window):
    """Exact source indices (q*s+t)*p+j of one table in the manuscript."""
    sizes = lane_lengths(n, p)
    if not 2 <= s <= 18 or not 0 <= lane < p or window < 0:
        raise ValueError("invalid window")
    bits = min(s, max(0, sizes[lane]-window*s))
    return tuple((window*s+t)*p+lane for t in range(bits))


def window_batches(n, s, p):
    """QROM address widths by round, preserving the original lane order.

    Split the N original bit-controlled addends into p strided lanes FIRST.
    Then split each lane into windows of at most s bits. A table may contain
    both P and Q addends; no artificial boundary at the middle of N is used.
    Only active lanes appear in each tail round (a prefix of the lane list).
    """
    sizes = lane_lengths(n, p)
    if not 2 <= s <= 18:
        raise ValueError("require 2 <= s <= 18")
    return [tuple(min(s, size-q*s) for size in sizes if size > q*s)
            for q in range((sizes[0]+s-1)//s)]


def tables(n, s, p):
    return [bits for batch in window_batches(n, s, p) for bits in batch]


def schedule(n, s, p):
    batches = window_batches(n, s, p)
    result = [("initialization", batches[0])]
    result.extend(("tail", batch) for batch in batches[1:])
    items = p
    while items > 1:
        result.append(("reduction", items // 2))
        items = (items + 1) // 2
    return result


def arithmetic_depth(n, k, objective="td"):
    factor = 1 if architecture(objective) == "retained" else 2
    return factor * (ceil_log2(n - 1) + 2 * ceil_log2(k) + 2)


def arithmetic_count(n, k, objective="td"):
    m, inversion_blocks = PARAMETERS[n]
    factor = 1 if architecture(objective) == "retained" else 2
    return factor * m * (inversion_blocks + 5 * k - 3)


def liveness(n, k, kind, objective="td"):
    if kind not in ("tail", "reduction") or k < 1:
        raise ValueError("invalid layer")
    m, inversion_blocks = PARAMETERS[n]
    targets = m * (inversion_blocks + 5 * k - 3)
    clean = 4 * (m - n) * k + n * (k // 2)
    # Retain the entire square workspace, including any unused zero lines.
    retained = targets + 4 * n + (2 * n * k if kind == "reduction" else 0)
    if architecture(objective) == "balanced":
        # Montgomery, lambda, and y multiplications compute-copy-uncompute.
        # Their target pool and square workspace return to zero. The copied
        # inverse, lambda, and y output DO NOT: keep all three n*k blocks live.
        clean += m * (inversion_blocks + 3 * (k - 1)) + 4 * n
        retained = (4 if kind == "reduction" else 3) * n * k
    return dict(input_qubits=4 * n * k, clean_qubits=clean,
                retained_new_qubits=retained, qubits=4 * n * k + clean + retained,
                output_qubits=2 * n * k,
                output_alias_input_qubits=n * k if kind == "tail" else 0)


def qrom(n, b, pair=False):
    entries = 1 << b
    h = 2 * n
    lookup_t = 2 * (entries - 2)
    lookup_c = entries - 2 + 3 * h * entries - 2 * entries
    lookup_d = 3 * entries - 4 + entries * (2 * ceil_log2(h) + 1)
    if not pair:
        return dict(toffoli=lookup_t, cnot=lookup_c, toffoli_depth=lookup_t,
                    full_depth=lookup_d, scratch_qubits=2 * n + b - 2)
    pair_t = 4 * entries - 2 * b - 6
    return dict(toffoli=pair_t, cnot=2 * lookup_c + 4 * n - 2 * (h - 1),
                toffoli_depth=pair_t, full_depth=2 * lookup_d - 2 * b - 2 * ceil_log2(h) + 3,
                scratch_qubits=4 * n + b - 2)


def required_layers(n, s, p):
    return {(n, kind, len(payload) if kind == "tail" else payload)
            for kind, payload in schedule(n, s, p) if kind != "initialization"}


def predicted_depth(n, s, p, objective="td"):
    depth = 0
    for kind, payload in schedule(n, s, p):
        if kind == "initialization":
            depth += max(qrom(n, b)["toffoli_depth"] for b in payload)
        elif kind == "tail":
            depth += max(qrom(n, b, True)["toffoli_depth"] for b in payload)
            depth += arithmetic_depth(n, len(payload), objective)
        else:
            depth += arithmetic_depth(n, payload, objective)
    return 2 * depth


def predicted_width(n, s, p, objective="td"):
    """Exact peak under this allocator policy, not an optimized width bound."""
    live = peak = 4 * n + 2
    for kind, payload in schedule(n, s, p):
        if kind == "initialization":
            live += 2 * n * len(payload)
            peak = max(peak, live + sum(qrom(n, b)["scratch_qubits"] for b in payload))
            continue
        k = len(payload) if kind == "tail" else payload
        if kind == "tail":
            live += 2 * n * k
            peak = max(peak, live + sum(qrom(n, b, True)["scratch_qubits"] for b in payload))
        layer = liveness(n, k, kind, objective)
        live += layer["retained_new_qubits"]
        peak = max(peak, live + layer["clean_qubits"])
    return peak


def grid(n, s_min=2, s_max=18):
    if not 2 <= s_min <= s_max <= 18:
        raise ValueError("require 2 <= s_min <= s_max <= 18")
    for s in range(s_min, s_max + 1):
        for p in range(1, 2*n+3):
            yield s, p


def formula_minimizers(n, s_min=2, s_max=18, objective="td"):
    architecture(objective)
    scored = []
    for s, p in grid(n, s_min, s_max):
        depth = predicted_depth(n, s, p, objective)
        width = predicted_width(n, s, p, objective) if objective == "tdw" else None
        score = depth if objective == "td" else depth * width
        scored.append((score, depth, width, s, p))
    best = min(row[0] for row in scored)
    result = []
    for score, depth, width, s, p in scored:
        if score != best:
            continue
        row = dict(n=n, s=s, p=p, w=p, predicted_toffoli_depth=depth)
        if objective == "tdw":
            row.update(predicted_width=width, predicted_tdw=score)
        result.append(row)
    return result


def validate_record(record, n, kind, k, objective="td"):
    if (record.get("model_version"), record.get("n"), record.get("kind"), record.get("lanes")) != (layer_version(objective), n, kind, k):
        raise ValueError("foreign or mismatched layer record")
    for field, value in liveness(n, k, kind, objective).items():
        if record.get(field) != value:
            raise ValueError(f"liveness mismatch {n}/{kind}/{k}/{field}")
    if record.get("toffoli") != arithmetic_count(n, k, objective):
        raise ValueError("gate stream does not match selected arithmetic Toffoli count formula")
    if record.get("toffoli_depth") != arithmetic_depth(n, k, objective):
        raise ValueError("gate stream does not match selected arithmetic depth formula; do not certify the formula scan")
    for key in ("cnot", "full_depth", "current_depth"):
        if type(record.get(key)) is not int or record[key] < 0:
            raise ValueError(f"invalid {key}")


def evaluate(n, s, p, records, objective="td"):
    missing = sorted(required_layers(n, s, p) - records.keys())
    stages = schedule(n, s, p)
    row = dict(n=n, s=s, p=p, w=p, original_addends=2*n+2,
               addends=sum(len(payload) for kind, payload in stages if kind != "reduction"),
               accumulation_layers=sum(kind == "tail" for kind, _ in stages),
               tree_layers=sum(kind == "reduction" for kind, _ in stages),
               windowing_model="lane_first", objective=objective, point_addition=architecture(objective),
               predicted_toffoli_depth=predicted_depth(n, s, p, objective), available=not missing,
               missing_layers=";".join(f"{kind}:{k}" for _, kind, k in missing))
    row["predicted_width"] = predicted_width(n, s, p, objective)
    row["predicted_tdw"] = row["predicted_toffoli_depth"] * row["predicted_width"]
    if missing:
        return row
    # Reserve address and final output-copy registers for the entire network.
    live = 2 * (n + 1) + 2 * n
    peak = live
    t = c = td = fd = 0
    for kind, payload in stages:
        if kind == "initialization":
            parts = [qrom(n, b) for b in payload]
            live += 2 * n * len(payload)
            peak = max(peak, live + sum(x["scratch_qubits"] for x in parts))
            t += sum(x["toffoli"] for x in parts)
            c += sum(x["cnot"] for x in parts)
            td += max(x["toffoli_depth"] for x in parts)
            fd += max(x["full_depth"] for x in parts)
            continue
        k = len(payload) if kind == "tail" else payload
        if kind == "tail":
            parts = [qrom(n, b, True) for b in payload]
            live += 2 * n * k  # fresh z,v, retained until the outer inverse
            peak = max(peak, live + sum(x["scratch_qubits"] for x in parts))
            t += sum(x["toffoli"] for x in parts)
            c += sum(x["cnot"] for x in parts)
            td += max(x["toffoli_depth"] for x in parts)
            fd += max(x["full_depth"] for x in parts)
        layer = records[n, kind, k]
        validate_record(layer, n, kind, k, objective)
        live += layer["retained_new_qubits"]
        peak = max(peak, live + layer["clean_qubits"])
        t += layer["toffoli"]; c += layer["cnot"]
        td += layer["toffoli_depth"]; fd += layer["full_depth"]
    assert 2 * td == row["predicted_toffoli_depth"]
    assert peak == row["predicted_width"]
    row.update(toffoli=2*t, cnot=2*c+2*n, width=peak, toffoli_depth=2*td,
               nct_depth_barrier=2*fd+2, dw=(2*fd+2)*peak, tdw=2*td*peak,
               forward_retained_qubits=live, qrom_model="h_u=2n",
               depth_model="module_barriers", width_model=("retain_all_square_workspace" if objective == "td"
                   else "balanced_clean_arithmetic_retain_copied_outputs"))
    return row
