"""Single controlled affine in-place point addition from [JSB+25] Algorithm 1."""

from __future__ import annotations

import os

try:
    from .circuit_library import CNOT_n, Toffoli_gate, mul, mul_matrix
    from .config import BINARY_ECC, CONFIGS, INVERSION_MODES, OPTIMAL_DEPTH
    from .inversion import Inverison_Itoh_Tsujii_based
    from .optimal_inversion import OPTIMAL_INVERSION_FUNCTIONS
except ImportError:  # Support direct execution of a script in this directory.
    from circuit_library import CNOT_n, Toffoli_gate, mul, mul_matrix
    from config import BINARY_ECC, CONFIGS, INVERSION_MODES, OPTIMAL_DEPTH
    from inversion import Inverison_Itoh_Tsujii_based
    from optimal_inversion import OPTIMAL_INVERSION_FUNCTIONS


def Point_addition(
    gm,
    matrix_loader,
    q,
    x,
    y,
    lambda_register,
    sqr1,
    sqr2,
    ancilla,
    toffoli_qubits,
    n,
    curve_a,
    x2,
    y2,
    path_config,
    base_path,
    *,
    inversion_mode=BINARY_ECC,
    sqr3=None,
    sqr4=None,
    side_ancilla=None,
):
    """Emit Algorithm 1: ``|q>|P>|0> -> |q>|P+qP2>|0>``.

    ``curve_a``, ``x2``, and ``y2`` must already be converted to the L-basis used
    by the AC multiplier. ``inversion_mode='binary_ecc'`` selects the low-width
    serial inversion. ``inversion_mode='optimal_depth'`` selects the copied
    minimum-Toffoli-depth parallel inversion and additionally requires ``sqr3``,
    ``sqr4``, and ``side_ancilla``.

    Both divisions use compute-copy-uncompute, so every workspace register is
    restored to zero before this function returns.
    """

    if n not in CONFIGS:
        raise ValueError(f"unsupported field degree n={n}; choose from {tuple(CONFIGS)}")
    if inversion_mode not in INVERSION_MODES:
        raise ValueError(
            f"unsupported inversion mode {inversion_mode!r}; choose from {INVERSION_MODES}"
        )
    if not (
        len(x)
        == len(y)
        == len(lambda_register)
        == len(sqr1)
        == len(sqr2)
        == n
    ):
        raise ValueError("x, y, lambda, sqr1, and sqr2 must all be n-bit registers")
    cfg = CONFIGS[n]
    block_size = cfg["block_size"]
    mul_size = cfg["multiplication_targets"]
    inversion_multiplications = cfg["inversion_multiplications"]
    if len(ancilla) < 2 * block_size:
        raise ValueError("the AC multiplication ancilla register is too short")
    if len(toffoli_qubits) < (inversion_multiplications + 1) * mul_size:
        raise ValueError("the Toffoli-target register is too short")
    if len(ancilla) < n - 1:
        raise ValueError("the control fan-out ancilla register is too short")
    if inversion_mode == OPTIMAL_DEPTH:
        if sqr3 is None or sqr4 is None or side_ancilla is None:
            raise ValueError(
                "optimal_depth mode requires sqr3, sqr4, and side_ancilla"
            )
        if len(sqr3) != n or len(sqr4) != n:
            raise ValueError("all four optimal_depth squaring registers must be n bits")
        if len(side_ancilla) < 2 * block_size:
            raise ValueError("the optimal_depth side-branch multiplication workspace is too short")

    # Step 1: x <- x + x2.
    for bit in range(n):
        if (x2 >> bit) & 1:
            gm.add_X(x[bit])

    # Steps 2--3: y <- y + q*y2; then reverse the binary fan-out exactly.
    fanout_start = gm.current_pointer()
    controls = [q]
    ancilla_index = 0
    while len(controls) < n:
        current_sources = controls[:]
        for source in current_sources:
            if len(controls) >= n:
                break
            target = ancilla[ancilla_index]
            ancilla_index += 1
            gm.add_CNOT(source, target)
            controls.append(target)
    fanout_end = gm.current_pointer()
    for bit in range(n):
        if (y2 >> bit) & 1:
            gm.add_CNOT(controls[bit], y[bit])
    gm.replay_reverse(fanout_start, fanout_end)

    # Step 4: lambda <- lambda + y/x.
    division_start = gm.current_pointer()
    if inversion_mode == BINARY_ECC:
        inverse, _, next_count = Inverison_Itoh_Tsujii_based(
            gm,
            matrix_loader,
            list(x),
            n,
            list(sqr1),
            list(sqr2),
            0,
            list(ancilla),
            list(toffoli_qubits),
            path_config,
            base_path,
        )
    else:
        inverse = OPTIMAL_INVERSION_FUNCTIONS[n](
            gm,
            matrix_loader,
            list(x),
            n,
            list(sqr1),
            list(sqr2),
            list(sqr3),
            list(sqr4),
            0,
            list(ancilla),
            list(side_ancilla),
            list(toffoli_qubits),
            path_config,
            base_path,
        )
        next_count = inversion_multiplications * mul_size
    quotient_targets = toffoli_qubits[next_count : next_count + mul_size]
    quotient = mul(
        gm,
        matrix_loader,
        list(y) + list(ancilla[:block_size]),
        list(inverse) + list(ancilla[block_size : 2 * block_size]),
        quotient_targets,
        n,
        path_config,
    )
    division_end = gm.current_pointer()
    CNOT_n(gm, quotient, lambda_register, n)
    gm.replay_reverse(division_start, division_end)

    # Step 5: y <- y + x*lambda.
    multiplication_start = gm.current_pointer()
    product = mul(
        gm,
        matrix_loader,
        list(x) + list(ancilla[:block_size]),
        list(lambda_register) + list(ancilla[block_size : 2 * block_size]),
        toffoli_qubits[:mul_size],
        n,
        path_config,
    )
    multiplication_end = gm.current_pointer()
    CNOT_n(gm, product, y, n)
    gm.replay_reverse(multiplication_start, multiplication_end)

    # Steps 6--11: f=lambda^2+lambda+a+x2; x <- x+q*f; uncompute f.
    affine_start = gm.current_pointer()
    square_plus_path = os.path.join(
        base_path,
        f"square_{n}",
        "square_Matrix_Squaring_plus.txt",
    )
    affine_registers = mul_matrix(
        gm,
        matrix_loader,
        list(lambda_register) + list(y),
        2 * n,
        square_plus_path,
    )
    if len(affine_registers) != 2 * n:
        raise RuntimeError("the squaring-plus-input matrix did not return 2n logical wires")
    transformed_y = affine_registers[n:]
    affine_constant = curve_a ^ x2
    for bit in range(n):
        if (affine_constant >> bit) & 1:
            gm.add_X(transformed_y[bit])
    affine_end = gm.current_pointer()

    fanout_start = gm.current_pointer()
    controls = [q]
    ancilla_index = 0
    while len(controls) < n:
        current_sources = controls[:]
        for source in current_sources:
            if len(controls) >= n:
                break
            target = ancilla[ancilla_index]
            ancilla_index += 1
            gm.add_CNOT(source, target)
            controls.append(target)
    fanout_end = gm.current_pointer()
    for bit in range(n):
        Toffoli_gate(gm, controls[bit], transformed_y[bit], x[bit])
    gm.replay_reverse(fanout_start, fanout_end)
    gm.replay_reverse(affine_start, affine_end)

    # Step 12: y <- y + x*lambda.
    multiplication_start = gm.current_pointer()
    product = mul(
        gm,
        matrix_loader,
        list(x) + list(ancilla[:block_size]),
        list(lambda_register) + list(ancilla[block_size : 2 * block_size]),
        toffoli_qubits[:mul_size],
        n,
        path_config,
    )
    multiplication_end = gm.current_pointer()
    CNOT_n(gm, product, y, n)
    gm.replay_reverse(multiplication_start, multiplication_end)

    # Step 13: lambda <- lambda + y/x, thereby clearing lambda.
    division_start = gm.current_pointer()
    if inversion_mode == BINARY_ECC:
        inverse, _, next_count = Inverison_Itoh_Tsujii_based(
            gm,
            matrix_loader,
            list(x),
            n,
            list(sqr1),
            list(sqr2),
            0,
            list(ancilla),
            list(toffoli_qubits),
            path_config,
            base_path,
        )
    else:
        inverse = OPTIMAL_INVERSION_FUNCTIONS[n](
            gm,
            matrix_loader,
            list(x),
            n,
            list(sqr1),
            list(sqr2),
            list(sqr3),
            list(sqr4),
            0,
            list(ancilla),
            list(side_ancilla),
            list(toffoli_qubits),
            path_config,
            base_path,
        )
        next_count = inversion_multiplications * mul_size
    quotient_targets = toffoli_qubits[next_count : next_count + mul_size]
    quotient = mul(
        gm,
        matrix_loader,
        list(y) + list(ancilla[:block_size]),
        list(inverse) + list(ancilla[block_size : 2 * block_size]),
        quotient_targets,
        n,
        path_config,
    )
    division_end = gm.current_pointer()
    CNOT_n(gm, quotient, lambda_register, n)
    gm.replay_reverse(division_start, division_end)

    # Steps 14--17: x <- x+x2; y <- y+q*y2+q*x.
    for bit in range(n):
        if (x2 >> bit) & 1:
            gm.add_X(x[bit])

    fanout_start = gm.current_pointer()
    controls = [q]
    ancilla_index = 0
    while len(controls) < n:
        current_sources = controls[:]
        for source in current_sources:
            if len(controls) >= n:
                break
            target = ancilla[ancilla_index]
            ancilla_index += 1
            gm.add_CNOT(source, target)
            controls.append(target)
    fanout_end = gm.current_pointer()
    for bit in range(n):
        if (y2 >> bit) & 1:
            gm.add_CNOT(controls[bit], y[bit])
    for bit in range(n):
        Toffoli_gate(gm, controls[bit], x[bit], y[bit])
    gm.replay_reverse(fanout_start, fanout_end)

    return x, y, lambda_register
