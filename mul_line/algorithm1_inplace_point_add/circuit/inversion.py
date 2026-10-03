"""Inversion chains from *New Quantum Cryptanalysis of Binary Elliptic Curves*."""

from __future__ import annotations

import os

try:
    from .circuit_library import mul, mul_matrix
    from .config import CONFIGS
except ImportError:  # Support direct execution of a script in this directory.
    from circuit_library import mul, mul_matrix
    from config import CONFIGS


_PACKAGE_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_OPTIMIZED_FROBENIUS_ROOT = os.path.join(_PACKAGE_ROOT, "square")
_REQUIRED_FROBENIUS_POWERS = {
    163: (1, 2, 4, 8, 16, 32, 34, 64),
    233: (1, 2, 4, 8, 16, 32, 40, 64, 104),
    283: (1, 2, 4, 8, 10, 16, 26, 32, 64, 128),
    571: (1, 2, 4, 8, 10, 16, 26, 32, 58, 64, 128, 256),
}


def Inverison_Itoh_Tsujii_based(
    gm,
    matrix_loader,
    a,
    n,
    sqr1,
    sqr2,
    count,
    ancilla,
    toffoli_qubits,
    path_config,
    base_path,
):
    """Emit the serial Itoh--Tsujii inversion for four NIST binary fields.

    The function retains the existing ``Inverison_Itoh_Tsujii_based`` spelling used
    by Binary_ECC and ``recursive_karatsuba``. Every nonlinear multiplication calls
    the local ``circuit_library.mul`` directly, and the return value preserves the
    ``(inverse, ancilla, count)`` interface. ``base_path`` remains in the established
    call signature but is no longer used to locate squaring linear-layer files.
    """

    if n not in CONFIGS:
        raise ValueError(f"unsupported field degree n={n}; choose from {tuple(CONFIGS)}")
    cfg = CONFIGS[n]
    block_size = cfg["block_size"]
    mul_size = cfg["multiplication_targets"]
    expected_multiplications = cfg["inversion_multiplications"]
    if not (len(a) == len(sqr1) == len(sqr2) == n):
        raise ValueError("the inversion input and both squaring registers must be n bits")
    if len(ancilla) < 2 * block_size:
        raise ValueError("the AC multiplication ancilla register is too short")
    if path_config["mul_dim_2"] != mul_size:
        raise ValueError("path_config has the wrong multiplication-target count for this field")

    missing_square_files = []
    for power in _REQUIRED_FROBENIUS_POWERS[n]:
        file_path = os.path.join(
            _OPTIMIZED_FROBENIUS_ROOT,
            f"RES_n{n}",
            f"result_square_Matrix_2_{power}.txt",
        )
        if not os.path.isfile(file_path):
            missing_square_files.append(file_path)
    if missing_square_files:
        raise FileNotFoundError(
            "missing optimized squaring circuits required by the inversion chain: "
            + ", ".join(missing_square_files)
            + ". The squaring linear layer has no matrix-power fallback."
        )

    start_count = count

    def square(source, target, power):
        optimized_file = os.path.join(
            _OPTIMIZED_FROBENIUS_ROOT,
            f"RES_n{n}",
            f"result_square_Matrix_2_{power}.txt",
        )
        if not os.path.isfile(optimized_file):
            raise FileNotFoundError(
                f"missing optimized squaring circuit: {optimized_file}. "
                "The squaring linear layer reads only square/RES_n*/"
                "result_square_Matrix_2_*.txt and no longer derives matrix powers "
                "from the 2^1 map."
            )

        outputs = mul_matrix(
            gm,
            matrix_loader,
            list(source) + list(target),
            2 * n,
            optimized_file,
        )
        if len(outputs) != 2 * n or len(set(outputs)) != 2 * n:
            raise ValueError(
                f"optimized squaring circuit {optimized_file} lacks a complete logical-output map"
            )
        return outputs[:n], outputs[n:]

    def multiply(left, right):
        nonlocal count
        targets = toffoli_qubits[count : count + mul_size]
        if len(targets) != mul_size:
            raise ValueError("the inversion Toffoli-target register is too short")
        product = mul(
            gm,
            matrix_loader,
            list(left) + list(ancilla[:block_size]),
            list(right) + list(ancilla[block_size : 2 * block_size]),
            targets,
            n,
            path_config,
        )
        count += mul_size
        return product

    # Common prefix: 1 -> 2 -> 4 -> 8 -> 16.
    a, a_squared = square(a, sqr1, 1)
    a_2_0 = multiply(a, a_squared)
    a_2_0, a_2_2 = square(a_2_0, sqr2, 2)
    a, sqr1 = square(a, a_squared, 1)
    b_value = multiply(a_2_2, a_2_0)
    a_2_0, sqr2 = square(a_2_0, a_2_2, 2)
    b_value, b_2_4 = square(b_value, sqr1, 4)
    c_value = multiply(b_value, b_2_4)
    b_value, sqr1 = square(b_value, b_2_4, 4)
    c_value, c_2_8 = square(c_value, sqr2, 8)
    d_value = multiply(c_value, c_2_8)
    c_value, sqr2 = square(c_value, c_2_8, 8)
    d_value, d_2_16 = square(d_value, sqr1, 16)

    if n == 163:
        e_value = multiply(d_value, d_2_16)
        d_value, sqr1 = square(d_value, d_2_16, 16)
        e_value, e_2_32 = square(e_value, sqr2, 32)
        f_value = multiply(e_value, e_2_32)
        e_value, sqr2 = square(e_value, e_2_32, 32)
        e_value, e_2_2 = square(e_value, sqr2, 2)
        f_value, f_2_64 = square(f_value, sqr1, 64)
        g_value = multiply(f_value, f_2_64)
        f_value, sqr1 = square(f_value, f_2_64, 64)
        g_value, g_2_34 = square(g_value, sqr1, 34)
        partial = multiply(a_2_0, e_2_2)
        e_value, sqr2 = square(e_value, e_2_2, 2)
        last_product = multiply(partial, g_2_34)
        last_product, inverse = square(last_product, sqr2, 1)

    elif n == 233:
        e_value = multiply(d_value, d_2_16)
        d_value, sqr1 = square(d_value, d_2_16, 16)
        e_value, e_2_32 = square(e_value, sqr2, 32)
        f_value = multiply(e_value, e_2_32)
        e_value, sqr2 = square(e_value, e_2_32, 32)
        f_value, f_2_64 = square(f_value, sqr1, 64)
        e_value, e_2_8 = square(e_value, sqr2, 8)
        g_value = multiply(f_value, f_2_64)
        f_value, sqr1 = square(f_value, f_2_64, 64)
        f_value, f_2_40 = square(f_value, sqr1, 40)
        partial_one = multiply(c_value, e_2_8)
        e_value, sqr2 = square(e_value, e_2_8, 8)
        g_value, g_2_104 = square(g_value, sqr2, 104)
        partial_two = multiply(partial_one, f_2_40)
        f_value, sqr1 = square(f_value, f_2_40, 40)
        last_product = multiply(partial_two, g_2_104)
        last_product, inverse = square(last_product, sqr1, 1)

    elif n == 283:
        c_value, c_2_2 = square(c_value, sqr2, 2)
        e_value = multiply(d_value, d_2_16)
        d_value, sqr1 = square(d_value, d_2_16, 16)
        e_value, e_2_32 = square(e_value, sqr1, 32)
        partial_one = multiply(a_2_0, c_2_2)
        c_value, sqr2 = square(c_value, c_2_2, 2)
        d_value, d_2_10 = square(d_value, sqr2, 10)
        f_value = multiply(e_value, e_2_32)
        e_value, sqr1 = square(e_value, e_2_32, 32)
        f_value, f_2_64 = square(f_value, sqr1, 64)
        partial_two = multiply(partial_one, d_2_10)
        d_value, sqr2 = square(d_value, d_2_10, 10)
        g_value = multiply(f_value, f_2_64)
        f_value, sqr1 = square(f_value, f_2_64, 64)
        g_value, g_2_128 = square(g_value, sqr2, 128)
        h_value = multiply(g_value, g_2_128)
        g_value, sqr2 = square(g_value, g_2_128, 128)
        h_value, h_2_26 = square(h_value, sqr1, 26)
        last_product = multiply(partial_two, h_2_26)
        last_product, inverse = square(last_product, sqr2, 1)

    else:  # n == 571
        c_value, c_2_2 = square(c_value, sqr2, 2)
        e_value = multiply(d_value, d_2_16)
        d_value, sqr1 = square(d_value, d_2_16, 16)
        e_value, e_2_32 = square(e_value, sqr1, 32)
        partial_one = multiply(a_2_0, c_2_2)
        c_value, sqr2 = square(c_value, c_2_2, 2)
        d_value, d_2_10 = square(d_value, sqr2, 10)
        f_value = multiply(e_value, e_2_32)
        e_value, sqr1 = square(e_value, e_2_32, 32)
        f_value, f_2_64 = square(f_value, sqr1, 64)
        partial_two = multiply(partial_one, d_2_10)
        d_value, sqr2 = square(d_value, d_2_10, 10)
        e_value, e_2_26 = square(e_value, sqr2, 26)
        g_value = multiply(f_value, f_2_64)
        f_value, sqr1 = square(f_value, f_2_64, 64)
        g_value, g_2_128 = square(g_value, sqr1, 128)
        partial_three = multiply(partial_two, e_2_26)
        e_value, sqr2 = square(e_value, e_2_26, 26)
        h_value = multiply(g_value, g_2_128)
        g_value, sqr1 = square(g_value, g_2_128, 128)
        h_value, h_2_256 = square(h_value, sqr2, 256)
        i_value = multiply(h_value, h_2_256)
        h_value, sqr2 = square(h_value, h_2_256, 256)
        i_value, i_2_58 = square(i_value, sqr1, 58)
        last_product = multiply(partial_three, i_2_58)
        last_product, inverse = square(last_product, sqr2, 1)

    expected_count = start_count + expected_multiplications * mul_size
    if count != expected_count:
        raise AssertionError(
            f"GF(2^{n}) inversion used {count - start_count} multiplication targets; "
            f"expected {expected_count - start_count}"
        )
    return inverse, ancilla, count
