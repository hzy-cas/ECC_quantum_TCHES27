"""Deterministic ordinary affine-add fixtures in the bundled AC basis."""
from .curves import CURVES, gf_inverse, gf_multiply, gf_square, affine_add, on_curve
from .basis import poly_to_l


def make_fixture(n, lanes):
    curve = CURVES[n]
    generator = (curve.gx, curve.gy)
    x, y = generator
    slope = x ^ gf_multiply(y, gf_inverse(x, curve), curve)
    x2 = gf_square(slope, curve) ^ slope ^ curve.a
    y2 = gf_square(x, curve) ^ gf_multiply(slope ^ 1, x2, curve)
    points = [generator, (x2, y2)]
    while len(points) < 2 * lanes:
        points.append(affine_add(points[-1], generator, curve))
    lines = [f"{n} {lanes}"]
    for i in range(lanes):
        left, right = points[2*i:2*i+2]
        output = affine_add(left, right, curve)
        assert on_curve(*left, curve) and on_curve(*right, curve) and on_curve(*output, curve)
        values = (*left, *right, *output)
        lines.append(" ".join(hex(poly_to_l(v, n)) for v in values))
    return "\n".join(lines) + "\n"
