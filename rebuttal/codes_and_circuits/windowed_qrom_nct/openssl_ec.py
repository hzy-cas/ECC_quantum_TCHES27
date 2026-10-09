"""Small ctypes binding for deterministic table generation with OpenSSL 1.1.

The workspace image has the OpenSSL runtime but not its development headers.
Using the stable 1.1 C ABI keeps point generation in optimized C while avoiding
an undeclared Python elliptic-curve dependency.
"""

from __future__ import annotations

import ctypes
import ctypes.util
from collections.abc import Iterator
from contextlib import AbstractContextManager

from .curves import BinaryCurve


class OpenSSLError(RuntimeError):
    pass


def _load_crypto() -> ctypes.CDLL:
    candidates = ("libcrypto.so.1.1", ctypes.util.find_library("crypto"), "libcrypto.so")
    for candidate in candidates:
        if not candidate:
            continue
        try:
            return ctypes.CDLL(candidate)
        except OSError:
            pass
    raise OpenSSLError("could not load libcrypto")


_LIB = _load_crypto()
_P = ctypes.c_void_p


def _bind(name: str, result, *arguments):
    fn = getattr(_LIB, name)
    fn.restype = result
    fn.argtypes = list(arguments)
    return fn


OBJ_txt2nid = _bind("OBJ_txt2nid", ctypes.c_int, ctypes.c_char_p)
EC_GROUP_new_by_curve_name = _bind("EC_GROUP_new_by_curve_name", _P, ctypes.c_int)
EC_GROUP_free = _bind("EC_GROUP_free", None, _P)
EC_GROUP_get0_generator = _bind("EC_GROUP_get0_generator", _P, _P)
EC_POINT_new = _bind("EC_POINT_new", _P, _P)
EC_POINT_free = _bind("EC_POINT_free", None, _P)
EC_POINT_copy = _bind("EC_POINT_copy", ctypes.c_int, _P, _P)
EC_POINT_add = _bind("EC_POINT_add", ctypes.c_int, _P, _P, _P, _P, _P)
EC_POINT_mul = _bind("EC_POINT_mul", ctypes.c_int, _P, _P, _P, _P, _P, _P)
EC_POINT_is_at_infinity = _bind("EC_POINT_is_at_infinity", ctypes.c_int, _P, _P)
EC_POINTs_make_affine = _bind(
    "EC_POINTs_make_affine", ctypes.c_int, _P, ctypes.c_size_t, ctypes.POINTER(_P), _P
)
EC_POINT_get_affine_coordinates_GF2m = _bind(
    "EC_POINT_get_affine_coordinates_GF2m", ctypes.c_int, _P, _P, _P, _P, _P
)
EC_POINT_set_affine_coordinates_GF2m = _bind(
    "EC_POINT_set_affine_coordinates_GF2m", ctypes.c_int, _P, _P, _P, _P, _P
)
EC_POINT_is_on_curve = _bind("EC_POINT_is_on_curve", ctypes.c_int, _P, _P, _P)
BN_CTX_new = _bind("BN_CTX_new", _P)
BN_CTX_free = _bind("BN_CTX_free", None, _P)
BN_new = _bind("BN_new", _P)
BN_free = _bind("BN_free", None, _P)
BN_bin2bn = _bind("BN_bin2bn", _P, ctypes.c_void_p, ctypes.c_int, _P)
BN_bn2binpad = _bind("BN_bn2binpad", ctypes.c_int, _P, ctypes.c_void_p, ctypes.c_int)


def _check(ok: int, operation: str) -> None:
    if ok != 1:
        raise OpenSSLError(f"OpenSSL operation failed: {operation}")


class PointGenerator(AbstractContextManager["PointGenerator"]):
    """Generate affine points ``[(start + i*step)G]`` in bounded chunks."""

    def __init__(self, curve: BinaryCurve):
        self.curve = curve
        nid = OBJ_txt2nid(curve.openssl_name.encode("ascii"))
        if nid == 0:
            raise OpenSSLError(f"OpenSSL does not know {curve.openssl_name}")
        self.group = EC_GROUP_new_by_curve_name(nid)
        self.ctx = BN_CTX_new()
        if not self.group or not self.ctx:
            self.close()
            raise OpenSSLError("could not allocate EC group/BN context")
        self._closed = False

    def close(self) -> None:
        if getattr(self, "_closed", True):
            return
        if self.ctx:
            BN_CTX_free(self.ctx)
            self.ctx = None
        if self.group:
            EC_GROUP_free(self.group)
            self.group = None
        self._closed = True

    def __exit__(self, exc_type, exc, traceback) -> None:
        self.close()

    def _scalar_bn(self, scalar: int):
        width = max(1, (scalar.bit_length() + 7) // 8)
        raw = scalar.to_bytes(width, "big")
        buffer = ctypes.create_string_buffer(raw)
        result = BN_bin2bn(buffer, len(raw), None)
        if not result:
            raise OpenSSLError("BN_bin2bn failed")
        return result

    def scalar_point(self, scalar: int):
        point = EC_POINT_new(self.group)
        bn = self._scalar_bn(scalar % self.curve.order)
        if not point:
            BN_free(bn)
            raise OpenSSLError("EC_POINT_new failed")
        try:
            _check(EC_POINT_mul(self.group, point, bn, None, None, self.ctx), "EC_POINT_mul")
        finally:
            BN_free(bn)
        return point

    def affine_point(self, x_coordinate: int, y_coordinate: int):
        """Allocate and validate a finite affine point on this curve."""
        if min(x_coordinate, y_coordinate) < 0:
            raise ValueError("affine coordinates must be nonnegative")
        if x_coordinate.bit_length() > self.curve.n or y_coordinate.bit_length() > self.curve.n:
            raise ValueError("affine coordinate exceeds the field width")
        point = EC_POINT_new(self.group)
        x = self._scalar_bn(x_coordinate)
        y = self._scalar_bn(y_coordinate)
        if not point:
            BN_free(x)
            BN_free(y)
            raise OpenSSLError("EC_POINT_new failed")
        try:
            _check(
                EC_POINT_set_affine_coordinates_GF2m(
                    self.group, point, x, y, self.ctx
                ),
                "EC_POINT_set_affine_coordinates_GF2m",
            )
            if EC_POINT_is_at_infinity(self.group, point) == 1:
                raise ValueError("the public point must be finite")
            if EC_POINT_is_on_curve(self.group, point, self.ctx) != 1:
                raise ValueError("the public point is not on the selected curve")
            return point
        except Exception:
            EC_POINT_free(point)
            raise
        finally:
            BN_free(x)
            BN_free(y)

    def multiply_point(self, point, scalar: int):
        result = EC_POINT_new(self.group)
        bn = self._scalar_bn(scalar % self.curve.order)
        if not result:
            BN_free(bn)
            raise OpenSSLError("EC_POINT_new failed")
        try:
            _check(
                EC_POINT_mul(self.group, result, None, point, bn, self.ctx),
                "EC_POINT_mul with an affine base",
            )
        except Exception:
            EC_POINT_free(result)
            raise
        finally:
            BN_free(bn)
        return result

    def validate_prime_subgroup_point(self, x_coordinate: int, y_coordinate: int) -> None:
        """Require a finite on-curve point in the subgroup of known order."""
        point = self.affine_point(x_coordinate, y_coordinate)
        multiple = None
        try:
            multiple = self.multiply_point(point, self.curve.order)
            if EC_POINT_is_at_infinity(self.group, multiple) != 1:
                raise ValueError("the public point is outside the prime-order subgroup")
        finally:
            if multiple:
                EC_POINT_free(multiple)
            EC_POINT_free(point)

    def multiple_coordinates(
        self, x_coordinate: int, y_coordinate: int, scalar: int
    ) -> tuple[int, int] | None:
        """Return ``[scalar](x,y)``, using ``None`` for the point at infinity."""
        point = self.affine_point(x_coordinate, y_coordinate)
        result = None
        try:
            result = self.multiply_point(point, scalar)
            if EC_POINT_is_at_infinity(self.group, result) == 1:
                return None
            return self.coordinates(result)
        finally:
            if result:
                EC_POINT_free(result)
            EC_POINT_free(point)

    def sum_coordinates(
        self, points: Iterator[tuple[int, int]]
    ) -> tuple[int, int] | None:
        """Add finite affine points, returning ``None`` for the identity."""
        accumulator = EC_POINT_new(self.group)
        if not accumulator:
            raise OpenSSLError("EC_POINT_new failed")
        initialized = False
        try:
            for x_coordinate, y_coordinate in points:
                point = self.affine_point(x_coordinate, y_coordinate)
                try:
                    if not initialized:
                        _check(EC_POINT_copy(accumulator, point), "EC_POINT_copy")
                        initialized = True
                    else:
                        _check(
                            EC_POINT_add(
                                self.group, accumulator, accumulator, point, self.ctx
                            ),
                            "EC_POINT_add",
                        )
                finally:
                    EC_POINT_free(point)
            if not initialized or EC_POINT_is_at_infinity(self.group, accumulator) == 1:
                return None
            return self.coordinates(accumulator)
        finally:
            EC_POINT_free(accumulator)

    def coordinates(self, point) -> tuple[int, int]:
        if EC_POINT_is_at_infinity(self.group, point) == 1:
            raise OpenSSLError("cannot encode point at infinity as affine coordinates")
        x = BN_new()
        y = BN_new()
        if not x or not y:
            if x:
                BN_free(x)
            if y:
                BN_free(y)
            raise OpenSSLError("BN_new failed")
        width = self.curve.coordinate_bytes
        xb = ctypes.create_string_buffer(width)
        yb = ctypes.create_string_buffer(width)
        try:
            _check(
                EC_POINT_get_affine_coordinates_GF2m(self.group, point, x, y, self.ctx),
                "EC_POINT_get_affine_coordinates_GF2m",
            )
            if BN_bn2binpad(x, xb, width) != width or BN_bn2binpad(y, yb, width) != width:
                raise OpenSSLError("BN_bn2binpad failed")
            return int.from_bytes(xb.raw, "big"), int.from_bytes(yb.raw, "big")
        finally:
            BN_free(x)
            BN_free(y)

    def progression(
        self, start_scalar: int, step_scalar: int, count: int, *, chunk_size: int = 4096
    ) -> Iterator[tuple[int, int]]:
        if count < 0 or chunk_size <= 0:
            raise ValueError("count must be nonnegative and chunk_size positive")
        if count == 0:
            return
        current = self.scalar_point(start_scalar)
        step = self.scalar_point(step_scalar)
        try:
            remaining = count
            while remaining:
                size = min(chunk_size, remaining)
                points = [EC_POINT_new(self.group) for _ in range(size)]
                if any(not point for point in points):
                    for point in points:
                        if point:
                            EC_POINT_free(point)
                    raise OpenSSLError("EC_POINT_new failed in progression")
                try:
                    _check(EC_POINT_copy(points[0], current), "EC_POINT_copy")
                    for index in range(1, size):
                        _check(
                            EC_POINT_add(
                                self.group, points[index], points[index - 1], step, self.ctx
                            ),
                            "EC_POINT_add",
                        )
                    array = (_P * size)(*points)
                    _check(
                        EC_POINTs_make_affine(self.group, size, array, self.ctx),
                        "EC_POINTs_make_affine",
                    )
                    for point in points:
                        yield self.coordinates(point)
                    _check(
                        EC_POINT_add(self.group, current, points[-1], step, self.ctx),
                        "advance progression",
                    )
                finally:
                    for point in points:
                        EC_POINT_free(point)
                remaining -= size
        finally:
            EC_POINT_free(current)
            EC_POINT_free(step)

    def affine_progression(
        self,
        base_x: int,
        base_y: int,
        start_multiplier: int,
        step_multiplier: int,
        count: int,
        *,
        chunk_size: int = 4096,
    ) -> Iterator[tuple[int, int]]:
        """Generate ``[(start+i*step)B]`` for a supplied affine base ``B``."""
        if count < 0 or chunk_size <= 0:
            raise ValueError("count must be nonnegative and chunk_size positive")
        if count == 0:
            return
        base = self.affine_point(base_x, base_y)
        current = None
        step = None
        try:
            current = self.multiply_point(base, start_multiplier)
            step = self.multiply_point(base, step_multiplier)
            if EC_POINT_is_at_infinity(self.group, current) == 1:
                raise ValueError("shifted table starts at the point at infinity")
            if EC_POINT_is_at_infinity(self.group, step) == 1:
                raise ValueError("window step is the point at infinity")
            remaining = count
            while remaining:
                size = min(chunk_size, remaining)
                points = [EC_POINT_new(self.group) for _ in range(size)]
                if any(not point for point in points):
                    for point in points:
                        if point:
                            EC_POINT_free(point)
                    raise OpenSSLError("EC_POINT_new failed in affine progression")
                try:
                    _check(EC_POINT_copy(points[0], current), "EC_POINT_copy")
                    for index in range(1, size):
                        _check(
                            EC_POINT_add(
                                self.group, points[index], points[index - 1], step, self.ctx
                            ),
                            "EC_POINT_add",
                        )
                        if EC_POINT_is_at_infinity(self.group, points[index]) == 1:
                            raise ValueError("shifted table contains the point at infinity")
                    array = (_P * size)(*points)
                    _check(
                        EC_POINTs_make_affine(self.group, size, array, self.ctx),
                        "EC_POINTs_make_affine",
                    )
                    for point in points:
                        yield self.coordinates(point)
                    _check(
                        EC_POINT_add(self.group, current, points[-1], step, self.ctx),
                        "advance affine progression",
                    )
                finally:
                    for point in points:
                        EC_POINT_free(point)
                remaining -= size
        finally:
            if current:
                EC_POINT_free(current)
            if step:
                EC_POINT_free(step)
            EC_POINT_free(base)
