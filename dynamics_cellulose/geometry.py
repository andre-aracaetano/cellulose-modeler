from __future__ import annotations

import math

from .model import Vec3


def add(a: Vec3, b: Vec3) -> Vec3:
    return (a[0] + b[0], a[1] + b[1], a[2] + b[2])


def sub(a: Vec3, b: Vec3) -> Vec3:
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def scale(v: Vec3, factor: float) -> Vec3:
    return (v[0] * factor, v[1] * factor, v[2] * factor)


def dot(a: Vec3, b: Vec3) -> float:
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def cross(a: Vec3, b: Vec3) -> Vec3:
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def norm(v: Vec3) -> float:
    return math.sqrt(dot(v, v))


def normalize(v: Vec3) -> Vec3:
    length = norm(v)
    if length == 0:
        raise ValueError("cannot normalize a zero-length vector")
    return scale(v, 1.0 / length)


def distance(a: Vec3, b: Vec3) -> float:
    return norm(sub(a, b))


def frame(origin: Vec3, axis1: Vec3, plane: Vec3) -> tuple[Vec3, Vec3, Vec3]:
    e1 = normalize(sub(axis1, origin))
    projected = sub(sub(plane, origin), scale(e1, dot(sub(plane, origin), e1)))
    e2 = normalize(projected)
    return e1, e2, normalize(cross(e1, e2))


def to_local(point: Vec3, origin: Vec3, axes: tuple[Vec3, Vec3, Vec3]) -> Vec3:
    relative = sub(point, origin)
    return tuple(dot(relative, axis) for axis in axes)  # type: ignore[return-value]


def to_global(local: Vec3, origin: Vec3, axes: tuple[Vec3, Vec3, Vec3]) -> Vec3:
    return add(
        origin,
        add(scale(axes[0], local[0]), add(scale(axes[1], local[1]), scale(axes[2], local[2]))),
    )
