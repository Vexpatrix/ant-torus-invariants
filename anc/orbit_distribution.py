"""Enumerate the reduced blank-class orbits on 4x4, 4x6, and 6x4.

Run with Python 3.10 or later alongside count_blank_class.py and verify_ant.py.
"""

from collections import Counter, defaultdict
from fractions import Fraction
from itertools import product
from math import expm1

from count_blank_class import blank_margins, count_blank_class
from verify_ant import DIRS, is_symmetric, step, symmetries


def row_blocks(residues, rows_by_residue):
    """Generate half-grids, keyed by their column sums modulo four."""
    for rows in product(*(rows_by_residue[residue] for residue in residues)):
        yield tuple(sum(column) % 4 for column in zip(*rows)), rows


def fibre_states(m, n):
    """Enumerate the original north-facing blank label before quotienting."""
    rows_by_residue = [[] for _ in range(4)]
    for row in product((0, 1), repeat=n):
        rows_by_residue[sum(row) % 4].append(row)
    for x, y, heading, rows, columns in blank_margins(m, n):
        upper = defaultdict(list)
        for sums, block in row_blocks(rows[:m // 2], rows_by_residue):
            upper[sums].append(block)
        for sums, lower in row_blocks(rows[m // 2:], rows_by_residue):
            complement = tuple((target - partial) % 4 for target, partial in zip(columns, sums))
            for block in upper.get(complement, ()):
                grid = tuple(bit ^ ((i + j) % 2)
                             for j, row in enumerate(block + lower)
                             for i, bit in enumerate(row))
                yield grid, x, y, heading


def canonicalizer(m, n):
    """Restrict the minimization to symmetries sending the ant to a fixed pose."""
    group = symmetries(m, n)
    empty = (0,) * (m * n)
    choices = {}
    for x, y, heading in product(range(n), range(m), range(4)):
        target = (0, 0, 0 if m == n else heading % 2)
        choices[x, y, heading] = [
            symmetry for symmetry in group
            if symmetry((empty, x, y, heading))[1:] == target
        ]
        assert len(choices[x, y, heading]) == 2
    return lambda state: min(symmetry(state) for symmetry in choices[state[1:]])


def covers_width_and_back(m, n, start):
    """Check a 4n-step walk on the integer lattice before reducing by symmetries."""
    grid, x, y, heading = start
    grid = list(grid)
    lift_x = [0]
    lift_y = [0]
    for _ in range(4 * n):
        x, y, heading = step(m, n, grid, x, y, heading)
        dx, dy = DIRS[heading]
        lift_x.append(lift_x[-1] + dx)
        lift_y.append(lift_y[-1] + dy)
    return ((tuple(grid), x, y, heading) == start
            and lift_x[-1] == lift_y[-1] == 0
            and max(lift_x) - min(lift_x) == n
            and max(lift_y) - min(lift_y) == 1)


def distribution(m, n):
    canonical = canonicalizer(m, n)
    raw_fibre = set(fibre_states(m, n))
    reduced = {canonical(state) for state in raw_fibre}
    fixed = {state for state in reduced if is_symmetric(m, n, state)}
    volume, half_turn, reflection = count_blank_class(m, n)
    assert len(raw_fibre) == volume
    assert Fraction(len(reduced), len(fixed)) == Fraction(volume, half_turn + reflection)

    unseen = reduced.copy()
    symmetric_periods = Counter()
    other_periods = Counter()
    width_walks = 0
    while unseen:
        start = state = unseen.pop()
        length = 0
        fixed_points = 0
        while True:
            length += 1
            fixed_points += state in fixed
            grid, x, y, heading = state
            grid = list(grid)
            pose = step(m, n, grid, x, y, heading)
            state = canonical((tuple(grid), *pose))
            if state == start:
                break
            unseen.remove(state)
        assert fixed_points in (0, 2)
        periods = symmetric_periods if fixed_points else other_periods
        periods[length] += 1
        if fixed_points and length == 4 * n and covers_width_and_back(m, n, start):
            width_walks += 1

    symmetric_volume = sum(length * count for length, count in symmetric_periods.items())
    assert 2 * sum(symmetric_periods.values()) == len(fixed)
    phi = Fraction(symmetric_volume, len(reduced))
    scale = Fraction(2 * len(reduced), len(fixed))
    mean = Fraction(symmetric_volume, sum(symmetric_periods.values()))
    assert mean == phi * scale
    return len(reduced), len(fixed), phi, scale, mean, symmetric_periods, other_periods, width_walks


def exponential_gap(periods, scale):
    """Maximum CDF gap from Exp(z), checked on both sides of every jump."""
    total = sum(periods.values())
    previous = gap = 0
    for length, count in sorted(periods.items()):
        expected = -expm1(-float(length / scale))
        following = previous + count
        gap = max(gap, abs(previous / total - expected), abs(following / total - expected))
        previous = following
    return gap


def main():
    for m, n in ((4, 4), (4, 6), (6, 4)):
        size, fixed, phi, scale, mean, symmetric, other, width_walks = distribution(m, n)
        print(f"{m}x{n}: reduced={size}, fixed={fixed}, phi={phi}, z={scale}, mean={mean}, "
              f"Exp(z) CDF gap={exponential_gap(symmetric, scale):.6f}")
        print("  Symmetric cycles (period: count):", dict(sorted(symmetric.items())))
        print("  Other cycles (period: count):", dict(sorted(other.items())))
        print("  Symmetric 4n-step cycles covering one width and back in two rows:", width_walks)


if __name__ == "__main__":
    main()
