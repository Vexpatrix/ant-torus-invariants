"""Count the exact blank-class scale of a rectangular Langton ant torus."""

import argparse
from collections import defaultdict
from fractions import Fraction
from functools import cache
from itertools import product

from verify_ant import DIRS, all_invariants, is_symmetric


ROOTS = ((1, 0), (0, 1), (-1, 0), (0, -1))
ONE_PLUS_ROOT = ((2, 0), (1, 1), (0, 0), (1, -1))


def add(left, right):
    return left[0] + right[0], left[1] + right[1]


def multiply(left, right):
    a, b = left
    c, d = right
    return a * c - b * d, a * d + b * c


def power(base, exponent):
    result = (1, 0)
    for _ in range(exponent):
        result = multiply(result, base)
    return result


@cache
def row_polynomial(frequencies, residue, modulus):
    """Sum the Fourier monomials for rows with the prescribed residue."""
    result = (0, 0)
    for index in range(0, 4, 4 // modulus):
        term = ROOTS[-index * residue % 4]
        for root, multiplicity in enumerate(frequencies):
            term = multiply(term, power(ONE_PLUS_ROOT[(root + index) % 4], multiplicity))
        result = add(result, term)
    assert result[0] % modulus == result[1] % modulus == 0
    return result[0] // modulus, result[1] // modulus


@cache
def count_matrices(rows, targets, mode, free_columns=0, *, row_modulus=4, fixed_column_parity=None):
    """Count rows of given residues with constraints down the columns.

    In 'sum' mode targets are single-column sums mod 4. In 'parity' mode
    they are single-column sums mod 2. In 'difference' mode each target
    is a difference of sums for a disjoint pair of columns mod 4.
    In 'pair_sum' mode each target is a pair's sum plus two free bits.
    A fixed column, when present, has its own parity target.
    """
    coefficients = {(free_columns, 0, 0, 0): (1, 0)}
    if fixed_column_parity is not None:
        coefficients = {(1, 0, 0, 0): (1, 0), (0, 0, 1, 0): ROOTS[-2 * fixed_column_parity % 4]}
    roots = (0, 2) if mode == "parity" else range(4)
    for target in targets:
        updated = defaultdict(lambda: (0, 0))
        phases = {root: ROOTS[-root * target % 4] for root in roots}
        if mode == "pair_sum":
            phases = {root: multiply(phase, power(ONE_PLUS_ROOT[root], 2))
                      for root, phase in phases.items()}
        for frequencies, coefficient in coefficients.items():
            for root, phase in phases.items():
                updated_frequencies = list(frequencies)
                updated_frequencies[root] += 1
                if mode == "difference":
                    updated_frequencies[-root % 4] += 1
                elif mode == "pair_sum":
                    updated_frequencies[root] += 1
                key = tuple(updated_frequencies)
                updated[key] = add(updated[key], multiply(coefficient, phase))
        coefficients = updated

    result = (0, 0)
    for frequencies, term in coefficients.items():
        for residue in rows:
            term = multiply(term, row_polynomial(frequencies, residue, row_modulus))
        result = add(result, term)
    denominator = len(roots) ** len(targets) * (2 if fixed_column_parity is not None else 1)
    assert result[1] == 0 and result[0] % denominator == 0 and result[0] >= 0
    return result[0] // denominator


def alternating_sums(adjacent_sums, parities):
    """Recover each of the up to two signed-margin vectors mod 4."""
    for first in range(parities[0], 4, 2):
        values = [first]
        for adjacent in adjacent_sums[:-1]:
            values.append((adjacent - values[-1]) % 4)
        if (values[-1] + first) % 4 != adjacent_sums[-1]:
            continue
        if all(value % 2 == parity for value, parity in zip(values, parities)):
            yield values


def reflected_pairs(length, offset):
    """Pairs of distinct indices under reversal with the given offset."""
    partners = [(offset - index) % length for index in range(length)]
    return [(a, b) for a, b in enumerate(partners) if a < b]


def complementary_row_half_turn(rows, columns, x, row_pairs, row_modulus):
    """Count a half-turn that complements bits between the paired rows."""
    m, n = len(rows), len(columns)
    if any((rows[a] + rows[b] - n) % row_modulus for a, b in row_pairs):
        return 0
    partners = [(2 * x - column) % n for column in range(n)]
    if any((columns[a] + columns[b] - m) % 4 for a, b in enumerate(partners)):
        return 0
    if any((columns[a] - m // 2) % 4 for a, b in enumerate(partners) if a == b):
        return 0
    differences = [(columns[a] - m // 2) % 4 for a, b in enumerate(partners) if a < b]
    targets = tuple(sorted(min(value, -value % 4) for value in differences))
    chosen_rows = tuple(sorted(rows[a] for a, _ in row_pairs))
    return count_matrices(chosen_rows, targets, "difference", n - 2 * len(targets),
                          row_modulus=row_modulus)


def symmetric_counts(rows, columns, x, y, heading):
    """Count half-turn and complement-reflection cases for fixed margins."""
    dx, dy = DIRS[heading]
    if dx:
        # Transpose the transformed grid to put the heading on the row axis.
        rows, columns, x, y, dy = columns, rows, y, x, dx
    row_pairs = reflected_pairs(len(rows), 2 * y - dy)
    independent_rows = tuple(sorted(rows[a] for a, _ in row_pairs))
    half_turn = complementary_row_half_turn(rows, columns, x, row_pairs, 4)

    reflection = 0
    if all(rows[a] == rows[b] for a, b in row_pairs) and all(value % 2 == 0 for value in columns):
        targets = tuple(sorted(value // 2 for value in columns))
        reflection = count_matrices(independent_rows, targets, "parity")
    return half_turn, reflection


def color_targets(invariants, blank, m, n, blank_heading=0):
    """Yield the required grid contributions to the blank state's full label."""
    label = {name: evaluate(blank, 0, 0, blank_heading) for name, _, evaluate in invariants}
    for x, y, heading in product(range(n), range(m), range(4)):
        delta = {name: (label[name] - evaluate(blank, x, y, heading)) % modulus
                 for name, modulus, evaluate in invariants}
        if delta.get("align", 0) == 0:
            yield x, y, heading, delta


def blank_margins(m, n):
    """Yield the poses and transformed margins in the blank fibre."""
    if m < 2 or n < 2 or m % 2 or n % 2:
        raise ValueError("both dimensions must be even and at least two")
    invariants = all_invariants(m, n)
    blank = [0] * (m * n)
    for x, y, heading, delta in color_targets(invariants, blank, m, n):
        column_sums = [delta[f"P{i}"] for i in range(n)]
        column_parities = [delta[f"Q{i}"] for i in range(n)]
        # The transposed row invariant sees the complement, reversing the signs.
        row_sums = [-delta[f"row:P{i}"] % 4 for i in range(m)]
        row_parities = [delta[f"row:Q{i}"] for i in range(m)]
        for a, b in product(alternating_sums(column_sums, column_parities),
                            alternating_sums(row_sums, row_parities)):
            columns = [(m // 2 + (-1) ** i * value) % 4 for i, value in enumerate(a)]
            rows = [(n // 2 + (-1) ** i * value) % 4 for i, value in enumerate(b)]
            if (sum(columns) - sum(rows)) % 4:
                continue
            black = m * n // 2 + sum(columns) - 2 * (sum(columns[1::2]) + sum(rows[1::2]))
            if black % 4 != delta["J"]:
                continue
            yield x, y, heading, rows, columns


def mixed_half_turn(rows, columns, x, y, heading):
    """Count half-turn states on an even-height, odd-width board."""
    m, n = len(rows), len(columns)
    dx, dy = DIRS[heading]
    row_pairs = reflected_pairs(m, 2 * y - dy)
    if dy:
        return complementary_row_half_turn(rows, columns, x, row_pairs, 2)

    if any(rows[a] != rows[b] for a, b in row_pairs):
        return 0
    column_offset = 2 * x - dx
    column_pairs = reflected_pairs(n, column_offset)
    if any(columns[a] != columns[b] for a, b in column_pairs):
        return 0
    fixed_col = (column_offset * ((n + 1) // 2)) % n
    # A fixed row's bit at the fixed column is its row parity.
    remainder = (columns[fixed_col] - rows[y] - rows[(y + m // 2) % m]) % 4
    if remainder % 2:
        return 0
    chosen_rows = tuple(sorted(rows[a] for a, _ in row_pairs))
    targets = tuple(sorted(columns[a] for a, _ in column_pairs))
    return count_matrices(chosen_rows, targets, "pair_sum", row_modulus=2,
                          fixed_column_parity=remainder // 2)


def count_mixed(m, n, blank_heading):
    """Count the mixed case after rotating to an even height."""
    invariants = all_invariants(m, n)
    blank = [0] * (m * n)
    volume = half_turn = 0
    for x, y, heading, delta in color_targets(invariants, blank, m, n, blank_heading):
        columns = [(m // 2 + delta[f"U{j}"]) % 4 for j in range(n)]
        rows = [(delta[f"row:Q{i}"] + i) % 2 for i in range(m)]
        black = m * n // 2 + sum(columns) - 2 * sum(rows[1::2])
        if (sum(rows) - sum(columns)) % 2 or black % 4 != delta["J"]:
            continue
        volume += count_matrices(tuple(sorted(rows)), tuple(sorted(columns)), "sum", row_modulus=2)
        half_turn += mixed_half_turn(rows, columns, x, y, heading)
    return volume, half_turn, 0


def count_blank_class(m, n):
    """Count volume and both types of symmetric points in the north-blank fibre."""
    if m < 2 or n < 2:
        raise ValueError("both dimensions must be at least two")
    if m % 2 and n % 2:
        free_bits = (m - 1) * (n - 1)
        return m * n * 2 ** (free_bits + 1), m * 2 ** (free_bits // 2 + 1), 0
    if m % 2:
        return count_mixed(n, m, 1)
    if n % 2:
        return count_mixed(m, n, 0)
    volume = half_turn = reflection = 0
    for x, y, heading, rows, columns in blank_margins(m, n):
        volume += count_matrices(tuple(sorted(rows)), tuple(sorted(columns)), "sum")
        rot, ref = symmetric_counts(rows, columns, x, y, heading)
        half_turn += rot
        reflection += ref
    return volume, half_turn, reflection


def brute_force(m, n):
    """Enumerate all grids of a small board using the original labels."""
    size = m * n
    if size > 16:
        raise ValueError("the brute-force check is restricted to 16 cells")
    invariants = all_invariants(m, n)
    blank = [0] * size
    targets = defaultdict(list)
    for x, y, heading, delta in color_targets(invariants, blank, m, n):
        targets[tuple(delta[name] for name, _, _ in invariants)].append((x, y, heading))

    # Compute grid contributions directly from the affine invariants.
    weights = []
    changes = []
    for cell in range(size):
        grid = blank.copy()
        grid[cell] = 1
        changes.append([(evaluate(grid, 0, 0, 0) - evaluate(blank, 0, 0, 0)) % modulus
                        for _, modulus, evaluate in invariants])
    for index, (_, modulus, _) in enumerate(invariants):
        masks = [0] * modulus
        for cell, values in enumerate(changes):
            masks[values[index]] |= 1 << cell
        weights.append((modulus, masks))

    volume = symmetric = 0
    for bits in range(1 << size):
        delta = tuple(sum(k * (bits & mask).bit_count() for k, mask in enumerate(masks) if k) % modulus
                      for modulus, masks in weights)
        positions = targets.get(delta)
        if positions:
            grid = tuple((bits >> cell) & 1 for cell in range(size))
            volume += len(positions)
            symmetric += sum(is_symmetric(m, n, (grid, x, y, h)) for x, y, h in positions)
    return volume, symmetric


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("rows", type=int, nargs="?", default=10)
    parser.add_argument("columns", type=int, nargs="?", default=10)
    parser.add_argument("--check", action="store_true", help="compare with full enumeration on tiny tori")
    args = parser.parse_args()
    if args.check:
        boards = ((2, 2), (2, 4), (4, 4), (3, 3), (3, 5), (5, 3),
                  (2, 5), (5, 2), (3, 4), (4, 3))
        for m, n in boards:
            volume, half, reflection = count_blank_class(m, n)
            assert (volume, half + reflection) == brute_force(m, n), (m, n, volume, half, reflection)
        print("Brute force agreed on " + ", ".join(f"{m}x{n}" for m, n in boards) + ".")
    try:
        volume, half_turn, reflection = count_blank_class(args.rows, args.columns)
    except ValueError as error:
        parser.error(str(error))
    scale = Fraction(2 * volume, half_turn + reflection)
    print(f"{args.rows}x{args.columns}: volume={volume}, half_turn={half_turn}, "
          f"complement_reflection={reflection}")
    try:
        approximation = f" ({float(scale):.12g})"
    except OverflowError:
        approximation = ""
    print(f"z={scale}{approximation}")


if __name__ == "__main__":
    main()
