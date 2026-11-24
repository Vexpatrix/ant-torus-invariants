"""Independent verifier for the ant-recurrence results (Python standard library only).

Checks for "Affine invariants and exact recurrence scales for Langton's ant on a torus":
  A. Blank reduced periods on small tori, including the manuscript table.
  B. Mod-2 invariant dimensions from the affine classification (Theorem 1).
  C. Z/4 module types from the same theorem, by reduction of the conservation matrix.
  D. Conservation of the explicit invariants, by checking every local transition.
  E. The mean period identity on 3x3, where every reduced orbit is symmetric.
  F. The exact blank-class scales on 3x3, 3x4, and 4x3, using the explicit invariants as labels.
  G. The blank label, half-turn symmetry, and period of the 24-step construction on 6x6 and 8x6.

    python verify_ant.py [--quick]
"""
import argparse
import itertools
import sys

DIRS = [(0, 1), (1, 0), (0, -1), (-1, 0)]  # N E S W; a heading is the direction of the last move


def step(m, n, g, x, y, h):
    """One step in place on the grid list g (index y*n + x). Returns new (x, y, h)."""
    c = g[y * n + x]
    h = (h + 1) % 4 if c == 0 else (h + 3) % 4
    g[y * n + x] = 1 - c
    dx, dy = DIRS[h]
    return (x + dx) % n, (y + dy) % m, h


# ---------------------------------------------------------------- A. blank periods
def blank_reduced_period(m, n, limit):
    g = [0] * (m * n)
    x = y = h = 0
    black = 0
    for t in range(1, limit + 1):
        c = g[y * n + x]
        black += 1 - 2 * c
        x, y, h = step(m, n, g, x, y, h)
        if black in (0, m * n) and (m == n or h % 2 == 0):
            return t
    return None


# ---------------------------------------------------------------- linear-invariant equations
def equations(m, n):
    """Rows (as dicts var -> coefficient) of the conservation system for Q = sum w g + f.
    Variables: w(c) for cells c, then f(p, h) at index size + 4p + h."""
    size = m * n
    rows = []
    for y in range(m):
        for x in range(n):
            p = y * n + x
            for h in range(4):
                for turn, colour in ((1, 0), (3, 1)):  # right on white, left on black
                    h2 = (h + turn) % 4
                    dx, dy = DIRS[h2]
                    q = ((y + dy) % m) * n + (x + dx) % n
                    row = {}
                    row[size + 4 * q + h2] = row.get(size + 4 * q + h2, 0) + 1
                    row[size + 4 * p + h] = row.get(size + 4 * p + h, 0) - 1
                    row[p] = row.get(p, 0) + (1 if colour == 0 else -1)  # cell rises / falls
                    rows.append(row)
    return rows, 5 * size


def gf2_nullity(m, n):
    rows, nv = equations(m, n)
    basis = {}  # pivot bit -> row
    for r in rows:
        v = 0
        for k, c in r.items():
            if c % 2:
                v ^= 1 << k
        while v:
            top = v.bit_length() - 1
            if top in basis:
                v ^= basis[top]
            else:
                basis[top] = v
                break
    return nv - len(basis)


def z4_module(m, n):
    """Counts (#Z/2, #Z/4) of the solution module of the system over Z/4 (Smith normal form)."""
    rows, nv = equations(m, n)
    A = [[0] * nv for _ in rows]
    for i, r in enumerate(rows):
        for k, c in r.items():
            A[i][k] = c % 4
    nr = len(A)
    used_rows, pivots = set(), {}
    for target in (1, 2):  # unit pivots, then pivots equal to 2
        while True:
            found = None
            for i in range(nr):
                if i in used_rows:
                    continue
                Ai = A[i]
                for j in range(nv):
                    if j in pivots:
                        continue
                    a = Ai[j]
                    if a and ((a % 2 == 1) if target == 1 else a == 2):
                        found = (i, j)
                        break
                if found:
                    break
            if not found:
                break
            i, j = found
            if A[i][j] == 3:
                A[i] = [(-a) % 4 for a in A[i]]
            for k in range(nr):
                if k != i and A[k][j]:
                    f = A[k][j] if target == 1 else A[k][j] // 2
                    Ak, Ai = A[k], A[i]
                    A[k] = [(Ak[t] - f * Ai[t]) % 4 for t in range(nv)]
            # clear the rest of row i by column operations (does not change the module type)
            used_rows.add(i)
            pivots[j] = target
    z2 = sum(1 for v in pivots.values() if v == 2)
    free = nv - len(pivots)
    return z2, free  # Z/2 summands (pivot 2), Z/4 summands (free columns)


def expected_module(m, n):
    if m % 2 and n % 2:
        return m + n - 2, 2
    if m % 2 == 0 and n % 2 == 0:
        return 2, m + n - 1
    if m % 2 == 0:
        return m - 1, n + 1
    return n - 1, m + 1


# ---------------------------------------------------------------- D. explicit invariants
def column_invariants(m, n):
    """List of (name, modulus, function(g, x, y, h)) on the m-by-n torus (m rows, n columns)."""
    inv = []

    def sgn(v):
        return 1 if v % 2 == 0 else -1

    inv.append(('J', 4, lambda g, x, y, h: (sum(g) - h) % 4))
    for x0 in range(n):
        inv.append((f'Q{x0}', 2, lambda g, x, y, h, x0=x0:
                    (sum(g[yy * n + x0] for yy in range(m)) + (x == x0 and h % 2 == 0)) % 2))
    if m % 2 == 0 and n % 2 == 0:
        inv.append(('align', 2, lambda g, x, y, h: (x + y + (h % 2 == 0)) % 2))
    if m % 2 == 0:
        for x0 in range(n):
            x1 = (x0 + 1) % n

            def pair(g, x, y, h, x0=x0, x1=x1):
                s = sum(sgn(yy) * (g[yy * n + x0] + g[yy * n + x1]) for yy in range(m))
                a, b, sy = int(x == x0), int(x == x1), sgn(y)
                f = {1: 2 * b, 3: 2 * a, 0: sy * (b - a), 2: sy * (a - b)}[h]
                return (s + f) % 4
            inv.append((f'P{x0}', 4, pair))
            if n % 2 == 1:
                def single(g, x, y, h, x0=x0):
                    s = sum(sgn(yy) * g[yy * n + x0] for yy in range(m))
                    c = (x - x0 - 1) % n
                    a, sy = int(x == x0), sgn(y)
                    loc = {0: -sy * a, 2: sy * a, 3: 2 * a, 1: 0}[h]
                    return (s + 2 * (c + y + (h % 2 == 0)) + loc) % 4
                inv.append((f'U{x0}', 4, single))
    return inv


REFLECT = [1, 0, 3, 2]


def all_invariants(m, n):
    inv = column_invariants(m, n)
    for name, k, fn in column_invariants(n, m):
        if name in ('J', 'align'):
            continue

        def via_transpose(g, x, y, h, fn=fn):
            g2 = [0] * (m * n)
            for yy in range(m):
                for xx in range(n):
                    g2[xx * m + yy] = 1 - g[yy * n + xx]
            return fn(g2, y, x, REFLECT[h])
        inv.append(('row:' + name, k, via_transpose))
    return inv


def exhaustive_transitions(m, n, rng):
    inv = all_invariants(m, n)
    bad = set()
    for y in range(m):
        for x in range(n):
            for h in range(4):
                for c in (0, 1):
                    g = [rng.getrandbits(1) for _ in range(m * n)]
                    g[y * n + x] = c
                    before = [fn(g, x, y, h) for _, _, fn in inv]
                    x2, y2, h2 = step(m, n, g, x, y, h)
                    after = [fn(g, x2, y2, h2) for _, _, fn in inv]
                    bad.update(inv[i][0] for i in range(len(inv)) if before[i] != after[i])
    return len(inv), bad


# ---------------------------------------------------------------- E, F. enumeration on tiny tori
def symmetries(m, n):
    """Torus symmetries as functions on states (g tuple, x, y, h)."""
    lin = [((1, 0, 0, 1), 0), ((-1, 0, 0, -1), 0), ((-1, 0, 0, 1), 1), ((1, 0, 0, -1), 1)]
    if m == n:
        lin += [((0, 1, -1, 0), 0), ((0, -1, 1, 0), 0), ((0, 1, 1, 0), 1), ((0, -1, -1, 0), 1)]
    out = []
    for (a, b, c, d), flip in lin:
        for tx in range(n):
            for ty in range(m):
                def apply(state, a=a, b=b, c=c, d=d, flip=flip, tx=tx, ty=ty):
                    g, x, y, h = state
                    g2 = [0] * (m * n)
                    for yy in range(m):
                        for xx in range(n):
                            X, Y = (a * xx + b * yy + tx) % n, (c * xx + d * yy + ty) % m
                            g2[Y * n + X] = g[yy * n + xx] ^ flip
                    hx, hy = DIRS[h]
                    h2 = DIRS.index((a * hx + b * hy, c * hx + d * hy))
                    return (tuple(g2), (a * x + b * y + tx) % n, (c * x + d * y + ty) % m, h2)
                out.append(apply)
    return out


def is_symmetric(m, n, state):
    """Fixed by gamma o M for the half-turn or the complementing reflection allowed at the ant."""
    g, x, y, h = state
    dx, dy = DIRS[h]
    vx, vy = 2 * x - dx, 2 * y - dy
    if all(g[yy * n + xx] == g[((vy - yy) % m) * n + (vx - xx) % n] for yy in range(m) for xx in range(n)):
        return True
    cells = []
    for yy in range(m):
        for xx in range(n):
            e = ((vy - yy) % m) * n + xx if h % 2 == 0 else yy * n + (vx - xx) % n
            cells.append((yy * n + xx, e))
    return all(i != e and g[i] != g[e] for i, e in cells)


def enumerate_states(m, n):
    for bits in range(2 ** (m * n)):
        g = tuple((bits >> k) & 1 for k in range(m * n))
        for y in range(m):
            for x in range(n):
                for h in range(4):
                    yield (g, x, y, h)


def canonical(syms, state):
    return min(s(state) for s in syms)


def check_identity_3x3():
    m = n = 3
    syms = symmetries(m, n)
    reduced = set()
    sym_reps = set()
    for st in enumerate_states(m, n):
        c = canonical(syms, st)
        reduced.add(c)
        if is_symmetric(m, n, st):
            sym_reps.add(c)
    total = 0
    for c in sym_reps:
        g, x, y, h = c
        gl = list(g)
        t = 0
        while True:
            x, y, h = step(m, n, gl, x, y, h)
            t += 1
            if canonical(syms, (tuple(gl), x, y, h)) == c:
                break
        total += t
    return len(reduced), len(sym_reps), total / len(sym_reps), 2 * len(reduced) / len(sym_reps)


def exact_class_scale(m, n):
    inv = all_invariants(m, n)
    label = lambda st: tuple(fn(list(st[0]), st[1], st[2], st[3]) for _, _, fn in inv)
    blank = (tuple([0] * (m * n)), 0, 0, 0)
    target = {label(s(blank)) for s in symmetries(m, n)}
    vol = sym = 0
    for st in enumerate_states(m, n):
        if label(st) in target:
            vol += 1
            sym += is_symmetric(m, n, st)
    return vol, sym, 2 * vol / sym


def check_short_orbit(m):
    """Check the 24-step construction of Proposition 2 on an m-by-6 torus."""
    n = 6
    white = (0,) * n
    rows = {
        0: (0, 0, 0, 0, 1, 1),
        1: (0, 0, 0, 0, 1, 1),
        m - 2: (0, 1, 1, 0, 0, 0),
        m - 1: (0, 1, 1, 0, 0, 0),
    }
    grid = [colour for y in range(m) for colour in rows.get(y, white)]
    start = (tuple(grid), 0, 0, 0)
    blank = ((0,) * (m * n), 0, 0, 0)
    same_label = all(fn(*start) == fn(*blank) for _, _, fn in all_invariants(m, n))
    half_turn = all(
        grid[y * n + x] == grid[((-1 - y) % m) * n + (-x) % n]
        for y in range(m) for x in range(n)
    )
    images = {symmetry(start) for symmetry in symmetries(m, n)}
    x = y = h = 0
    reduced_period = None
    for elapsed in range(1, 25):
        x, y, h = step(m, n, grid, x, y, h)
        state = (tuple(grid), x, y, h)
        if reduced_period is None and state in images:
            reduced_period = elapsed
    return same_label, half_turn, reduced_period, state == start


def main():
    import random
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('--quick', action='store_true', help='skip the slowest checks')
    args = parser.parse_args()
    ok = True

    def report(label, passed, detail=''):
        nonlocal ok
        ok &= passed
        print(f"{'PASS' if passed else 'FAIL'}  {label}  {detail}")

    periods = {
        (3, 3): 22, (4, 4): 48, (4, 6): 696, (6, 4): 96,
        (5, 5): 2342, (6, 6): 4592, (5, 6): 16804, (6, 5): 4184,
        (8, 5): 15800,
    }
    if not args.quick:
        periods.update({(7, 7): 9166514, (8, 8): 5751232})
    for (m, n), expected in periods.items():
        t = blank_reduced_period(m, n, 10 ** 8)
        report(f'A blank reduced period {m}x{n}', t == expected, f'{t}')

    for m in range(2, 9):
        for n in range(2, 9):
            d = gf2_nullity(m, n)
            e = m + n + (1 if m % 2 == 0 and n % 2 == 0 else 0)
            report(f'B mod-2 dimension {m}x{n}', d == e, f'{d}')

    top = 6 if args.quick else 7
    for m in range(2, top + 1):
        for n in range(2, top + 1):
            got = z4_module(m, n)
            report(f'C Z/4 module {m}x{n}', got == expected_module(m, n), f'(Z/2, Z/4) = {got}')

    rng = random.Random(2026)
    for m in range(2, 9):
        for n in range(2, 9):
            count, bad = exhaustive_transitions(m, n, rng)
            report(f'D explicit invariants {m}x{n}', not bad, f'{count} invariants' + (f' violated {sorted(bad)}' if bad else ''))

    red, sym, mean, z = check_identity_3x3()
    report('E exact identity 3x3', abs(mean - z) < 1e-12, f'reduced {red}, symmetric {sym}, mean {mean}, 2N/g {z}')

    for (m, n), expected in (((3, 3), 24), ((3, 4), 104), ((4, 3), 52)):
        vol, sym, z = exact_class_scale(m, n)
        report(f'F exact blank class scale {m}x{n}', z == expected, f'class {vol}, symmetric {sym}, z {z}')

    for m in (6, 8):
        label, half_turn, period, exact_return = check_short_orbit(m)
        report(f'G 24-step construction {m}x6', label and half_turn and period == 24 and exact_return,
               f'blank label {label}, half-turn {half_turn}, reduced period {period}, exact return {exact_return}')

    print('ALL PASSED' if ok else 'FAILURES')
    sys.exit(0 if ok else 1)


if __name__ == '__main__':
    main()
