// Langton's ant on an m-by-n torus: recurrence of the blank state.
//
// Conventions
//   Cells (x, y) with 0 <= x < n (columns) and 0 <= y < m (rows); index y*n + x.
//   Headings 0 = N (y+1), 1 = E (x+1), 2 = S (y-1), 3 = W (x-1).
//   One step: read colour c of the ant's cell; turn right if c = 0 (white), left if c = 1;
//   flip the cell; move one cell forward. The ant starts at (0, 0) facing N.
//
// Usage
//   node ant.js blank M N [maxSteps]        orbit of the blank torus
'use strict';

function makeTorus(m, n) {
  const size = m * n;
  const nbr = new Int32Array(size * 4);
  for (let y = 0; y < m; y++) for (let x = 0; x < n; x++) {
    const i = y * n + x;
    nbr[4 * i + 0] = ((y + 1) % m) * n + x;
    nbr[4 * i + 1] = y * n + ((x + 1) % n);
    nbr[4 * i + 2] = ((y + m - 1) % m) * n + x;
    nbr[4 * i + 3] = y * n + ((x + n - 1) % n);
  }
  return { m, n, size, nbr };
}

// Linear parts of the symmetry group acting on headings and cell offsets.
// A symmetry is (L, b, flip): position p -> L p + b, heading h -> L h, colours complemented if flip.
// Rotations keep colours; reflections swap left and right turns, so they are paired with
// complementing every colour. On a non-square torus only the four maps preserving the axes exist.
const LINEAR = [
  // [a, b, c, d] for (x, y) -> (a x + b y, c x + d y); flip = det < 0
  { name: 'id', M: [1, 0, 0, 1], flip: 0 },
  { name: 'rot90', M: [0, 1, -1, 0], flip: 0 },   // N -> E (clockwise)
  { name: 'rot180', M: [-1, 0, 0, -1], flip: 0 },
  { name: 'rot270', M: [0, -1, 1, 0], flip: 0 },
  { name: 'refx', M: [-1, 0, 0, 1], flip: 1 },     // x -> -x, keeps N
  { name: 'refy', M: [1, 0, 0, -1], flip: 1 },     // y -> -y, N -> S
  { name: 'diag', M: [0, 1, 1, 0], flip: 1 },      // swap axes
  { name: 'anti', M: [0, -1, -1, 0], flip: 1 },
];
const HEADVEC = [[0, 1], [1, 0], [0, -1], [-1, 0]];
function headOf(v) { return HEADVEC.findIndex(h => h[0] === v[0] && h[1] === v[1]); }
function applyM(M, v) { return [M[0] * v[0] + M[1] * v[1], M[2] * v[0] + M[3] * v[1]]; }

// Group element mapping the blank start state (white, (0,0), N) to (colour u, position p, heading d).
function elementFor(T, u, p, d) {
  for (const L of LINEAR) {
    if (L.flip !== u) continue;
    if (T.m !== T.n && (L.M[1] !== 0)) continue; // axis swaps need a square torus
    if (headOf(applyM(L.M, [0, 1])) !== d) continue;
    return { L, b: p };
  }
  return null;
}

// Order of the element g on the blank state: least k with g^k(B) = B.
function orderOnBlank(T, g) {
  let p = [0, 0], h = [0, 1], u = 0;
  for (let k = 1; k <= 64 * T.size; k++) {
    p = applyM(g.L.M, p); p = [p[0] + g.b[0], p[1] + g.b[1]];
    p = [((p[0] % T.n) + T.n) % T.n, ((p[1] % T.m) + T.m) % T.m];
    h = applyM(g.L.M, h); u ^= g.L.flip;
    if (u === 0 && p[0] === 0 && p[1] === 0 && h[0] === 0 && h[1] === 1) return k;
  }
  return -1;
}

function blankOrbit(m, n, maxSteps) {
  const T = makeTorus(m, n);
  const { size, nbr } = T;
  const g = new Uint8Array(size);
  let i = 0, d = 0, black = 0, maxBlack = 0, t = 0;
  const square = m === n;
  const t0 = Date.now();
  const quasi = [];
  while (t < maxSteps) {
    const c = g[i];
    d = (d + 1 + 2 * c) & 3;
    g[i] = c ^ 1;
    black += 1 - 2 * c;
    i = nbr[4 * i + d];
    t++;
    if (black > maxBlack) maxBlack = black;
    if (black === 0 || black === size) {
      if (square || d === 0 || d === 2) {
        const u = black === 0 ? 0 : 1;
        const p = [i % n, Math.floor(i / n)];
        const el = elementFor(T, u, p, d);
        quasi.push({ t, colour: u, x: p[0], y: p[1], heading: 'NESW'[d], element: el.L.name });
        const k = orderOnBlank(T, el);
        const secs = (Date.now() - t0) / 1000;
        return { m, n, reducedPeriod: t, fullPeriod: t * k, multiplier: k, maxBlack, element: { L: el.L.name, b: el.b }, state: quasi[0], seconds: secs, stepsPerSecond: Math.round(t / Math.max(secs, 1e-3)) };
      }
    }
  }
  return { m, n, reducedPeriod: null, stepsSimulated: t, maxBlack };
}

if (require.main === module) {
  const [mode, ...a] = process.argv.slice(2);
  if (mode === 'blank') {
    const m = +a[0], n = +a[1], maxSteps = a[2] ? Number(a[2]) : 1e12;
    console.log(JSON.stringify(blankOrbit(m, n, maxSteps)));
  } else {
    console.error('usage: node ant.js blank M N [maxSteps]');
    process.exit(2);
  }
}

module.exports = { makeTorus, blankOrbit, LINEAR, elementFor, orderOnBlank };
