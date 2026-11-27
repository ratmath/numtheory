// Exact linear algebra on small BigInt matrices for the number-field code:
// Hermite normal forms of lattices, row spaces and kernels modulo a prime,
// determinants and rational inverses. Rows are vectors throughout.

import { gcd, inverse } from "./int.js";

export const modp = (a, p) => {
  a %= p;
  return a < 0n ? a + p : a;
};
export const unitVec = (i, n) => Array.from({ length: n }, (_, j) => (i === j ? 1n : 0n));
export const identity = (n) => Array.from({ length: n }, (_, i) => unitVec(i, n));
export const isZeroVec = (v) => v.every((x) => x === 0n);
const fdivr = (a, b) => {
  // floor(a / b) for b > 0
  const q = a / b;
  return a % b !== 0n && a < 0n ? q - 1n : q;
};

/**
 * The lower-triangular Hermite normal form of the lattice spanned by integer
 * rows of length n (the rows must span a lattice of full rank n): n rows, row
 * i zero beyond column i, a positive diagonal, and each entry below the
 * diagonal reduced into [0, the diagonal entry of its column).
 */
export function hnfLower(rows, n) {
  const A = rows.filter((r) => !isZeroVec(r)).map((r) => r.slice());
  const out = new Array(n);
  let free = A;
  for (let col = n - 1; col >= 0; col--) {
    // gcd the rows that are zero beyond `col` … every free row is, by induction
    let piv = null;
    const rest = [];
    for (const r of free) {
      if (r[col] === 0n) {
        rest.push(r);
        continue;
      }
      if (!piv) {
        piv = r;
        continue;
      }
      // Euclid on (piv, r) in column col
      let a = piv;
      let b = r;
      while (b[col] !== 0n) {
        const q = fdivr(a[col], b[col] < 0n ? -b[col] : b[col]) * (b[col] < 0n ? -1n : 1n);
        const t = a.map((x, j) => x - q * b[j]);
        a = b;
        b = t;
      }
      piv = a;
      if (!isZeroVec(b)) rest.push(b);
    }
    if (!piv) throw new RangeError("the rows do not span a lattice of full rank");
    if (piv[col] < 0n) piv = piv.map((x) => -x);
    out[col] = piv;
    // keep the free rows small: reduce them by the pivot is not needed (they are zero in col)
    free = rest;
  }
  for (let j = 1; j < n; j++) {
    for (let i = j - 1; i >= 0; i--) {
      const q = fdivr(out[j][i], out[i][i]);
      if (q !== 0n) for (let k = 0; k <= i; k++) out[j][k] -= q * out[i][k];
    }
  }
  return out;
}

/** The reduced row echelon form mod p of the span of the rows: { rows, pivots }. A canonical form of the subspace. */
export function spanMod(vectors, p) {
  const A = vectors.map((r) => r.map((x) => modp(x, p)));
  const c = A.length ? A[0].length : 0;
  const pivots = [];
  let rank = 0;
  for (let col = 0; col < c && rank < A.length; col++) {
    let piv = -1;
    for (let i = rank; i < A.length; i++) {
      if (A[i][col] !== 0n) {
        piv = i;
        break;
      }
    }
    if (piv < 0) continue;
    [A[rank], A[piv]] = [A[piv], A[rank]];
    const inv = inverse(A[rank][col], p);
    A[rank] = A[rank].map((x) => (x * inv) % p);
    for (let i = 0; i < A.length; i++) {
      if (i === rank || A[i][col] === 0n) continue;
      const f = A[i][col];
      A[i] = A[i].map((x, j) => modp(x - f * A[rank][j], p));
    }
    pivots.push(col);
    rank++;
  }
  return { rows: A.slice(0, rank), pivots, p };
}
/** v reduced modulo a subspace from spanMod: zero exactly when v lies in it. */
export function reduceMod(v, S) {
  const p = S.p;
  let w = v.map((x) => modp(x, p));
  S.rows.forEach((row, k) => {
    const f = w[S.pivots[k]];
    if (f !== 0n) w = w.map((x, j) => modp(x - f * row[j], p));
  });
  return w;
}
export const sameSpace = (S, T) => S.rows.length === T.rows.length && S.rows.every((r, i) => r.every((x, j) => x === T.rows[i][j]));

/** A basis of { x : x·M ≡ 0 mod p } for a matrix M given by its rows. */
export function kernelMod(M, p) {
  const r = M.length;
  if (!r) return [];
  const c = M[0].length;
  const A = M.map((row, i) => [...row.map((x) => modp(x, p)), ...unitVec(i, r)]);
  let rank = 0;
  for (let col = 0; col < c && rank < r; col++) {
    let piv = -1;
    for (let i = rank; i < r; i++) {
      if (A[i][col] !== 0n) {
        piv = i;
        break;
      }
    }
    if (piv < 0) continue;
    [A[rank], A[piv]] = [A[piv], A[rank]];
    const inv = inverse(A[rank][col], p);
    for (let i = rank + 1; i < r; i++) {
      if (A[i][col] === 0n) continue;
      const f = (A[i][col] * inv) % p;
      A[i] = A[i].map((x, j) => modp(x - f * A[rank][j], p));
    }
    rank++;
  }
  return A.slice(rank).map((row) => row.slice(c));
}

/** The determinant of a square integer matrix (Bareiss, fraction-free). */
export function det(M) {
  const n = M.length;
  if (n === 0) return 1n;
  const A = M.map((r) => r.slice());
  let sign = 1n;
  let prev = 1n;
  for (let k = 0; k < n - 1; k++) {
    if (A[k][k] === 0n) {
      const s = A.findIndex((r, i) => i > k && r[k] !== 0n);
      if (s < 0) return 0n;
      [A[k], A[s]] = [A[s], A[k]];
      sign = -sign;
    }
    for (let i = k + 1; i < n; i++) {
      for (let j = k + 1; j < n; j++) A[i][j] = (A[i][j] * A[k][k] - A[i][k] * A[k][j]) / prev;
      A[i][k] = 0n;
    }
    prev = A[k][k];
  }
  return sign * A[n - 1][n - 1];
}

/** The inverse of a square integer matrix as { inv, den }: inv·M = den·I, den > 0, entries and den without a common factor. */
export function inverseInt(M) {
  const n = M.length;
  // Gauss–Jordan on fractions kept as [numerator, denominator]
  const A = M.map((r, i) => [...r.map((x) => [x, 1n]), ...unitVec(i, n).map((x) => [x, 1n])]);
  const red = ([a, b]) => {
    if (a === 0n) return [0n, 1n];
    const g = gcd(a, b);
    return b < 0n ? [-a / g, -b / g] : [a / g, b / g];
  };
  const sub = ([a, b], [c, d]) => red([a * d - c * b, b * d]);
  const mul = ([a, b], [c, d]) => red([a * c, b * d]);
  const div = ([a, b], [c, d]) => red([a * d, b * c]);
  for (let c = 0; c < n; c++) {
    const s = A.findIndex((r, i) => i >= c && r[c][0] !== 0n);
    if (s < 0) throw new RangeError("the matrix is singular");
    [A[c], A[s]] = [A[s], A[c]];
    const pv = A[c][c];
    A[c] = A[c].map((x) => div(x, pv));
    for (let i = 0; i < n; i++) {
      if (i === c || A[i][c][0] === 0n) continue;
      const f = A[i][c];
      A[i] = A[i].map((x, j) => sub(x, mul(f, A[c][j])));
    }
  }
  let den = 1n;
  for (const row of A) for (let j = n; j < 2 * n; j++) den = (den * row[j][1]) / gcd(den, row[j][1]);
  return { inv: A.map((row) => row.slice(n).map(([a, b]) => a * (den / b))), den };
}

/** row · matrix for integer vectors. */
export function rowTimes(v, M) {
  const c = M[0].length;
  const out = new Array(c).fill(0n);
  for (let i = 0; i < v.length; i++) {
    if (v[i] === 0n) continue;
    for (let j = 0; j < c; j++) out[j] += v[i] * M[i][j];
  }
  return out;
}
export const content = (v) => v.reduce((g, x) => gcd(g, x), 0n);
