// Diophantine equations: linear equations in any number of unknowns, the
// general Pell equation x² − D·y² = N, the Frobenius coin problem, and
// Egyptian fractions.

import { Rational } from "@ratmath/core";
import { big, egcd, fdiv, gcd, isqrt, isSquare, rat, refuse } from "./int.js";
import { pell } from "./cf.js";
import { divisors } from "./arith.js";

// ---- linear ----

/**
 * All integer solutions of a₁x₁ + … + aₙxₙ = c. Returns
 *   { particular, basis, g }   every solution is particular + Σ tᵢ·basisᵢ (tᵢ any integers)
 * or { solvable: false, g, reason } when g = gcd(a₁, …, aₙ) does not divide c.
 * The basis spans the integer solutions of the homogeneous equation; it is
 * found by unimodular column operations (the extended Euclidean algorithm on
 * many numbers), so it is a lattice basis, not just a set of solutions.
 */
export function linearDiophantine(coefficients, c) {
  const a = coefficients.map((x) => big(x, "a coefficient"));
  const C = big(c, "c");
  const n = a.length;
  if (!n) refuse("bad-argument", "give at least one coefficient");
  // U starts as the identity; row ops keep a·U's first entry the running gcd
  const U = a.map((_, i) => a.map((__, j) => (i === j ? 1n : 0n)));
  const v = [...a];
  const col = (j) => U.map((r) => r[j]);
  const setCol = (j, xs) => xs.forEach((x, i) => (U[i][j] = x));
  for (let j = 1; j < n; j++) {
    // combine columns 0 and j so that v[0] = gcd(v[0], v[j]), v[j] = 0
    const { g, x, y } = egcd(v[0], v[j]);
    if (g === 0n) continue;
    const c0 = col(0);
    const cj = col(j);
    const s = v[0] / g;
    const t = v[j] / g;
    setCol(0, c0.map((u, i) => x * u + y * cj[i]));
    setCol(j, c0.map((u, i) => -t * u + s * cj[i]));
    v[0] = g;
    v[j] = 0n;
  }
  let g = v[0];
  if (g < 0n) {
    g = -g;
    setCol(0, col(0).map((u) => -u));
  }
  if (g === 0n) {
    if (C !== 0n) return { solvable: false, g, reason: "every coefficient is 0 but c is not" };
    return { solvable: true, g, particular: a.map(() => 0n), basis: a.map((_, i) => a.map((__, j) => (i === j ? 1n : 0n))) };
  }
  if (C % g !== 0n) return { solvable: false, g, reason: `gcd = ${g} does not divide ${C}` };
  const k = C / g;
  return { solvable: true, g, particular: col(0).map((u) => u * k), basis: Array.from({ length: n - 1 }, (_, j) => col(j + 1)) };
}

/**
 * The non-negative solutions of a·x + b·y = c (a, b > 0): { solutions } in
 * increasing x. There are about c/(ab) of them.
 */
export function nonNegativeSolutions(a, b, c) {
  const A = big(a, "a");
  const B = big(b, "b");
  const C = big(c, "c");
  if (A <= 0n || B <= 0n) refuse("bad-argument", "a and b must be positive");
  const r = linearDiophantine([A, B], C);
  if (!r.solvable) return { solutions: [], reason: r.reason };
  const [x0, y0] = r.particular;
  const [dx, dy] = r.basis[0].map((t) => (r.basis[0][0] < 0n ? -t : t)); // dx > 0, dy < 0
  // x = x0 + t·dx ≥ 0 and y = y0 + t·dy ≥ 0
  const tMin = -fdiv(x0, dx);
  const tMax = fdiv(y0, -dy);
  if (tMax - tMin > 100000n) refuse("too-large", `there are ${tMax - tMin + 1n} solutions`);
  const out = [];
  for (let t = tMin; t <= tMax; t++) out.push([x0 + t * dx, y0 + t * dy]);
  return { solutions: out };
}

// ---- Frobenius ----

/**
 * The Frobenius number of coin values with gcd 1: the largest amount that
 * cannot be paid. Two coins: ab − a − b (Sylvester). More: shortest paths on
 * the residues modulo the smallest coin (each residue's cheapest payable
 * amount; the answer is the largest of those, minus that coin).
 */
export function frobenius(coins) {
  const cs = coins.map((x) => big(x, "a coin")).sort((a, b) => (a < b ? -1 : 1));
  if (cs.some((c) => c <= 0n)) refuse("bad-argument", "coins must be positive");
  if (cs.reduce(gcd, 0n) !== 1n) refuse("bad-argument", "the coins have a common factor, so infinitely many amounts are impossible");
  if (cs[0] === 1n) return -1n;
  if (cs.length === 2) return cs[0] * cs[1] - cs[0] - cs[1];
  const m = cs[0];
  if (m > 1000000n) refuse("too-large", "the smallest coin must be at most 10⁶");
  const M = Number(m);
  const dist = new Array(M).fill(null);
  dist[0] = 0n;
  // Dijkstra on residues: from r, a coin c leads to (r + c) mod m at cost c
  const done = new Uint8Array(M);
  for (;;) {
    let best = -1;
    for (let r = 0; r < M; r++) if (!done[r] && dist[r] !== null && (best < 0 || dist[r] < dist[best])) best = r;
    if (best < 0) break;
    done[best] = 1;
    for (const c of cs.slice(1)) {
      const r = (best + Number(c % m)) % M;
      const d = dist[best] + c;
      if (dist[r] === null || d < dist[r]) dist[r] = d;
    }
  }
  return dist.reduce((a, b) => (b > a ? b : a), 0n) - m;
}
/** The number of ways to pay an amount with unlimited coins of the given values (order ignored). */
export function coinWays(amount, coins) {
  const A = Number(big(amount, "amount"));
  if (A > 10000000) refuse("too-large", "the amount is limited to 10⁷");
  const ways = new Array(A + 1).fill(0n);
  ways[0] = 1n;
  for (const c of coins.map((x) => Number(big(x, "a coin")))) for (let s = c; s <= A; s++) ways[s] += ways[s - c];
  return ways[A];
}

// ---- Pell ----

/**
 * The integer solutions of x² − D·y² = N.
 *
 * D > 0 and not a square: the solutions fall into finitely many classes,
 * each the orbit of one fundamental solution under multiplication by the
 * fundamental solution (x₁, y₁) of x² − Dy² = 1. Fundamental solutions are
 * found in the classical bounded range for y (Nagell): for N > 0,
 * 0 ≤ y ≤ √(N(x₁ − 1)/(2D)); for N < 0, √(|N|/D) ≤ y ≤ √(|N|(x₁ + 1)/(2D)).
 * Returns { kind: "infinite", fundamental, unit, solutions(count) }.
 *
 * D < 0, D = 0 or D a square: finitely many, all listed:
 * { kind: "finite", solutions }.
 */
export function pellSolutions(D, N) {
  const d = big(D, "D");
  const n = big(N, "N");
  const norm = (x, y) => x * x - d * y * y;
  if (d <= 0n || isSquare(d)) {
    if (n === 0n && d >= 0n && isSquare(d) && d !== 0n) return { kind: "infinite", reason: `x = ${isqrt(d)}·y gives a solution for every y`, solutions: (k) => Array.from({ length: k }, (_, i) => [isqrt(d) * BigInt(i), BigInt(i)]) };
    if (d === 0n) {
      if (!isSquare(n)) return { kind: "finite", solutions: [] };
      return { kind: "infinite", reason: "y is free when D = 0", solutions: (k) => Array.from({ length: k }, (_, i) => [isqrt(n), BigInt(i)]) };
    }
    const sols = [];
    if (d < 0n) {
      // x² + |D|y² = N: bounded
      if (n < 0n) return { kind: "finite", solutions: [] };
      for (let y = 0n; -d * y * y <= n; y++) {
        const r = n + d * y * y;
        if (isSquare(r)) {
          const x = isqrt(r);
          for (const sx of x ? [x, -x] : [0n]) for (const sy of y ? [y, -y] : [0n]) sols.push([sx, sy]);
        }
        if (y > 10000000n) refuse("too-large", "the search for y passed 10⁷");
      }
    } else {
      // D = s²: (x − sy)(x + sy) = N
      const s = isqrt(d);
      if (n === 0n) return { kind: "infinite", reason: "x = ±s·y", solutions: (k) => Array.from({ length: k }, (_, i) => [s * BigInt(i), BigInt(i)]) };
      const an = n < 0n ? -n : n;
      for (const e of divisors(an)) for (const sg of [1n, -1n]) {
        const u = sg * e;
        const v = n / u; // x − sy = u, x + sy = v
        if ((u + v) % 2n !== 0n || (v - u) % (2n * s) !== 0n) continue;
        sols.push([(u + v) / 2n, (v - u) / (2n * s)]);
      }
    }
    const key = new Set();
    const uniq = sols.filter(([x, y]) => (key.has(`${x},${y}`) ? false : key.add(`${x},${y}`)));
    return { kind: "finite", solutions: uniq.sort((a, b) => (a[0] < b[0] ? -1 : a[0] > b[0] ? 1 : a[1] < b[1] ? -1 : 1)) };
  }
  const { x: x1, y: y1 } = pell(d);
  const fundamental = [];
  if (n === 0n) return { kind: "finite", solutions: [[0n, 0n]] };
  const an = n < 0n ? -n : n;
  let lo;
  let hi;
  if (n > 0n) {
    lo = 0n;
    hi = isqrt((n * (x1 - 1n)) / (2n * d));
  } else {
    lo = isqrt(an / d);
    hi = isqrt((an * (x1 + 1n)) / (2n * d));
  }
  if (hi - lo > 10000000n) refuse("too-large", `the search range for fundamental solutions has ${hi - lo} values; D or N is too large for this method`);
  for (let y = lo; y <= hi; y++) {
    const t = n + d * y * y;
    if (t < 0n || !isSquare(t)) continue;
    const x = isqrt(t);
    for (const sx of x ? [x, -x] : [0n]) {
      // (x, y) and (x, −y) may be in different classes; (−x, −y) is (x, y) times −1, its own class
      fundamental.push([sx, y]);
      if (y) fundamental.push([sx, -y]);
    }
  }
  // keep one representative per class: two solutions are in one class when (x x' − D y y')/N and (x y' − x' y)/N are integers
  const classes = [];
  for (const s of fundamental) {
    if (classes.some((c) => sameClass(c, s, d, n))) continue;
    classes.push(s);
  }
  for (const s of classes) if (norm(s[0], s[1]) !== n) throw new Error("internal check failed");
  const step = ([x, y]) => [x * x1 + d * y * y1, x * y1 + y * x1];
  return {
    kind: "infinite",
    D: d,
    N: n,
    unit: [x1, y1],
    fundamental: classes,
    /** The first `count` solutions with x, y ≥ 0, in increasing order of x. */
    solutions(count = 10) {
      const k = Number(count);
      if (!classes.length) return [];
      const inv = ([x, y]) => [x * x1 - d * y * y1, y * x1 - x * y1];
      const out = new Map();
      for (const c of classes) {
        let up = c;
        let down = c;
        for (let i = 0; i <= k + 1; i++) {
          for (const [x, y] of [up, down]) {
            const p = [x < 0n ? -x : x, y < 0n ? -y : y];
            out.set(`${p[0]},${p[1]}`, p);
          }
          up = step(up);
          down = inv(down);
        }
      }
      return [...out.values()].sort((a, b) => (a[0] < b[0] ? -1 : a[0] > b[0] ? 1 : a[1] < b[1] ? -1 : 1)).slice(0, k);
    },
  };
}
function sameClass([x, y], [u, v], d, n) {
  return (x * u - d * y * v) % n === 0n && (x * v - y * u) % n === 0n;
}
/** The norm equation N(x + y√D) = x² − D·y² = n in the order Z[√D]; the same as pellSolutions(D, n). */
export const normEquation = (D, n) => pellSolutions(D, n);

// ---- Egyptian fractions ----

const posRat = (q) => {
  const r = rat(q, "the fraction");
  if (r.numerator <= 0n) refuse("bad-argument", "the fraction must be positive");
  return r;
};
/**
 * The greedy (Fibonacci–Sylvester) expansion of a positive fraction as a
 * sum of distinct unit fractions: take the largest 1/n that fits, repeat.
 * Returns the denominators. A fraction ≥ 1 starts with 1/1, 1/2, … greedily.
 */
export function egyptianGreedy(q) {
  let r = posRat(q);
  const out = [];
  let next = 1n;
  while (r.numerator !== 0n) {
    const p = r.numerator;
    const d = r.denominator;
    let n = (d + p - 1n) / p; // ⌈d/p⌉
    if (n < next) n = next;
    out.push(n);
    r = r.subtract(new Rational(1n, n));
    next = n + 1n;
    if (out.length > 10000) refuse("too-large", "more than 10⁴ terms");
  }
  return out;
}
/**
 * A shortest expansion as distinct unit fractions, by search with a bound on
 * the denominators (default 10⁶); among the shortest, the one with the
 * smallest largest denominator. Returns the denominators, or null when none
 * exists within the bound and length 8.
 */
export function egyptianShortest(q, { maxDenominator = 1000000n, maxTerms = 8 } = {}) {
  const r = posRat(q);
  const B = big(maxDenominator, "maxDenominator");
  for (let k = 1; k <= maxTerms; k++) {
    let best = null;
    const search = (rem, from, left, acc) => {
      if (rem.numerator === 0n) return left === 0;
      if (left === 0) return false;
      const p = rem.numerator;
      const d = rem.denominator;
      if (left === 1) {
        if (p === 1n && d >= from && d <= B) {
          const cand = [...acc, d];
          if (!best || d < best[best.length - 1]) best = cand;
          return true;
        }
        return false;
      }
      // 1/n ≤ rem and left/n ≥ rem: ⌈d/p⌉ ≤ n ≤ ⌊left·d/p⌋
      let lo = (d + p - 1n) / p;
      if (lo < from) lo = from;
      const hi = (BigInt(left) * d) / p;
      let found = false;
      for (let n = lo; n <= hi && n <= B; n++) {
        if (best && n >= best[best.length - 1]) break;
        if (search(rem.subtract(new Rational(1n, n)), n + 1n, left - 1, [...acc, n])) found = true;
      }
      return found;
    };
    search(r, 1n, k, []);
    if (best) return best;
  }
  return null;
}
