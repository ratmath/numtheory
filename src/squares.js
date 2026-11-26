// Sums of squares: two squares (Fermat, with every representation through
// Gaussian integers), three squares (Legendre's condition), four squares
// (Lagrange), the counting functions r₂ and r₄ by their formulas, rₖ by direct
// counting, and Pythagorean triples.

import { sqrtMod } from "@ratmath/algebraic";
import { big, isqrt, refuse } from "./int.js";
import { factor } from "./primes.js";
import { divisors } from "./arith.js";

const nat = (n) => {
  const N = big(n, "n");
  if (N < 0n) throw new RangeError("n must not be negative");
  return N;
};
const sortPair = ([a, b]) => {
  a = a < 0n ? -a : a;
  b = b < 0n ? -b : b;
  return a >= b ? [a, b] : [b, a];
};

/** A prime p ≡ 1 (mod 4) as a² + b² (a > b > 0): a square root of −1 mod p, then Euclid's algorithm until the remainder drops below √p. */
export function primeAsTwoSquares(p) {
  const P = big(p, "p");
  if (P === 2n) return [1n, 1n];
  if (P % 4n !== 1n) refuse("bad-argument", `${P} is not 1 mod 4`);
  const x = sqrtMod(P - 1n, P);
  let [r0, r1] = [P, x];
  const s = isqrt(P);
  while (r1 > s) [r0, r1] = [r1, r0 % r1];
  const b = isqrt(P - r1 * r1);
  if (r1 * r1 + b * b !== P) throw new Error(`no representation found for ${P}: is it prime?`);
  return sortPair([r1, b]);
}

const gmul = ([a, b], [c, d]) => [a * c - b * d, a * d + b * c];
const gpow = (z, k) => {
  let r = [1n, 0n];
  for (let i = 0; i < k; i++) r = gmul(r, z);
  return r;
};

/**
 * Every way to write n as a² + b² with a ≥ b ≥ 0, sorted by a descending
 * (empty when there is none). By unique factorisation in the Gaussian
 * integers: each prime p ≡ 1 (mod 4) splits as (a + bi)(a − bi), and the
 * representations come from choosing how many factors of each kind to take.
 */
export function twoSquaresAll(n) {
  const N = nat(n);
  if (N === 0n) return [[0n, 0n]];
  let base = [[1n, 0n]];
  for (const { prime, exponent } of factor(N, { prove: false })) {
    if (prime === 2n) base = base.map((z) => gmul(z, gpow([1n, 1n], exponent)));
    else if (prime % 4n === 3n) {
      if (exponent % 2) return [];
      base = base.map((z) => gmul(z, [prime ** BigInt(exponent / 2), 0n]));
    } else {
      const [a, b] = primeAsTwoSquares(prime);
      const next = [];
      for (const z of base) for (let k = 0; k <= exponent; k++) next.push(gmul(z, gmul(gpow([a, b], k), gpow([a, -b], exponent - k))));
      base = next;
    }
  }
  const seen = new Map();
  for (const z of base) {
    const p = sortPair(z);
    seen.set(`${p[0]},${p[1]}`, p);
  }
  return [...seen.values()].sort((x, y) => (x[0] > y[0] ? -1 : 1));
}
/** One way to write n as a² + b² (a ≥ b ≥ 0), or null when there is none (twoSquaresReason says why). */
export function twoSquares(n) {
  const all = twoSquaresAll(n);
  return all.length ? all[0] : null;
}
/** Why n is or is not a sum of two squares: { possible, reason }. */
export function twoSquaresReason(n) {
  const N = nat(n);
  if (N === 0n) return { possible: true, reason: "0 = 0² + 0²" };
  const bad = factor(N, { prove: false }).find((f) => f.prime % 4n === 3n && f.exponent % 2);
  return bad
    ? { possible: false, reason: `${bad.prime} ≡ 3 (mod 4) divides ${N} to the odd power ${bad.exponent}` }
    : { possible: true, reason: "every prime ≡ 3 (mod 4) appears to an even power" };
}
/** r₂(n): the number of (x, y) in Z² with x² + y² = n, signs and order counted: 4·(d₁(n) − d₃(n)). */
export function r2(n) {
  const N = nat(n);
  if (N === 0n) return 1n;
  let r = 4n;
  for (const { prime, exponent } of factor(N, { prove: false })) {
    if (prime % 4n === 3n) {
      if (exponent % 2) return 0n;
    } else if (prime % 4n === 1n) r *= BigInt(exponent + 1);
  }
  return r;
}

const excluded = (N) => {
  let m = N;
  while (m > 0n && m % 4n === 0n) m /= 4n;
  return m % 8n === 7n;
};
/** n as a² + b² + c² (a ≥ b ≥ c ≥ 0), or null exactly when n = 4^a(8b + 7) (Legendre and Gauss). */
export function threeSquares(n) {
  const N = nat(n);
  if (excluded(N)) return null;
  for (let x = isqrt(N); x >= 0n; x--) {
    const t = twoSquares(N - x * x);
    if (t) return [x, ...t].sort((a, b) => (a > b ? -1 : a < b ? 1 : 0));
  }
  throw new Error("unreachable: Legendre's three-square theorem");
}
/** n as a sum of four squares (Lagrange: always possible), a ≥ b ≥ c ≥ d ≥ 0. */
export function fourSquares(n) {
  const N = nat(n);
  if (!excluded(N)) return [...threeSquares(N), 0n];
  // n = 4^a·m with m ≡ 7 (mod 8): m − 1 ≡ 6 (mod 8) is a sum of three squares, then scale by 2^a
  let m = N, s = 1n;
  while (m % 4n === 0n) [m, s] = [m / 4n, s * 2n];
  return [1n, ...threeSquares(m - 1n)].map((v) => v * s).sort((a, b) => (a > b ? -1 : a < b ? 1 : 0));
}
/** r₄(n) = 8·Σ d over the divisors d of n not divisible by 4 (Jacobi). */
export function r4(n) {
  const N = nat(n);
  if (N === 0n) return 1n;
  return 8n * divisors(N).filter((d) => d % 4n !== 0n).reduce((s, d) => s + d, 0n);
}
/** r_k(n): the number of (x₁, …, x_k) in Z^k with Σ xᵢ² = n, by counting (n up to 10⁵, k up to 12). */
export function rk(n, k) {
  const N = Number(nat(n));
  const K = Number(big(k, "k"));
  if (N > 100000 || K > 12 || K < 1) refuse("too-large", "rk counts directly: n ≤ 10⁵ and 1 ≤ k ≤ 12");
  // coefficients of θ(q)^k up to q^n, θ = 1 + 2Σ q^(m²)
  const theta = new Array(N + 1).fill(0n);
  theta[0] = 1n;
  for (let m = 1; m * m <= N; m++) theta[m * m] += 2n;
  let acc = theta.slice();
  for (let j = 1; j < K; j++) {
    const next = new Array(N + 1).fill(0n);
    for (let i = 0; i <= N; i++) {
      if (!acc[i]) continue;
      for (let m = 0; i + m * m <= N; m++) next[i + m * m] += acc[i] * (m === 0 ? 1n : 2n);
    }
    acc = next;
  }
  return acc[N];
}

// ---- Pythagorean triples ----

/**
 * The primitive Pythagorean triples [a, b, c] (a² + b² = c², gcd 1, a odd)
 * with c ≤ limit, by Euclid's formula a = m² − n², b = 2mn, c = m² + n² for
 * coprime m > n > 0 of opposite parity. Sorted by c, then a.
 */
export function pythagoreanTriples(limit) {
  const L = nat(limit);
  if (L > 100000000n) refuse("too-large", "the limit is 10⁸");
  const out = [];
  for (let m = 2n; m * m + 1n <= L; m++) {
    for (let k = 1n + (m % 2n); k < m; k += 2n) {
      const c = m * m + k * k;
      if (c > L) break;
      let [x, y] = [m, k];
      while (y) [x, y] = [y, x % y];
      if (x === 1n) out.push([m * m - k * k, 2n * m * k, c]);
    }
  }
  return out.sort((u, v) => (u[2] < v[2] ? -1 : u[2] > v[2] ? 1 : u[0] < v[0] ? -1 : 1));
}
/** Every triple [a, b, c] with hypotenuse c (a < b), from the representations of c² as a sum of two squares. */
export function triplesWithHypotenuse(c) {
  const C = nat(c);
  return twoSquaresAll(C * C)
    .filter(([, b]) => b > 0n)
    .map(([a, b]) => [b, a, C])
    .sort((u, v) => (u[0] < v[0] ? -1 : 1));
}
/** Every triple [a, b, c] with a leg equal to a, from the factor pairs of a² of equal parity. */
export function triplesWithLeg(a) {
  const A = nat(a);
  if (A < 3n) return [];
  const sq = A * A;
  const out = [];
  for (const d of divisors(sq)) {
    const e = sq / d;
    if (d >= e || (e - d) % 2n) continue;
    out.push([A, (e - d) / 2n, (e + d) / 2n]);
  }
  return out.sort((u, v) => (u[1] < v[1] ? -1 : 1));
}
