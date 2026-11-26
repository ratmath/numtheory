// Congruences: linear congruences and systems of them, the Chinese remainder
// theorem for any moduli, polynomial congruences (roots mod p, lifted to p^k
// by Hensel's lemma, combined by the Chinese remainder theorem), primitive
// roots modulo any n that has one, and indices (discrete logarithms).

import { discreteLogMod, multiplicativeOrder } from "@ratmath/algebraic";
import { big, egcd, gcd, inverse, mod, powmod, refuse, val } from "./int.js";
import { factor } from "./primes.js";
import { totient } from "./arith.js";

const pos = (m, what = "the modulus") => {
  const n = big(m, what);
  if (n < 1n) throw new RangeError(`${what} must be at least 1`);
  return n;
};

/**
 * All x mod m with a·x ≡ b (mod m). Returns { solutions, modulus, g } where
 * g = gcd(a, m); there are g solutions, spaced m/g apart. When g does not
 * divide b there are none, and the result says why.
 */
export function solveLinear(a, b, m) {
  const M = pos(m);
  const A = mod(big(a, "a"), M);
  const B = mod(big(b, "b"), M);
  const { g, x } = egcd(A, M);
  if (B % g !== 0n) return { solutions: [], modulus: M, g, reason: `gcd(${A}, ${M}) = ${g} does not divide ${B}` };
  const step = M / g;
  const x0 = mod((x * (B / g)), step);
  return { solutions: Array.from({ length: Number(g) }, (_, k) => x0 + BigInt(k) * step), modulus: M, g, reduced: { x: x0, modulus: step } };
}

/**
 * The Chinese remainder theorem for any moduli: x ≡ rᵢ (mod mᵢ) for all i.
 * Returns { x, modulus } with modulus the lcm, or null when the congruences
 * contradict each other (with `conflict` naming two that do).
 */
export function crt(residues, moduli) {
  if (residues.length !== moduli.length) refuse("bad-argument", "give one modulus per residue");
  let x = 0n;
  let m = 1n;
  for (let i = 0; i < residues.length; i++) {
    const r = big(residues[i], "a residue");
    const n = pos(moduli[i]);
    const { g, x: u } = egcd(m, n);
    if (mod(r - x, g) !== 0n) return null;
    const l = (m / g) * n;
    x = mod(x + m * mod(((r - x) / g) * u, n / g), l);
    m = l;
  }
  return { x, modulus: m };
}

/** Solve a system of linear congruences aᵢ·x ≡ bᵢ (mod mᵢ): each reduced, then combined. { x, modulus } or null. */
export function solveLinearSystem(equations) {
  const rs = [];
  const ms = [];
  for (const [a, b, m] of equations) {
    const s = solveLinear(a, b, m);
    if (!s.solutions.length) return null;
    rs.push(s.reduced.x);
    ms.push(s.reduced.modulus);
  }
  return crt(rs, ms);
}

// ---- polynomial congruences ----

const coeffs = (poly) => poly.map((c) => big(c, "a coefficient"));
/** Evaluate an integer polynomial (constant term first) at x, mod m. */
export function polyMod(poly, x, m) {
  let r = 0n;
  for (let i = poly.length - 1; i >= 0; i--) r = mod(r * x + big(poly[i]), m);
  return r;
}
const deriv = (p) => p.slice(1).map((c, i) => c * BigInt(i + 1));

/**
 * Hensel's lemma: a root r of f mod p with f'(r) ≢ 0 (mod p) lifts to a unique
 * root mod p^k, by Newton's step r ← r − f(r)/f'(r). Returns the root mod p^k.
 * A root where f'(r) ≡ 0 is refused here; rootsMod handles those by trying
 * every lift.
 */
export function hensel(poly, r, p, k) {
  const f = coeffs(poly);
  const P = big(p, "p");
  const K = Number(big(k, "k"));
  let x = mod(big(r, "r"), P);
  if (polyMod(f, x, P) !== 0n) refuse("bad-argument", `${x} is not a root of the polynomial mod ${P}`);
  const df = deriv(f);
  const d = polyMod(df, x, P);
  if (d === 0n) refuse("singular", `f'(${x}) ≡ 0 (mod ${P}): the root does not lift uniquely`);
  let pk = P;
  for (let j = 1; j < K; j++) {
    pk *= P;
    x = mod(x - polyMod(f, x, pk) * inverse(polyMod(df, x, pk), pk), pk);
  }
  return x;
}

/** The roots of f mod a prime power p^k, including roots where f' vanishes (each lifted by trying every digit). */
function rootsModPrimePower(f, p, k) {
  if (p > 2000000n) refuse("too-large", `finding roots mod ${p} by search is limited to primes below 2·10⁶`);
  let roots = [];
  for (let x = 0n; x < p; x++) if (polyMod(f, x, p) === 0n) roots.push(x);
  let pk = p;
  for (let j = 1; j < k; j++) {
    const next = pk * p;
    const lifted = new Set();
    for (const r of roots) for (let t = 0n; t < p; t++) {
      const y = r + t * pk;
      if (polyMod(f, y, next) === 0n) lifted.add(y);
    }
    roots = [...lifted].sort((a, b) => (a < b ? -1 : 1));
    pk = next;
  }
  return roots;
}

/**
 * Every root of the integer polynomial f (coefficients, constant first) modulo
 * m: roots mod each prime power of m, combined by the Chinese remainder theorem.
 * Returns the sorted roots mod m.
 */
export function rootsMod(poly, m) {
  const f = coeffs(poly);
  const M = pos(m);
  if (M === 1n) return [0n];
  let acc = [{ x: 0n, modulus: 1n }];
  for (const { prime, exponent } of factor(M, { prove: false })) {
    const rs = rootsModPrimePower(f, prime, exponent);
    const pk = prime ** BigInt(exponent);
    const next = [];
    for (const a of acc) for (const r of rs) next.push(crt([a.x, r], [a.modulus, pk]));
    acc = next;
    if (!acc.length) return [];
  }
  return acc.map((a) => a.x).sort((a, b) => (a < b ? -1 : a > b ? 1 : 0));
}

// ---- primitive roots and indices ----

/** True when n has a primitive root: n = 1, 2, 4, p^k or 2p^k for an odd prime p. */
export function hasPrimitiveRoot(n) {
  const N = pos(n, "n");
  if (N <= 4n) return true;
  const f = factor(N, { prove: false });
  const odd = f.filter((x) => x.prime !== 2n);
  const two = f.find((x) => x.prime === 2n)?.exponent ?? 0;
  return odd.length === 1 && two <= 1;
}
/** The least primitive root modulo n (a generator of the units mod n), or null when there is none. */
export function primitiveRootMod(n) {
  const N = pos(n, "n");
  if (!hasPrimitiveRoot(N)) return null;
  if (N <= 2n) return N - 1n === 0n ? 0n : 1n;
  if (N === 4n) return 3n;
  const phi = totient(N);
  const qs = factor(phi, { prove: false }).map((x) => x.prime);
  for (let g = 2n; g < N; g++) if (gcd(g, N) === 1n && qs.every((q) => powmod(g, phi / q, N) !== 1n)) return g;
  return null;
}
/** The order of a modulo n: the least k > 0 with a^k ≡ 1. */
export function order(a, n) {
  return multiplicativeOrder(big(a, "a"), pos(n, "n"));
}
/**
 * The index (discrete logarithm) of h to the base g modulo n: the least
 * x ≥ 0 with g^x ≡ h, or null. Baby-step giant-step; n prime uses
 * algebraic's discreteLogMod.
 */
export function discreteLog(g, h, n) {
  const N = pos(n, "n");
  const G = mod(big(g, "g"), N);
  const H = mod(big(h, "h"), N);
  if (gcd(G, N) !== 1n) refuse("bad-argument", `${G} is not a unit mod ${N}`);
  const o = multiplicativeOrder(G, N);
  if (o > 1n << 44n) refuse("too-large", "the order is too large for baby-step giant-step here");
  const m = BigInt(Math.ceil(Math.sqrt(Number(o)))) + 1n;
  const table = new Map();
  let e = 1n;
  for (let j = 0n; j < m; j++) {
    if (!table.has(e)) table.set(e, j);
    e = (e * G) % N;
  }
  const step = inverse(powmod(G, m, N), N);
  let y = H;
  for (let i = 0n; i <= m; i++) {
    const j = table.get(y);
    if (j !== undefined) return mod(i * m + j, o);
    y = (y * step) % N;
  }
  return null;
}
export { discreteLogMod };

/** v_p(n), the exponent of p in n (Infinity for 0). */
export function valuation(n, p) {
  const P = big(p, "p");
  if (P < 2n) throw new RangeError("p must be at least 2");
  const x = big(n, "n");
  return val(x, P);
}
