// Arithmetic functions: divisors, σ, τ, μ, φ, λ, ω, Ω and their relatives,
// Dirichlet convolution and Möbius inversion, perfect, abundant and amicable
// numbers, and aliquot sequences. All exact, on BigInts.

import { big, gcd, lcm, refuse } from "./int.js";
import { factor } from "./primes.js";

const fs = (n, what = "n") => {
  const m = big(n, what);
  if (m < 1n) throw new RangeError(`${what} must be a positive whole number, not ${m}`);
  return factor(m, { prove: false });
};

/** The positive divisors of n, in increasing order. */
export function divisors(n) {
  let ds = [1n];
  for (const { prime, exponent } of fs(n)) {
    const next = [];
    for (const d of ds) {
      let pk = 1n;
      for (let k = 0; k <= exponent; k++, pk *= prime) next.push(d * pk);
    }
    ds = next;
  }
  return ds.sort((a, b) => (a < b ? -1 : a > b ? 1 : 0));
}
/** σ_k(n) = Σ_{d | n} d^k (k = 1 by default: the sum of divisors). */
export function sigma(n, k = 1) {
  const K = big(k, "k");
  if (K < 0n) throw new RangeError("k must be at least 0");
  let out = 1n;
  for (const { prime, exponent } of fs(n)) {
    if (K === 0n) out *= BigInt(exponent + 1);
    else {
      const q = prime ** K;
      out *= (q ** BigInt(exponent + 1) - 1n) / (q - 1n);
    }
  }
  return out;
}
/** τ(n), the number of divisors. */
export function tau(n) {
  return sigma(n, 0);
}
/** μ(n): 0 if a square divides n, else (−1)^(number of prime factors). */
export function mobius(n) {
  const f = fs(n);
  if (f.some((x) => x.exponent > 1)) return 0n;
  return f.length % 2 ? -1n : 1n;
}
/** Euler's φ(n): how many of 1…n are prime to n. */
export function totient(n) {
  let r = big(n, "n");
  for (const { prime } of fs(n)) r = (r / prime) * (prime - 1n);
  return r;
}
/** Jordan's totient J_k(n) = n^k Π (1 − p^−k). */
export function jordan(n, k) {
  const K = big(k, "k");
  let r = big(n, "n") ** K;
  for (const { prime } of fs(n)) r = (r / prime ** K) * (prime ** K - 1n);
  return r;
}
/** Carmichael's λ(n): the exponent of the group of units mod n, the least m with a^m ≡ 1 for every a prime to n. */
export function carmichael(n) {
  let r = 1n;
  for (const { prime, exponent } of fs(n)) {
    const e = BigInt(exponent);
    let l = prime ** (e - 1n) * (prime - 1n);
    if (prime === 2n && e >= 3n) l /= 2n;
    r = lcm(r, l);
  }
  return r;
}
/** ω(n): the number of distinct prime factors. */
export function omega(n) {
  return BigInt(fs(n).length);
}
/** Ω(n): the number of prime factors counted with multiplicity. */
export function bigOmega(n) {
  return fs(n).reduce((s, f) => s + BigInt(f.exponent), 0n);
}
/** Liouville's λ(n) = (−1)^Ω(n). */
export function liouville(n) {
  return bigOmega(n) % 2n ? -1n : 1n;
}
/** rad(n): the product of the distinct primes dividing n. */
export function radical(n) {
  return fs(n).reduce((s, f) => s * f.prime, 1n);
}
export function isSquarefree(n) {
  return fs(n).every((f) => f.exponent === 1);
}
/** s(n) = σ(n) − n, the sum of the proper divisors. */
export function aliquotSum(n) {
  return sigma(n) - big(n, "n");
}
/** "perfect", "abundant" or "deficient", by comparing σ(n) with 2n. */
export function classify(n) {
  const s = aliquotSum(n);
  const m = big(n, "n");
  return s === m ? "perfect" : s > m ? "abundant" : "deficient";
}
export const isPerfect = (n) => classify(n) === "perfect";
/**
 * The aliquot sequence n, s(n), s(s(n)), …, until it reaches 1 (then 0), enters
 * a cycle, or the step limit. Returns { terms, ending } with ending
 * "terminates", "perfect", "amicable", "cycle of length k" or "unfinished".
 */
export function aliquotSequence(n, { steps = 100, max = 10n ** 30n } = {}) {
  const terms = [big(n, "n")];
  const seen = new Map([[terms[0], 0]]);
  while (terms.length <= steps) {
    const x = terms[terms.length - 1];
    if (x === 0n) return { terms, ending: "terminates" };
    if (x > max) return { terms, ending: "unfinished" };
    const y = x === 1n ? 0n : aliquotSum(x);
    if (seen.has(y)) {
      const len = terms.length - seen.get(y);
      terms.push(y);
      return { terms, ending: len === 1 ? "perfect" : len === 2 ? "amicable" : `cycle of length ${len}` };
    }
    seen.set(y, terms.length);
    terms.push(y);
  }
  return { terms, ending: "unfinished" };
}
/** The amicable pairs [a, b], a < b ≤ limit (s(a) = b, s(b) = a), by sieving σ up to the limit (at most 10⁷). */
export function amicablePairs(limit) {
  const L = Number(big(limit, "limit"));
  if (L > 10000000) refuse("too-large", "the limit is 10⁷");
  const s = new Float64Array(L + 1);
  for (let d = 1; d <= L / 2; d++) for (let m = 2 * d; m <= L; m += d) s[m] += d;
  const out = [];
  for (let a = 2; a <= L; a++) {
    const b = s[a];
    if (b > a && b <= L && s[b] === a) out.push([BigInt(a), BigInt(b)]);
  }
  return out;
}
/**
 * The even perfect numbers 2^(p−1)(2^p − 1) for the first `count` Mersenne
 * primes 2^p − 1 (Euclid and Euler: these are all the even ones). Mersenne
 * primes are proven by the Lucas–Lehmer test.
 */
export function perfectNumbers(count) {
  const k = Number(big(count, "count"));
  if (k > 20) refuse("too-large", "ask for at most 20");
  const out = [];
  for (let p = 2; out.length < k; p++) {
    if (!isPrimeSmall(p)) continue;
    if (lucasLehmer(p)) out.push((1n << BigInt(p - 1)) * ((1n << BigInt(p)) - 1n));
  }
  return out;
}
const isPrimeSmall = (p) => {
  if (p < 2) return false;
  for (let d = 2; d * d <= p; d++) if (p % d === 0) return false;
  return true;
};
/** The Lucas–Lehmer test: for an odd prime p, 2^p − 1 is prime exactly when s_(p−2) ≡ 0, with s_0 = 4, s_(k+1) = s_k² − 2. */
export function lucasLehmer(p) {
  const P = Number(big(p, "p"));
  if (P === 2) return true;
  if (!isPrimeSmall(P)) return false;
  const M = (1n << BigInt(P)) - 1n;
  let s = 4n;
  for (let i = 0; i < P - 2; i++) s = (s * s - 2n) % M;
  return s === 0n;
}

// ---- Dirichlet convolution ----

/** The Dirichlet convolution f * g: n ↦ Σ_{d | n} f(d)·g(n/d). Functions take and return BigInts. */
export function dirichlet(f, g) {
  return (n) => {
    const m = big(n, "n");
    return divisors(m).reduce((s, d) => s + big(f(d)) * big(g(m / d)), 0n);
  };
}
/** Möbius inversion: given F(n) = Σ_{d | n} f(d), recover f(n) = Σ_{d | n} μ(d)·F(n/d). */
export function mobiusInversion(F) {
  return dirichlet(mobius, F);
}
/** The Dirichlet identity ε, the constant 1, and the identity function, for building convolutions. */
export const unit = (n) => (big(n) === 1n ? 1n : 0n);
export const one = () => 1n;
export const id = (n) => big(n);
/** True when f(mn) = f(m)f(n) for all coprime m, n ≤ limit (a check on examples, not a proof). */
export function isMultiplicative(f, limit = 60) {
  const L = BigInt(limit);
  for (let m = 1n; m <= L; m++) for (let n = 1n; n <= L; n++) if (gcd(m, n) === 1n && big(f(m * n)) !== big(f(m)) * big(f(n))) return false;
  return true;
}
