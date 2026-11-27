// π(x) exactly, far beyond what a sieve can list: Meissel's formula
//   π(x) = φ(x, a) + a − 1 − P₂(x, a),   a = π(x^(1/3)),
// where φ(x, a) counts the n ≤ x with no prime factor among the first a
// primes (Legendre's sum, computed by its recurrence), and P₂ counts the
// n ≤ x with exactly two prime factors, both above the a-th prime.

import { refuse, big } from "./int.js";

const LIMIT = 10n ** 12n;

/** Odd-only sieve up to n with fast π(m) for any m ≤ n. */
function countingSieve(n) {
  const size = Math.floor((n + 1) / 2); // index i stands for 2i + 1
  const words = new Uint32Array((size >> 5) + 1);
  // bit set = composite (or 1)
  words[0] |= 1;
  for (let i = 1; ; i++) {
    const p = 2 * i + 1;
    if (p * p > n) break;
    if (words[i >> 5] & (1 << (i & 31))) continue;
    for (let j = (p * p - 1) / 2; j < size; j += p) words[j >> 5] |= 1 << (j & 31);
  }
  // prefix counts of primes per word
  const prefix = new Uint32Array(words.length + 1);
  const pop = (v) => {
    v = v - ((v >>> 1) & 0x55555555);
    v = (v & 0x33333333) + ((v >>> 2) & 0x33333333);
    return (((v + (v >>> 4)) & 0x0f0f0f0f) * 0x01010101) >>> 24;
  };
  for (let w = 0; w < words.length; w++) prefix[w + 1] = prefix[w] + 32 - pop(words[w]);
  const pi = (m) => {
    if (m < 2) return 0;
    if (m > n) throw new RangeError("π lookup beyond the sieve");
    const i = (m - 1) >> 1; // the last odd number ≤ m has index i
    const w = i >> 5;
    const b = i & 31;
    const mask = b === 31 ? 0xffffffff : (1 << (b + 1)) - 1;
    return 1 + prefix[w] + (b + 1) - pop((words[w] & mask) >>> 0); // +1 for the prime 2
  };
  const isPrime = (m) => m === 2 || (m > 2 && m % 2 === 1 && !(words[(m >> 1) >> 5] & (1 << ((m >> 1) & 31))));
  return { pi, isPrime, size };
}

/**
 * Legendre's φ(x, a): how many n in 1 … x have no prime factor among the
 * first a primes. φ(x, 0) = x and φ(x, a) = φ(x, a − 1) − φ(⌊x/p_a⌋, a − 1).
 */
export function legendrePhi(x, a) {
  const X = Number(big(x, "x"));
  const A = Number(big(a, "a"));
  if (X < 0 || A < 0) throw new RangeError("legendrePhi needs x ≥ 0 and a ≥ 0");
  if (X > 1e9) refuse("too-large", "legendrePhi on its own stops at 10⁹; primePi goes further");
  const ps = [];
  const sv = countingSieve(Math.max(100, Math.ceil(A * (Math.log(A + 2) + Math.log(Math.log(A + 3)) + 2))));
  for (let m = 2; ps.length < A; m++) if (sv.isPrime(m)) ps.push(m);
  const rec = (n, b) => {
    if (b === 0) return n;
    if (n < ps[b - 1]) return n >= 1 ? 1 : 0;
    let s = n;
    for (let i = 1; i <= b; i++) {
      const m = Math.floor(n / ps[i - 1]);
      if (m === 0) break;
      s -= rec(m, i - 1);
    }
    return s;
  };
  return BigInt(rec(X, A));
}

/**
 * π(x) by Meissel's formula, with its pieces: { pi, y, a, phi, P2 } where
 * y = ⌊x^(1/3)⌋ and a = π(y). For x up to 10¹².
 */
export function meissel(x) {
  const n = big(x, "x");
  if (n > LIMIT) refuse("too-large", `primePi stops at ${LIMIT}`);
  const X = Number(n);
  if (X < 2) return { pi: 0n, y: 0n, a: 0n, phi: 0n, P2: 0n };
  let y = Math.floor(Math.cbrt(X));
  while ((y + 1) ** 3 <= X) y++;
  while (y ** 3 > X) y--;
  let sq = Math.floor(Math.sqrt(X));
  while ((sq + 1) * (sq + 1) <= X) sq++;
  while (sq * sq > X) sq--;
  const top = Math.max(Math.floor(X / Math.max(y, 1)), sq, 30); // x^(2/3): the largest argument π is asked for
  const sv = countingSieve(top);
  const ps = [0]; // 1-based: ps[i] is the i-th prime
  for (let m = 2; m <= sq; m++) if (sv.isPrime(m)) ps.push(m);
  const a = sv.pi(y);
  // P2 = Σ_{y < p ≤ √x} (π(x/p) − π(p) + 1)
  let P2 = 0;
  for (let i = a + 1; i < ps.length; i++) P2 += sv.pi(Math.floor(X / ps[i])) - i + 1;
  // φ(x, a), with the first C primes handled by a table
  const C = Math.min(a, 6);
  let prim = 1;
  for (let i = 1; i <= C; i++) prim *= ps[i];
  const tbl = new Int32Array(prim + 1);
  {
    const coprime = new Uint8Array(prim + 1).fill(1);
    for (let i = 1; i <= C; i++) for (let m = 0; m <= prim; m += ps[i]) coprime[m] = 0;
    for (let m = 1; m <= prim; m++) tbl[m] = tbl[m - 1] + coprime[m];
  }
  const totient = tbl[prim];
  const small = (m) => Math.floor(m / prim) * totient + tbl[m % prim];
  const phi = (m, b) => {
    if (b <= C) {
      if (b === C) return small(m);
      // fewer primes than the table uses: plain recurrence
      let s = m;
      for (let i = 1; i <= b; i++) s -= phi(Math.floor(m / ps[i]), i - 1);
      return s;
    }
    const pb = ps[b];
    if (m < pb) return m >= 1 ? 1 : 0;
    if (m < pb * pb && m <= top) return sv.pi(m) - b + 1;
    let s = small(m);
    for (let i = C + 1; i <= b; i++) {
      const p = ps[i];
      const t = Math.floor(m / p);
      if (t < p) {
        // every remaining term is 1 while p_i ≤ m
        const last = Math.min(b, sv.pi(m));
        if (last >= i) s -= last - i + 1;
        break;
      }
      s -= phi(t, i - 1);
    }
    return s;
  };
  const ph = phi(X, a);
  return { pi: BigInt(ph + a - 1 - P2), y: BigInt(y), a: BigInt(a), phi: BigInt(ph), P2: BigInt(P2) };
}
