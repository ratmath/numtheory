// Primes: sieves, next and previous primes, counting, primality with proofs,
// and factoring by trial division, Pollard's rho and Pollard's p − 1.
//
// The house rule is that an answer is proven or says it is not. The strong
// Baillie–PSW test (algebraic's isProbablePrime) is known to be correct below
// 2^64, so below that `isPrime` is a proof. Above it, `isPrime` asks for a
// certificate (provePrime) and refuses to answer when none can be found,
// rather than return a probable answer as if it were certain.

import { meissel } from "./primecount.js";
import { isProbablePrime } from "@ratmath/algebraic";
import { abs, big, gcd, iroot, isqrt, mod, powmod, refuse } from "./int.js";
import { ecm } from "./ecm.js";
import { ecppStep } from "./ecpp.js";

export const PROVEN_BELOW = 1n << 64n;
const SIEVE_LIMIT = 100000000;

// ---- sieves ----

/** An odd-only sieve: bit i set when 2i + 1 is composite. */
function oddSieve(limit) {
  const n = Math.floor((limit - 1) / 2) + 1; // indices 0…n−1 stand for 1, 3, 5, …
  const comp = new Uint8Array(n);
  comp[0] = 1;
  for (let i = 1; ; i++) {
    const p = 2 * i + 1;
    if (p * p > limit) break;
    if (comp[i]) continue;
    for (let j = (p * p - 1) / 2; j < n; j += p) comp[j] = 1;
  }
  return comp;
}
const sieveSize = (x, what) => {
  const n = big(x, what);
  if (n > BigInt(SIEVE_LIMIT)) refuse("too-large", `${what} = ${n} is above the sieve's limit of ${SIEVE_LIMIT}`);
  return Number(n);
};

/** The primes up to a limit (at most 10⁷ here), as BigInts, by the sieve of Eratosthenes. */
export function primes(limit) {
  const L = sieveSize(limit, "limit");
  if (L > 10000000) refuse("too-large", "listing primes stops at 10⁷; primePi counts further");
  if (L < 2) return [];
  const comp = oddSieve(L);
  const out = [2n];
  for (let i = 1; i < comp.length; i++) if (!comp[i]) out.push(BigInt(2 * i + 1));
  return out;
}
/** π(x): how many primes are ≤ x. A sieve below 10⁷, Meissel's formula above, up to 10¹². */
export function primePi(x) {
  const n = big(x, "x");
  if (n < 2n) return 0n;
  if (n > 10000000n) return meissel(n).pi;
  const comp = oddSieve(Number(n));
  let c = 1;
  for (let i = 1; i < comp.length; i++) if (!comp[i]) c++;
  return BigInt(c);
}
/** The k-th prime (the 1st is 2), for k up to about 5·10⁶. */
export function nthPrime(k) {
  const m = big(k, "k");
  if (m < 1n) throw new RangeError("the first prime is nthPrime(1)");
  const K = Number(m);
  // p_k < k(ln k + ln ln k) for k ≥ 6 (Rosser)
  const bound = K < 6 ? 15 : Math.ceil(K * (Math.log(K) + Math.log(Math.log(K))));
  if (bound > SIEVE_LIMIT) refuse("too-large", `the ${K}th prime is beyond the sieve's reach`);
  if (K === 1) return 2n;
  const comp = oddSieve(bound);
  let c = 1;
  for (let i = 1; i < comp.length; i++) if (!comp[i] && ++c === K) return BigInt(2 * i + 1);
  throw new Error("unreachable: Rosser's bound");
}

// ---- primality ----

/**
 * True when n is prime. A proof below 2^64 (Baillie–PSW has no exceptions
 * there); above, a certificate is sought, and if none is found the call
 * throws with code "unproven" instead of guessing. Use isProbablePrime for a
 * fast answer that is not a proof.
 */
export function isPrime(n) {
  const m = big(n, "n");
  if (m < PROVEN_BELOW) return isProbablePrime(m);
  if (!isProbablePrime(m)) return false;
  const r = provePrime(m);
  if (r.proven) return r.prime;
  refuse("unproven", `${m} passes the Baillie–PSW test but no certificate was found (${r.reason}); isProbablePrime(n) gives the unproven answer`);
}
export { isProbablePrime };

/** The least prime > n. */
export function nextPrime(n) {
  let m = big(n, "n");
  if (m < 2n) return 2n;
  m += m % 2n === 0n ? 1n : 2n;
  while (!isProbablePrime(m)) m += 2n;
  if (m >= PROVEN_BELOW && !isPrime(m)) return nextPrime(m);
  return m;
}
/** The greatest prime < n (n > 2). */
export function prevPrime(n) {
  let m = big(n, "n");
  if (m <= 2n) refuse("bad-argument", "there is no prime below 2");
  if (m === 3n) return 2n;
  m -= m % 2n === 0n ? 1n : 2n;
  while (!isProbablePrime(m)) m -= 2n;
  if (m >= PROVEN_BELOW && !isPrime(m)) return prevPrime(m);
  return m;
}
/** The primes in [a, b], b < 2^64 and b − a ≤ 10⁷. */
export function primesBetween(a, b) {
  const lo = big(a, "a");
  const hi = big(b, "b");
  if (hi >= PROVEN_BELOW) refuse("too-large", "primesBetween works below 2^64, where every answer is proven");
  if (hi - lo > 10000000n) refuse("too-large", "the range is limited to 10⁷ numbers");
  const out = [];
  for (let m = lo < 2n ? 2n : lo; m <= hi; m += m === 2n ? 1n : m % 2n === 0n ? 1n : 2n) if (isProbablePrime(m)) out.push(m);
  return out;
}

// ---- factoring ----

/** A nontrivial factor of composite n by Pollard's rho (Brent's cycle finding), or null within the step budget. */
export function pollardRho(n, { steps = 1 << 20, seed = 1n } = {}) {
  n = big(n, "n");
  if (n % 2n === 0n) return 2n;
  for (let c = seed; c < seed + 20n; c++) {
    let y = 2n;
    let r = 1;
    let q = 1n;
    let g = 1n;
    let x = y;
    let ys = y;
    let used = 0;
    const f = (v) => (v * v + c) % n;
    while (g === 1n && used < steps) {
      x = y;
      for (let i = 0; i < r; i++) y = f(y);
      let k = 0;
      while (k < r && g === 1n) {
        ys = y;
        const m = Math.min(128, r - k);
        for (let i = 0; i < m; i++) {
          y = f(y);
          q = (q * abs(x - y)) % n;
        }
        g = gcd(q, n);
        k += m;
        used += m;
      }
      r *= 2;
    }
    if (g === n) {
      do {
        ys = f(ys);
        g = gcd(abs(x - ys), n);
      } while (g === 1n);
    }
    if (g !== 1n && g !== n) return g;
    if (used >= steps) return null;
  }
  return null;
}

/**
 * Pollard's p − 1: finds a prime factor p of n when p − 1 has only small
 * prime factors (all ≤ bound). Computes a^M with M the product of prime
 * powers up to the bound, then gcd(a^M − 1, n). Returns a factor or null.
 */
export function pollardPm1(n, bound = 100000) {
  n = big(n, "n");
  const B = Number(big(bound, "bound"));
  if (B > 10000000) refuse("too-large", "keep the bound at 10⁷ or below");
  let a = 2n;
  const ps = B >= 2 ? primes(B) : [];
  for (const p of ps) {
    let pk = p;
    while (pk * p <= BigInt(B)) pk *= p;
    a = powmod(a, pk, n);
  }
  const g = gcd(a - 1n, n);
  return g === 1n || g === n ? null : g;
}

/** Split off small prime factors (below `limit`); returns { found: Map(prime → exponent), rest }. */
export function trialDivide(n, limit = 100000) {
  let m = abs(big(n, "n"));
  const found = new Map();
  const add = (p) => found.set(p, (found.get(p) ?? 0) + 1);
  for (const p of [2n, 3n, 5n]) while (m % p === 0n) add(p), (m /= p);
  const L = BigInt(limit);
  // wheel of 30
  const steps = [4n, 2n, 4n, 2n, 4n, 6n, 2n, 6n];
  for (let p = 7n, i = 0; p < L && p * p <= m; p += steps[i++ & 7]) while (m % p === 0n) add(p), (m /= p);
  if (m > 1n && m < L * L) {
    // what is left has no factor below min(limit, √m): it is prime if m < limit²
    add(m);
    m = 1n;
  }
  return { found, rest: m };
}

/**
 * The prime factorisation of n ≠ 0 as [{ prime, exponent, proven }] in
 * increasing order (the sign is dropped). Trial division, then Pollard's rho,
 * p − 1 and the elliptic curve method with growing bounds. `proven` is true
 * for every factor below 2^64 and for larger factors with a certificate;
 * option prove: false skips the certificates, and effort (default 1) scales
 * the number of curves tried before giving up.
 */
export function factor(n, { prove = true, effort = 1 } = {}) {
  const m = big(n, "n");
  if (m === 0n) refuse("bad-argument", "0 has no factorisation");
  const { found, rest } = trialDivide(m);
  const stack = rest > 1n ? [rest] : [];
  while (stack.length) {
    const x = stack.pop();
    if (x === 1n) continue;
    if (isProbablePrime(x)) {
      found.set(x, (found.get(x) ?? 0) + 1);
      continue;
    }
    const r = iroot(x, 2);
    if (r * r === x) {
      stack.push(r, r);
      continue;
    }
    let d = pollardRho(x, { steps: 1 << 19 });
    if (!d) d = pollardPm1(x, 100000);
    for (const [B1, curves] of ECM_LEVELS) {
      if (d) break;
      d = ecm(x, { B1, curves: Math.ceil(curves * effort) });
    }
    if (!d) refuse("too-hard", `could not split ${x} with rho, p − 1 or the elliptic curve method; its prime factors are probably above 25 digits`);
    stack.push(d, x / d);
  }
  return [...found]
    .sort((a, b) => (a[0] < b[0] ? -1 : 1))
    .map(([prime, exponent]) => ({ prime, exponent, proven: prime < PROVEN_BELOW || (prove ? provePrime(prime).proven : false) }));
}
/** "2^3 · 3 · 5^2" for a factorisation. */
export function formatFactors(fs) {
  if (!fs.length) return "1";
  return fs.map(({ prime, exponent }) => (exponent === 1 ? `${prime}` : `${prime}^${exponent}`)).join(" · ");
}

const ECM_LEVELS = [[2000, 30], [11000, 90], [50000, 120]];

// ---- certificates ----

const SMALL = 1n << 32n;

/** Factor m as far as the budget allows: { primes: [[q, e]], rest } with every q passing BPSW. */
function partialFactor(m, budget) {
  const { found, rest } = trialDivide(m, 100000);
  const done = [...found];
  const stack = rest > 1n ? [rest] : [];
  let left = rest > 1n ? [] : [];
  while (stack.length) {
    const x = stack.pop();
    if (isProbablePrime(x)) {
      const k = done.findIndex(([q]) => q === x);
      if (k >= 0) done[k][1]++;
      else done.push([x, 1]);
      continue;
    }
    const d = pollardRho(x, { steps: budget }) ?? pollardPm1(x, 20000);
    if (d) stack.push(d, x / d);
    else left.push(x);
  }
  return { primes: done, rest: left.reduce((a, b) => a * b, 1n) };
}

/**
 * Try to prove n prime. Returns one of
 *   { n, prime: true, proven: true, certificate }   a certificate verifyCertificate checks
 *   { n, prime: false, proven: true, factor }       or with `witness` (a Miller–Rabin witness)
 *   { n, prime: true, proven: false, reason }       probably prime, but n − 1 could not be factored far enough
 * The certificate is Pocklington's (n − 1 factored past √n) or
 * Brillhart–Lehmer–Selfridge's (past the cube root), with each prime factor
 * of n − 1 certified in turn; small primes are certified by trial division.
 * When n − 1 cannot be factored far enough, an elliptic curve certificate is
 * sought instead (Atkin and Morain): a curve mod n with a point whose order
 * is a prime q large enough to rule out any factor of n, and q is certified
 * in turn. Options: budget (Pollard rho steps per factor of n − 1),
 * elliptic: false to use n − 1 only, maxDiscriminant and maxClassNumber for
 * the curves tried.
 */
export function provePrime(n, { budget = 1 << 18, elliptic = true, maxDiscriminant = 5000, maxClassNumber = 12 } = {}) {
  const N = big(n, "n");
  if (N < 2n) return { n: N, prime: false, proven: true, reason: "below 2" };
  const td = trialDivide(N, 1000);
  if (td.found.size && !(td.found.size === 1 && td.found.has(N))) return { n: N, prime: false, proven: true, factor: [...td.found.keys()][0] };
  if (N < SMALL) {
    const p = isProbablePrime(N);
    if (!p) return { n: N, prime: false, proven: true, factor: smallFactor(N) };
    return { n: N, prime: true, proven: true, certificate: { n: N, kind: "small" } };
  }
  if (!isProbablePrime(N)) {
    for (let a = 2n; a < 200n; a++) if (millerRabinWitness(N, a)) return { n: N, prime: false, proven: true, witness: a };
    const d = pollardRho(N, { steps: budget });
    if (d) return { n: N, prime: false, proven: true, factor: d };
    return { n: N, prime: false, proven: false, reason: "fails Baillie–PSW, but no Miller–Rabin witness below 200 and no factor was found" };
  }
  const state = elliptic ? { steps: 4000, backtrack: 4, rhoBudget: 1 << 12, maxDiscriminant, maxClassNumber } : null;
  const cert = certify(N, budget, 0, state);
  if (cert.certificate) return { n: N, prime: true, proven: true, certificate: cert.certificate };
  if (cert.witness) return { n: N, prime: false, proven: true, witness: cert.witness };
  if (cert.factor) return { n: N, prime: false, proven: true, factor: cert.factor };
  return { n: N, prime: true, proven: false, reason: cert.reason };
}

function smallFactor(n) {
  for (let p = 2n; p * p <= n; p += p === 2n ? 1n : 2n) if (n % p === 0n) return p;
  return n;
}
/** True when a proves n composite by the strong (Miller–Rabin) test. */
export function millerRabinWitness(n, a) {
  n = big(n, "n");
  a = big(a, "a");
  if (n % 2n === 0n) return n !== 2n;
  let d = n - 1n;
  let s = 0;
  while (d % 2n === 0n) {
    d /= 2n;
    s++;
  }
  let x = powmod(a, d, n);
  if (x === 1n || x === n - 1n || x === 0n) return false;
  for (let i = 1; i < s; i++) {
    x = (x * x) % n;
    if (x === n - 1n) return false;
  }
  return true;
}

/** A certificate for the probable prime n: by n − 1 first, then by an elliptic curve when `state` allows it. */
function certify(n, budget, depth, state) {
  if (n < SMALL) return isProbablePrime(n) ? { certificate: { n, kind: "small" } } : { reason: `${n} is not prime` };
  if (depth > 1500) return { reason: "the chain of certificates is too deep" };
  const r = certifyMinusOne(n, budget, depth, state);
  if (r.certificate || r.witness || r.factor || !state) return r;
  // inside an elliptic chain the n − 1 attempts get a small budget: most of them fail
  const e = ecppStep(n, { certify: (q) => certify(q, Math.min(budget, 1 << 12), depth + 1, state), partialFactor, state });
  if (e.certificate || e.factor) return e;
  return { reason: `${r.reason}; ${e.reason}` };
}
function certifyMinusOne(n, budget, depth, state) {
  const m = n - 1n;
  const { primes: qs } = partialFactor(m, budget);
  qs.sort((a, b) => (a[0] < b[0] ? -1 : 1));
  const usable = [];
  let F = 1n;
  for (const [q, e] of qs) {
    const sub = certify(q, budget, depth + 1, state);
    if (!sub.certificate) continue;
    // a base for q: a^(n−1) ≡ 1 and gcd(a^((n−1)/q) − 1, n) = 1
    let base = null;
    for (let a = 2n; a < 2000n; a++) {
      if (powmod(a, m, n) !== 1n) return { witness: a };
      if (gcd(powmod(a, m / q, n) - 1n, n) === 1n) {
        base = a;
        break;
      }
    }
    if (base === null) continue;
    usable.push({ q, e, a: base, certificate: sub.certificate });
    F *= q ** BigInt(e);
    if ((F + 1n) * (F + 1n) > n) break;
  }
  const R = m / F;
  if ((F + 1n) * (F + 1n) > n) return { certificate: { n, kind: "pocklington", F, R, factors: usable } };
  if (F * F * F > n) {
    const c2 = R / F;
    const c1 = R % F;
    const disc = c1 * c1 - 4n * c2;
    if (!(disc >= 0n && isqrt(disc) ** 2n === disc)) return { certificate: { n, kind: "bls", F, R, c1, c2, factors: usable } };
    // then n = (aF + 1)(bF + 1) with a, b the roots of x² − c1·x + c2
    const a = (c1 + isqrt(disc)) / 2n;
    return { factor: a * F + 1n };
  }
  return { reason: `only ${F.toString(2).length} of the ${n.toString(2).length} bits of n − 1 could be factored and certified (a third are needed)` };
}

export { mod };
