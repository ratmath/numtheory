// Quadratic residues and reciprocity: the Legendre, Jacobi and Kronecker
// symbols, the reciprocity computation written out step by step, Euler's
// criterion, Gauss's lemma, and square roots modulo any n.

import { sqrtMod } from "@ratmath/algebraic";
import { big, gcd, mod, powmod, refuse } from "./int.js";
import { isPrime, factor } from "./primes.js";
import { crt, hensel, rootsMod } from "./congruence.js";

const oddPos = (n) => {
  const N = big(n, "n");
  if (N <= 0n || N % 2n === 0n) throw new RangeError(`the Jacobi symbol needs an odd positive n, not ${N}`);
  return N;
};

/** The Jacobi symbol (a/n) for odd n > 0, by reciprocity: −1, 0 or 1 (as a number). */
export function jacobi(a, n) {
  let N = oddPos(n);
  let A = mod(big(a, "a"), N);
  let t = 1;
  while (A !== 0n) {
    while (A % 2n === 0n) {
      A /= 2n;
      const r = N % 8n;
      if (r === 3n || r === 5n) t = -t;
    }
    [A, N] = [N, A];
    if (A % 4n === 3n && N % 4n === 3n) t = -t;
    A %= N;
  }
  return N === 1n ? t : 0;
}
/** The Legendre symbol (a/p) for an odd prime p. Refuses a p that is not prime. */
export function legendre(a, p) {
  const P = big(p, "p");
  if (P === 2n || !isPrime(P)) throw new RangeError(`the Legendre symbol needs an odd prime, and ${P} is not one`);
  return jacobi(a, P);
}
/** The Kronecker symbol (a/n) for any whole n, extending Jacobi's to even and negative n. */
export function kronecker(a, n) {
  let A = big(a, "a");
  let N = big(n, "n");
  if (N === 0n) return A === 1n || A === -1n ? 1 : 0;
  let t = 1;
  if (N < 0n) {
    N = -N;
    if (A < 0n) t = -t;
  }
  while (N % 2n === 0n) {
    N /= 2n;
    if (A % 2n === 0n) return 0;
    const r = mod(A, 8n);
    if (r === 3n || r === 5n) t = -t;
  }
  return N === 1n ? t : t * jacobi(A, N) || 0;
}

/**
 * The Jacobi symbol computed with every step written out: reduce the top mod
 * the bottom, take out factors of 2 by the second supplement, flip by
 * quadratic reciprocity. Returns { value, steps } with each step a line of text.
 */
export function jacobiSteps(a, n) {
  let N = oddPos(n);
  let A = big(a, "a");
  let sign = 1;
  const steps = [];
  const sym = (x, y) => `${sign < 0 ? "−" : ""}(${x}/${y})`;
  steps.push(`${sym(A, N)}`);
  for (;;) {
    if (A < 0n || A >= N) {
      A = mod(A, N);
      steps.push(`= ${sym(A, N)}   reduce the top mod ${N}`);
    }
    if (A === 0n) {
      steps.push(`= 0   since ${N === 1n ? "…" : `gcd > 1`}`);
      return { value: N === 1n ? sign : 0, steps };
    }
    if (A === 1n) {
      steps.push(`= ${sign < 0 ? "−1" : "1"}   since (1/${N}) = 1`);
      return { value: sign, steps };
    }
    if (A % 2n === 0n) {
      let k = 0;
      while (A % 2n === 0n) {
        A /= 2n;
        k++;
      }
      const r = N % 8n;
      const two = r === 1n || r === 7n ? 1 : -1;
      if (k % 2 && two < 0) sign = -sign;
      steps.push(`= ${sym(A, N)}   take out 2${k > 1 ? `^${k}` : ""}: (2/${N}) = ${two < 0 ? "−1" : "1"} since ${N} ≡ ${r} (mod 8)${k % 2 ? "" : ", and an even power counts as 1"}`);
      continue;
    }
    if (gcd(A, N) !== 1n) {
      steps.push(`= 0   since ${A} and ${N} share the factor ${gcd(A, N)}`);
      return { value: 0, steps };
    }
    const flip = A % 4n === 3n && N % 4n === 3n;
    if (flip) sign = -sign;
    [A, N] = [N, A];
    steps.push(`= ${sym(A, N)}   reciprocity: ${flip ? `both are 3 mod 4, so the sign changes` : `one of them is 1 mod 4, so the sign stays`}`);
  }
}

/** Euler's criterion: a^((p−1)/2) mod p is 1 for a residue, p − 1 for a non-residue, 0 when p | a. */
export function eulerCriterion(a, p) {
  const P = big(p, "p");
  if (P === 2n || !isPrime(P)) throw new RangeError("Euler's criterion needs an odd prime");
  const power = powmod(big(a, "a"), (P - 1n) / 2n, P);
  return { power, value: power === P - 1n ? -1 : Number(power) };
}
/**
 * Gauss's lemma: (a/p) = (−1)^μ where μ counts the least positive residues of
 * a, 2a, …, ((p−1)/2)·a that exceed p/2. Returns { value, mu, large }.
 */
export function gaussLemma(a, p) {
  const P = big(p, "p");
  if (P === 2n || !isPrime(P)) throw new RangeError("Gauss's lemma needs an odd prime");
  const A = big(a, "a");
  if (A % P === 0n) refuse("bad-argument", "p divides a");
  const large = [];
  for (let k = 1n; k <= (P - 1n) / 2n; k++) {
    const r = mod(k * A, P);
    if (2n * r > P) large.push(r);
  }
  return { value: large.length % 2 ? -1 : 1, mu: large.length, large };
}
/** The quadratic residues mod n: the distinct squares x² mod n, sorted (n up to 10⁶). */
export function quadraticResidues(n) {
  const N = big(n, "n");
  if (N < 1n || N > 1000000n) refuse("too-large", "n must be between 1 and 10⁶");
  const s = new Set();
  for (let x = 0n; x <= N / 2n; x++) s.add((x * x) % N);
  return [...s].sort((a, b) => (a < b ? -1 : 1));
}

/** Every x mod n with x² ≡ a (mod n), sorted. */
export function sqrtsMod(a, n) {
  const N = big(n, "n");
  if (N < 1n) throw new RangeError("the modulus must be at least 1");
  const A = mod(big(a, "a"), N);
  let acc = [{ x: 0n, modulus: 1n }];
  for (const { prime, exponent } of factor(N, { prove: false })) {
    const pk = prime ** BigInt(exponent);
    let rs;
    if (prime < 2000000n || pk === prime || A % prime !== 0n) {
      if (prime >= 2000000n && A % prime !== 0n) {
        const r = sqrtMod(A % prime, prime);
        rs = r === null ? [] : [...new Set([r, prime - r].map((x) => hensel([-A, 0n, 1n], x, prime, exponent)))];
      } else if (prime >= 2000000n) rs = A % pk === 0n ? [0n] : [];
      else rs = rootsMod([-A, 0n, 1n], pk);
    } else refuse("too-large", `square roots of multiples of the large prime ${prime} modulo its powers are not handled`);
    const next = [];
    for (const s of acc) for (const r of rs) next.push(crt([s.x, r], [s.modulus, pk]));
    acc = next;
    if (!acc.length) return [];
  }
  return acc.map((s) => s.x).sort((a, b) => (a < b ? -1 : a > b ? 1 : 0));
}
export { sqrtMod };
