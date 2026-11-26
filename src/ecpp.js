// The search for elliptic curve primality certificates (Atkin and Morain).
//
// To prove n prime, find a curve E mod n whose number of points m would be
// k·q with q a prime above (n^(1/4) + 1)², and a point P with m·P = 0 and
// (m/q)·P ≠ 0. Then n is prime (Goldwasser and Kilian), provided q is, and q
// is proved the same way: a chain down to numbers small enough to check
// directly. verify.js checks such a chain without using anything here.
//
// Counting points on random curves is too slow, so curves are built with a
// known count by complex multiplication: when 4n = t² + |D|v² for a negative
// discriminant D, a curve whose j-invariant is a root mod n of the Hilbert
// class polynomial H_D has n + 1 ± t points. H_D is computed here from the
// values of j at the reduced forms of discriminant D, in fixed-point complex
// arithmetic, and rounded to integers.
//
// Nothing in this file needs to be trusted for a certificate to be sound: a
// wrong class polynomial or a failed square root only means that no
// certificate is found.

import { abs, gcd, isqrt, isSquare, mod, powmod } from "./int.js";
import { ecAdd, ecMul } from "./ecmod.js";
import { invMod, NotInvertible, oneRoot } from "./polymod.js";

// ---- reduced forms and class polynomials ----

/** The reduced primitive forms [a, b, c] of a negative discriminant D (as numbers). */
export function reducedTriples(D) {
  const out = [];
  const A = Math.floor(Math.sqrt(-D / 3));
  for (let a = 1; a <= A; a++) {
    for (let b = -a + 1; b <= a; b++) {
      if ((b - D) % 2 !== 0) continue;
      const num = b * b - D;
      if (num % (4 * a) !== 0) continue;
      const c = num / (4 * a);
      if (c < a || (c === a && b < 0)) continue;
      if (g3(a, b, c) !== 1) continue;
      out.push([a, b, c]);
    }
  }
  return out;
}
const g2 = (a, b) => {
  a = Math.abs(a);
  b = Math.abs(b);
  while (b) [a, b] = [b, a % b];
  return a;
};
const g3 = (a, b, c) => g2(g2(a, b), c);

function isFundamental(D) {
  const m = -D;
  let core;
  if (m % 4 === 3) core = m;
  else if (m % 4 === 0 && ((m / 4) % 4 === 1 || (m / 4) % 4 === 2)) core = m / 4;
  else return false;
  for (let d = 2; d * d <= core; d++) if (core % (d * d) === 0) return false;
  return true;
}

// fixed-point complex numbers: [re, im] scaled by 2^P
function fixedOps(P) {
  const one = 1n << P;
  // products are truncated towards zero, so that small terms of either sign die out
  const fm = (a, b) => {
    const x = a * b;
    return x < 0n ? -(-x >> P) : x >> P;
  };
  const cmul = ([a, b], [c, d]) => [fm(a, c) - fm(b, d), fm(a, d) + fm(b, c)];
  const cinv = ([a, b]) => {
    const nn = a * a + b * b;
    return [(a << (2n * P)) / nn, (-b << (2n * P)) / nn];
  };
  const atanInv = (k) => {
    let sum = 0n;
    let term = one / k;
    const k2 = k * k;
    for (let i = 0n; term !== 0n; i++) {
      sum += i % 2n === 0n ? term / (2n * i + 1n) : -(term / (2n * i + 1n));
      term /= k2;
    }
    return sum;
  };
  const pi = 16n * atanInv(5n) - 4n * atanInv(239n);
  const cexp = (z) => {
    const big = abs(z[0]) > abs(z[1]) ? abs(z[0]) : abs(z[1]);
    const s = BigInt(Math.max(0, big.toString(2).length - Number(P) + 6));
    const w = [z[0] / (1n << s), z[1] / (1n << s)];
    let sum = [one, 0n];
    let term = [one, 0n];
    for (let k = 1n; term[0] !== 0n || term[1] !== 0n; k++) {
      term = cmul(term, w);
      term = [term[0] / k, term[1] / k];
      sum = [sum[0] + term[0], sum[1] + term[1]];
    }
    for (let i = 0n; i < s; i++) sum = cmul(sum, sum);
    return sum;
  };
  return { one, fm, cmul, cinv, pi, cexp };
}

const classPolyCache = new Map();
/**
 * The Hilbert class polynomial H_D for a negative discriminant D: the monic
 * integer polynomial whose roots are the j-invariants of the curves with
 * complex multiplication by the order of discriminant D. Its degree is the
 * class number h(D). Coefficients are BigInts, constant term first; null
 * when the rounding to integers does not pass its check.
 */
export function classPolynomial(D) {
  D = Number(D);
  if (!(D < 0) || !Number.isInteger(D) || ((D % 4) + 4) % 4 > 1) throw new RangeError("the discriminant must be negative and 0 or 1 mod 4");
  if (D < -1000000) throw new RangeError("class polynomials are limited to |D| ≤ 10⁶");
  if (classPolyCache.has(D)) return classPolyCache.get(D);
  const forms = reducedTriples(D);
  const sumInv = forms.reduce((s, [a]) => s + 1 / a, 0);
  const P = BigInt(Math.ceil(((Math.PI * Math.sqrt(-D)) / Math.LN2) * (sumInv + 2)) + 160);
  const { one, cmul, cinv, pi, cexp } = fixedOps(P);
  const sqrtD = isqrt(BigInt(-D) << (2n * P));
  let poly = [[one, 0n]];
  for (const [a, b] of forms) {
    // q = exp(2πiτ), τ = (−b + i√|D|)/(2a)
    const q = cexp([-(pi * sqrtD) / (BigInt(a) << P), -(pi * BigInt(b)) / BigInt(a)]);
    // j = E4³/Δ with E4 = 1 + 240 Σ n³qⁿ/(1 − qⁿ), Δ = q·Π(1 − qⁿ)²⁴
    let e4 = [one, 0n];
    let prod = [one, 0n];
    let qn = q;
    for (let k = 1n; qn[0] !== 0n || qn[1] !== 0n; k++) {
      const om = [one - qn[0], -qn[1]];
      const t = cmul(qn, cinv(om));
      const c = 240n * k * k * k;
      e4 = [e4[0] + c * t[0], e4[1] + c * t[1]];
      prod = cmul(prod, om);
      qn = cmul(qn, q);
    }
    const p2 = cmul(prod, prod);
    const p4 = cmul(p2, p2);
    const p8 = cmul(p4, p4);
    const p24 = cmul(cmul(p8, p8), p8);
    const j = cmul(cmul(e4, cmul(e4, e4)), cinv(cmul(q, p24)));
    // poly ← poly·(X − j)
    const next = new Array(poly.length + 1);
    for (let i = 0; i <= poly.length; i++) {
      const up = i > 0 ? poly[i - 1] : [0n, 0n];
      const here = i < poly.length ? cmul(poly[i], j) : [0n, 0n];
      next[i] = [up[0] - here[0], up[1] - here[1]];
    }
    poly = next;
  }
  const half = one >> 1n;
  const tol = one >> 24n;
  const out = [];
  let ok = true;
  for (const [re, im] of poly) {
    const r = (re + half) >> P;
    if (abs(re - (r << P)) > tol || abs(im) > tol) ok = false;
    out.push(r);
  }
  const result = ok ? out : null;
  classPolyCache.set(D, result);
  return result;
}

// ---- arithmetic mod a probable prime ----

function jacobiSymbol(a, n) {
  a = mod(a, n);
  let t = 1;
  while (a !== 0n) {
    while (a % 2n === 0n) {
      a /= 2n;
      const r = n % 8n;
      if (r === 3n || r === 5n) t = -t;
    }
    [a, n] = [n, a];
    if (a % 4n === 3n && n % 4n === 3n) t = -t;
    a %= n;
  }
  return n === 1n ? t : 0;
}
/** A square root of a mod n when n is prime; null when none is found (n may be composite). */
function sqrtModN(a, n) {
  a = mod(a, n);
  if (a === 0n) return 0n;
  if (jacobiSymbol(a, n) !== 1) return null;
  let r;
  if (n % 4n === 3n) r = powmod(a, (n + 1n) / 4n, n);
  else {
    let q = n - 1n;
    let s = 0;
    while (q % 2n === 0n) {
      q /= 2n;
      s++;
    }
    let z = 2n;
    while (z < 1000n && jacobiSymbol(z, n) !== -1) z++;
    let c = powmod(z, q, n);
    let t = powmod(a, q, n);
    r = powmod(a, (q + 1n) / 2n, n);
    let m = s;
    while (t !== 1n) {
      let i = 0;
      let tt = t;
      while (tt !== 1n && i < m) {
        tt = (tt * tt) % n;
        i++;
      }
      if (i === 0 || i >= m) return null;
      const b = powmod(c, 1n << BigInt(m - i - 1), n);
      r = (r * b) % n;
      c = (b * b) % n;
      t = (t * c) % n;
      m = i;
    }
  }
  return (r * r) % n === a ? r : null;
}

/** [t, v] with 4n = t² + |D|v² (Cornacchia), or null. */
export function cornacchia4(n, D) {
  const d = BigInt(-D);
  if (jacobiSymbol(BigInt(D), n) !== 1) return null;
  let x0 = sqrtModN(BigInt(D), n);
  if (x0 === null) return null;
  if (mod(x0 - BigInt(D), 2n) !== 0n) x0 = n - x0;
  let [a, b] = [2n * n, x0];
  const L = isqrt(4n * n);
  while (b > L) [a, b] = [b, a % b];
  const rest = 4n * n - b * b;
  if (rest % d !== 0n) return null;
  const c = rest / d;
  if (!isSquare(c)) return null;
  return [b, isqrt(c)];
}

// ---- the search ----

const discriminants = [];
let scannedTo = 2;
/** Fundamental discriminants −3, −4, −7, … with class number at most hMax, down to −limit. */
function discriminantList(limit, hMax) {
  for (let m = scannedTo + 1; m <= limit; m++) {
    if (!isFundamental(-m)) continue;
    const h = reducedTriples(-m).length;
    discriminants.push({ D: -m, h });
  }
  if (limit > scannedTo) scannedTo = limit;
  return discriminants.filter((d) => d.D >= -limit && d.h <= hMax);
}

/** Curves [a, b] mod n with complex multiplication by D, both (or all) twists. */
function cmCurves(n, D) {
  let g = 2n;
  while (g < 1000n && jacobiSymbol(g, n) !== -1) g++;
  const out = [];
  if (D === -3) for (let b = 1n; b <= 40n; b++) out.push([0n, b]);
  else if (D === -4) for (let a = 1n; a <= 30n; a++) out.push([a, 0n]);
  else {
    const H = classPolynomial(D);
    if (!H) return out;
    const j = oneRoot(H.map((c) => mod(c, n)), n);
    if (j === null || j === 0n || j === 1728n % n) return out;
    const c = (j * invMod(1728n - j, n)) % n;
    const [a, b] = [(3n * c) % n, (2n * c) % n];
    out.push([a, b], [(a * g * g) % n, (((b * g * g) % n) * g) % n]);
  }
  return out;
}

/**
 * One step of the chain for the probable prime n. Returns
 *   { certificate }  an elliptic certificate, its prime q already certified;
 *   { factor }       n turned out composite;
 *   { reason }       nothing found within the limits.
 * `tools` supplies certify(q) (which may come back here) and
 * partialFactor(m) from primes.js, and a shared `state` with a step budget.
 */
export function ecppStep(n, tools) {
  const { certify, partialFactor, state } = tools;
  if (gcd(n, 6n) !== 1n) return { reason: "n is divisible by 2 or 3" };
  const tried = new Set();
  let dead = 0;
  try {
    for (const { D } of discriminantList(state.maxDiscriminant, state.maxClassNumber)) {
      if (state.steps <= 0) return { reason: "the search budget for elliptic curve certificates ran out" };
      const tv = cornacchia4(n, D);
      if (!tv) continue;
      const [t, v] = tv;
      const traces = [t, -t];
      if (D === -4) traces.push(2n * v, -2n * v);
      if (D === -3) traces.push((t + 3n * v) / 2n, -(t + 3n * v) / 2n, (t - 3n * v) / 2n, -(t - 3n * v) / 2n);
      let curves = null;
      for (const tr of traces) {
        const m = n + 1n - tr;
        if (tried.has(m)) continue;
        tried.add(m);
        const { primes } = partialFactor(m, state.rhoBudget);
        const q = primes.reduce((best, [p]) => (p > best ? p : best), 0n);
        if (q === m || q < 5n || (isqrt(q) - 1n) ** 4n <= n) continue;
        state.steps--;
        curves ??= cmCurves(n, D);
        const hit = findPoint(n, curves, m, q);
        if (!hit) continue;
        const sub = certify(q);
        if (sub.certificate) {
          return { certificate: { n, kind: "elliptic", discriminant: D, a: hit.a, b: hit.b, x: hit.x, y: hit.y, m, q, certificate: sub.certificate } };
        }
        if (++dead >= state.backtrack) return { reason: `no certificate for the primes q reached from ${n}` };
      }
    }
  } catch (e) {
    if (e instanceof NotInvertible && e.factor > 1n && e.factor < n) return { factor: e.factor };
    throw e;
  }
  return { reason: `no suitable curve with |D| ≤ ${state.maxDiscriminant} and class number ≤ ${state.maxClassNumber}` };
}

/** A curve among `curves` and a point P on it with m·P = 0 and (m/q)·P ≠ 0. */
function findPoint(n, curves, m, q) {
  const k = m / q;
  for (const [a, b] of curves) {
    if (gcd(mod(4n * a * a * a + 27n * b * b, n), n) !== 1n) continue;
    let found = 0;
    for (let x = 1n; x < 200n && found < 3; x++) {
      const rhs = mod(x * x * x + a * x + b, n);
      const y = sqrtModN(rhs, n);
      if (y === null || y === 0n) continue;
      found++;
      const U = ecMul(k, [x, y], a, n);
      if (U === null) continue;
      if (ecMul(q, U, a, n) !== null) break; // this curve has another order
      return { a, b, x, y };
    }
  }
  return null;
}
