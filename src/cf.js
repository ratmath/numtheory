// Continued fractions of exact numbers.
//
// For a rational: Euclid's algorithm. For any real algebraic number x, the
// root of p in (lo, hi): t = ⌊x⌋ is read off once the interval sits between
// two integers, then x = t + 1/y where y is the root of yᵈ·p(t + 1/y) in the
// image interval. Every term is exact, as far as you ask.
//
// Quadratic irrationals are eventually periodic (Lagrange). Written as
// (P + √D)/Q with Q | D − P², the next complete quotient is
// (P' + √D)/Q' with P' = aQ − P and Q' = (D − P'²)/Q, so the period is found
// exactly when a pair (P, Q) repeats. Pell's equation and the fundamental
// unit of a real quadratic field come out of that period.

//
// Moved here from @ratmath/algebraic (roadmap section 4.3): nothing in
// algebraic uses continued fractions, so they live in number theory.

import { AlgebraicReal, NumberField, fromRational, simplify, poly } from "@ratmath/algebraic";
import { reduce, sqrtFloor, abs, gcd } from "@ratmath/algebraic/q";
import { Rational } from "@ratmath/core";

const { primitive, degree } = poly;

// ---- three helpers from algebraic's radical.js, which is not exported ----

/** [s, d] with D = s²·d and d square-free (sign kept on d); trial division to 10⁵. */
function squarefreeSplit(D) {
  if (D === 0n) return [0n, 0n];
  const sg = D < 0n ? -1n : 1n;
  let n = abs(D);
  let s = 1n;
  let d = 1n;
  for (let p = 2n; p <= 100000n && p * p <= n; p += p === 2n ? 1n : 2n) {
    while (n % (p * p) === 0n) {
      n /= p * p;
      s *= p;
    }
    if (n % p === 0n) {
      n /= p;
      d *= p;
    }
  }
  const r = sqrtFloor(n);
  if (r * r === n) s *= r;
  else d *= n;
  return [s, sg * d];
}
const cmpPair = ([a, b], [c, d]) => {
  const x = a * d - c * b;
  return x > 0n ? 1 : x < 0n ? -1 : 0;
};
/** The root of c + b·x + a·x² on the given side, as (p + q√d)/r in lowest terms. */
function quadraticParts([c, b, a], side) {
  const D = b * b - 4n * a * c;
  const [s, d] = squarefreeSplit(D);
  let p = -b;
  let q = BigInt(side) * s;
  let r = 2n * a;
  if (r < 0n) {
    p = -p;
    q = -q;
    r = -r;
  }
  const g = gcd(gcd(p, q), r) || 1n;
  return { p: p / g, q: q / g, d, r: r / g };
}
/** Which root of a real quadratic the interval holds: +1 above the vertex (for a > 0), else −1. */
function realSide([, b, a], get, refine) {
  const v = reduce(-b, 2n * a);
  for (let bits = 8; bits <= 4096; bits *= 2) {
    const [lo, hi] = get();
    const up = cmpPair(lo, v) > 0;
    const down = cmpPair(hi, v) < 0;
    if (up || down) return (up ? 1 : -1) * (a > 0n ? 1 : -1);
    refine(bits);
  }
  throw new Error("could not place the root beside the vertex");
}

const floorDiv = (a, b) => {
  const q = a / b;
  return (a % b !== 0n) && ((a < 0n) !== (b < 0n)) ? q - 1n : q;
};
const floorPair = ([n, d]) => floorDiv(n, d);

/** Taylor shift: the coefficients of p(x + t). */
function shift(p, t) {
  const a = [...p];
  const n = a.length - 1;
  for (let i = 0; i < n; i += 1) for (let j = n - 1; j >= i; j -= 1) a[j] += t * a[j + 1];
  return a;
}

/** Strict bound on the positive roots of p (Cauchy). */
function bound(p) {
  const lc = abs(p[p.length - 1]);
  let m = 0n;
  for (let i = 0; i < p.length - 1; i += 1) if (abs(p[i]) > m) m = abs(p[i]);
  return m / lc + 2n;
}

const asReal = (x) => (x instanceof AlgebraicReal ? x : fromRational(x));

/**
 * The first `terms` partial quotients [a₀; a₁, a₂, …] of a real number
 * (rational or AlgebraicReal) as BigInts. A rational's expansion ends early.
 */
export function continuedFraction(x, terms = 20) {
  const a = asReal(x);
  const r = a.rational;
  if (r) {
    let [n, d] = r;
    const out = [];
    while (d !== 0n && out.length < terms) {
      const t = floorDiv(n, d);
      out.push(t);
      [n, d] = [d, n - t * d];
    }
    return out;
  }
  let cur = new AlgebraicReal(a.poly, a.lo, a.hi);
  const out = [];
  while (out.length < terms) {
    let t;
    for (let bits = 4; ; bits += 4) {
      cur.refine(bits);
      const f = floorPair(cur.lo);
      if (floorPair(cur.hi) === f || (cur.hi[1] === 1n && cur.hi[0] === f + 1n)) {
        t = f;
        break;
      }
    }
    out.push(t);
    // y = 1/(x − t), a root of yᵈ·p(t + 1/y)
    const q = primitive([...shift(cur.poly, t)].reverse());
    const lo = reduce(cur.lo[0] - t * cur.lo[1], cur.lo[1]);
    const hi = reduce(cur.hi[0] - t * cur.hi[1], cur.hi[1]);
    const newLo = reduce(hi[1], hi[0]);
    const newHi = lo[0] === 0n ? [bound(q), 1n] : reduce(lo[1], lo[0]);
    cur = new AlgebraicReal(q, newLo, newHi);
  }
  return out;
}

/** The convergents p/q of a list of partial quotients, as Core Rationals. */
export function convergents(terms) {
  let [p0, q0, p1, q1] = [1n, 0n, 0n, 1n];
  return terms.map((t) => {
    const a = BigInt(t);
    [p0, q0, p1, q1] = [a * p0 + p1, a * q0 + q1, p0, q0];
    return new Rational(p0, q0);
  });
}

/** The real number [a₀; a₁, …, aₙ] (finite) or [a₀; …, (period repeating)] when `period` is given. */
export function fromContinuedFraction(terms, period = null) {
  if (!period || !period.length) {
    const c = convergents(terms);
    return fromRational(c[c.length - 1]);
  }
  // y = [period, y] is the positive root of q·y² + (q' − p)·y − p' = 0, from the period's convergents
  const c = matrix(period);
  const [p, pp, q, qq] = c;
  const poly = primitive([-pp, qq - p, q]);
  let y = quadraticRoot(poly, 1);
  if (!(y.lo[0] > 0n)) y = quadraticRoot(poly, -1);
  // x = [terms…, y]
  const [a, b, cc, d] = matrix(terms);
  // x = (a·y + b)/(c·y + d)
  return mobius(y, a, b, cc, d);
}

function matrix(terms) {
  let [p, pp, q, qq] = [1n, 0n, 0n, 1n];
  for (const t of terms) [p, pp, q, qq] = [BigInt(t) * p + pp, p, BigInt(t) * q + qq, q];
  return [p, pp, q, qq];
}

function quadraticRoot(poly, side) {
  const { p, q, d, r } = quadraticParts(poly, side);
  return fromQuadratic(p, q, d, r);
}

/** (p + q√d)/r as an AlgebraicReal on its minimal polynomial. */
export function fromQuadratic(p, q, d, r) {
  p = BigInt(p); q = BigInt(q); d = BigInt(d); r = BigInt(r);
  if (q === 0n || d === 0n) return fromRational([p, r]);
  // root of (r x − p)² = q² d: r² x² − 2 p r x + p² − q² d
  const poly = primitive([p * p - q * q * d, -2n * p * r, r * r]);
  if (degree(poly) < 2) return fromRational(reduce(p, r));
  const s = sqrtFloor(q * q * d);
  const sg = q > 0n ? 1n : -1n;
  // q√d = sg·√(q²d), bracketed by s and s + 1, then divided by r
  const ends = [reduce(p + sg * s, r), reduce(p + sg * (s + 1n), r)];
  ends.sort((u, v) => (u[0] * v[1] < v[0] * u[1] ? -1 : 1));
  if (s * s === q * q * d) return fromRational(reduce(p + sg * s, r));
  return simplify(new AlgebraicReal(poly, ends[0], ends[1]));
}

function mobius(y, a, b, c, d) {
  // (a y + b)/(c y + d) for a quadratic y; stays quadratic
  const { p, q, d: D, r } = partsOf(y);
  // y = (p + q√D)/r
  const num = [a * p + b * r, a * q]; // (num0 + num1 √D)/r
  const den = [c * p + d * r, c * q];
  // divide: (n0 + n1√D)/(d0 + d1√D) = (n0 + n1√D)(d0 − d1√D)/(d0² − d1² D)
  const N0 = num[0] * den[0] - num[1] * den[1] * D;
  const N1 = num[1] * den[0] - num[0] * den[1];
  const M = den[0] * den[0] - den[1] * den[1] * D;
  return fromQuadratic(N0, N1, D, M);
}

function partsOf(y) {
  const side = realSide(y.poly, () => [y.lo, y.hi], (b) => y.refine(b));
  return quadraticParts(y.poly, side);
}

/** (P + √D)/Q with Q | D − P², for a quadratic irrational x. */
function pqd(x) {
  const a = asReal(x);
  if (a.rational || degree(simplify(a).poly) !== 2) throw new RangeError("needs a quadratic irrational");
  const s = simplify(a);
  const { p, q, d, r } = partsOf(s);
  let P = q > 0n ? p : -p;
  let Q = q > 0n ? r : -r;
  let D = q * q * d;
  if ((D - P * P) % Q !== 0n) {
    const k = abs(Q);
    P *= k;
    D *= Q * Q;
    Q *= k;
  }
  return { P, Q, D };
}

function* quadTerms({ P, Q, D }) {
  const s = sqrtFloor(D);
  for (;;) {
    const a = Q > 0n ? floorDiv(P + s, Q) : floorDiv(P + s + 1n, Q);
    yield { a, P, Q };
    P = a * Q - P;
    Q = (D - P * P) / Q;
  }
}

/**
 * The eventually periodic expansion of a quadratic irrational:
 * { preperiod: [a₀, …], period: [aₖ, …] } (BigInts), so x = [preperiod; period, period, …].
 */
export function periodicContinuedFraction(x) {
  const seen = new Map();
  const terms = [];
  for (const { a, P, Q } of quadTerms(pqd(x))) {
    const key = `${P},${Q}`;
    if (seen.has(key)) {
      const k = seen.get(key);
      return { preperiod: terms.slice(0, k), period: terms.slice(k) };
    }
    seen.set(key, terms.length);
    terms.push(a);
  }
  return null;
}

/**
 * The smallest positive solution of Pell's equation x² − D·y² = 1 (D a
 * positive non-square integer), and of x² − D·y² = −1 when it has one.
 */
export function pell(D) {
  D = BigInt(D);
  const s = sqrtFloor(D);
  if (D < 2n || s * s === D) throw new RangeError("pell needs a positive non-square D");
  const { preperiod, period } = periodicContinuedFraction(fromQuadratic(0n, 1n, D, 1n));
  const L = period.length;
  const terms = [...preperiod, ...period, ...period].slice(0, 2 * L);
  const c = convergents(terms);
  const at = (k) => ({ x: c[k].numerator, y: c[k].denominator });
  if (L % 2 === 0) return { x: at(L - 1).x, y: at(L - 1).y, negative: null };
  return { ...at(2 * L - 1), negative: at(L - 1) };
}

/**
 * The fundamental unit ε > 1 of the ring of integers of Q(√d) (d > 1
 * square-free), as an element of NumberField.quadratic(d), with its norm ±1.
 * It is the product of the complete quotients over one period of the reduced
 * number (⌊(√d − 1)/2⌋-shifted ω).
 */
export function fundamentalUnit(d) {
  d = BigInt(d);
  const [s] = squarefreeSplit(d);
  if (d < 2n || s !== 1n) throw new RangeError("fundamentalUnit needs a square-free d > 1");
  const K = NumberField.quadratic(d);
  // ω = (1 + √d)/2 when d ≡ 1 mod 4, else √d; θ = ω + k is reduced (θ > 1, −1 < θ̄ < 0)
  const one4 = d % 4n === 1n;
  const r = sqrtFloor(d);
  const omegaParts = one4 ? { p: 1n, r: 2n } : { p: 0n, r: 1n };
  // θ̄ = (p − √d)/r + k ∈ (−1, 0): k = ⌊(√d − p)/r⌋ … choose k so that θ̄ < 0 < θ̄ + 1
  const k = floorDiv(r - omegaParts.p, omegaParts.r);
  const theta = { P: omegaParts.p + k * omegaParts.r, Q: omegaParts.r, D: d };
  // normalise to Q | D − P²
  let { P, Q, D } = theta;
  if ((D - P * P) % Q !== 0n) {
    P *= Q;
    D *= Q * Q;
    Q *= Q;
  }
  const seen = new Set();
  let eps = K.one;
  const sq = Number(sqrtFloor(D / d)); // D = d·m²
  const m = BigInt(sq);
  for (const st of quadTerms({ P, Q, D })) {
    const key = `${st.P},${st.Q}`;
    if (seen.has(key)) break;
    seen.add(key);
    // complete quotient (P + m√d)/Q
    eps = eps.mul(K.of([[st.P, st.Q], [m, st.Q]].map(([n, dd]) => reduce(n, dd))));
  }
  if (eps.toAlgebraic().toNumber() < 1) eps = eps.inv();
  return { unit: eps, norm: Number(eps.norm().numerator) };
}
