// Values of the Riemann zeta function, the Hurwitz zeta function and
// Dirichlet L-functions: exact ones where a closed form exists (Bernoulli
// numbers at the non-positive integers, rational multiples of powers of π at
// the positive even ones), and proven enclosures anywhere else on the line or
// in the plane, by Euler–Maclaurin summation with an explicit remainder in
// ball arithmetic.

import { Rational, RationalInterval } from "@ratmath/core";
import { NumberField, cyclotomic } from "@ratmath/algebraic";
import { Ctx } from "./ball.js";
import { DirichletCharacter } from "./characters.js";
import { refuse, big, gcd } from "./int.js";

// ---- Bernoulli numbers ----

const BERN = [new Rational(1n)];
/** The Bernoulli number B_n, with B_1 = −1/2. */
export function bernoulli(n) {
  n = Number(big(n, "n"));
  if (n < 0) throw new RangeError("bernoulli needs n ≥ 0");
  for (let m = BERN.length; m <= n; m++) {
    if (m > 1 && m % 2 === 1) {
      BERN.push(new Rational(0n));
      continue;
    }
    let s = new Rational(0n);
    let c = 1n; // C(m+1, k)
    for (let k = 0; k < m; k++) {
      s = s.add(BERN[k].multiply(new Rational(c)));
      c = (c * BigInt(m + 1 - k)) / BigInt(k + 1);
    }
    BERN.push(s.divide(new Rational(BigInt(-(m + 1)))));
  }
  return BERN[n];
}
/** The Bernoulli polynomial B_n(x) at a rational x. */
export function bernoulliPolynomial(n, x) {
  n = Number(big(n, "n"));
  x = toRat(x);
  let s = new Rational(0n);
  let c = 1n;
  for (let k = 0; k <= n; k++) {
    // C(n, k)·B_k·x^(n−k)
    s = s.add(bernoulli(k).multiply(new Rational(c)).multiply(x.pow(n - k)));
    c = (c * BigInt(n - k)) / BigInt(k + 1);
  }
  return s;
}
function toRat(x) {
  if (x instanceof Rational) return x;
  if (typeof x === "bigint") return new Rational(x);
  if (typeof x === "number") {
    if (!Number.isInteger(x)) throw new TypeError(`write ${x} as an exact fraction, for example "1/2"`);
    return new Rational(BigInt(x));
  }
  if (typeof x === "string") return new Rational(x.trim());
  if (x && typeof x.numerator === "bigint") return new Rational(x.numerator, x.denominator);
  throw new TypeError(`expected an exact number, not ${String(x)}`);
}

// ---- exact special values ----

const cycFields = new Map();
function cyclotomicField(N) {
  if (!cycFields.has(N)) cycFields.set(N, new NumberField(cyclotomic(N), { name: `ζ${N}`, check: false }));
  return cycFields.get(N);
}
/**
 * The generalised Bernoulli number B_{k,χ} = q^(k−1)·Σ χ(a)·B_k(a/q): a
 * Rational for a real character, otherwise an element of the cyclotomic
 * field of the character's order (ζ_N = e(1/N)).
 */
export function generalizedBernoulli(k, chi) {
  k = Number(big(k, "k"));
  const q = chi.modulus;
  const N = chi.order;
  const coef = Array.from({ length: N }, () => new Rational(0n));
  for (let a = 1; a <= q; a++) {
    const e = chi._k(a);
    if (e === null) continue;
    coef[e] = coef[e].add(bernoulliPolynomial(k, new Rational(BigInt(a), BigInt(q))));
  }
  const scale = new Rational(BigInt(q)).pow(k - 1);
  const cs = coef.map((c) => c.multiply(scale));
  if (N === 1) return cs[0];
  if (N === 2) return cs[0].subtract(cs[1]);
  return cyclotomicField(N).of(cs.map((c) => [c.numerator, c.denominator]));
}
/**
 * L(n, χ) exactly at an integer n ≤ 0: −B_{1−n,χ}/(1−n). A Rational for a
 * real character, an element of a cyclotomic field otherwise.
 */
export function lValueExact(chi, n) {
  n = Number(big(n, "n"));
  if (n > 0) refuse("bad-argument", "exact L-values are given at the integers n ≤ 0; use lValue for an enclosure elsewhere");
  const k = 1 - n;
  const B = generalizedBernoulli(k, chi);
  const f = new Rational(-1n, BigInt(k));
  return B instanceof Rational ? B.multiply(f) : B.mul(B.field.of([[f.numerator, f.denominator]]));
}
/**
 * ζ(n) in closed form: a Rational for n ≤ 0 (−B_{1−n}/(1−n)), and for even
 * n > 0 the object { coefficient, piPower } standing for coefficient·π^n.
 */
export function zetaExact(n) {
  n = Number(big(n, "n"));
  if (n === 1) refuse("pole", "ζ has its pole at 1");
  if (n <= 0) return n === 0 ? new Rational(-1n, 2n) : bernoulli(1 - n).divide(new Rational(BigInt(-(1 - n))));
  if (n % 2) refuse("unsupported", `no closed form is known for ζ(${n}); zetaValue gives an enclosure`);
  // ζ(2k) = (−1)^(k+1)·B_2k·2^(2k−1)·π^(2k)/(2k)!
  let fact = 1n;
  for (let i = 2n; i <= BigInt(n); i++) fact *= i;
  const c = bernoulli(n)
    .multiply(new Rational(2n ** BigInt(n - 1), fact))
    .multiply(new Rational((n / 2) % 2 ? 1n : -1n));
  return { coefficient: c, piPower: n, toString: () => `${c.numerator === 1n ? "" : c.numerator}π${n === 1 ? "" : `^${n}`}${c.denominator === 1n ? "" : `/${c.denominator}`}` };
}

// ---- enclosures ----

function parseS(s) {
  if (Array.isArray(s)) return { re: toRat(s[0]), im: toRat(s[1]) };
  if (s && typeof s === "object" && !(s instanceof Rational) && "re" in s) return { re: toRat(s.re), im: toRat(s.im ?? 0) };
  return { re: toRat(s), im: new Rational(0n) };
}
const ratNum = (q) => Number(q.numerator) / Number(q.denominator);

/** ζ(s, u/v) as a complex ball; with regular = true the pole part 1/(s − 1) is left out (and s = 1 is allowed). */
function hurwitzBall(ctx, sx, u, v, regular) {
  const P = ctx.bits;
  const s = [ctx.frac(sx.re.numerator, sx.re.denominator), ctx.frac(sx.im.numerator, sx.im.denominator)];
  const sigma = ratNum(sx.re);
  const size = Math.abs(sigma) + Math.abs(ratNum(sx.im));
  const M = Math.ceil(P / 4) + 2 + Math.max(0, Math.ceil((2 - sigma) / 2));
  const N = Math.ceil(0.7 * (size + 2 * M)) + 2;
  const atOne = sx.im.numerator === 0n && sx.re.numerator === sx.re.denominator;
  if (atOne && !regular) refuse("pole", "the function has a pole at s = 1");
  const negS = [ctx.neg(s[0]), ctx.neg(s[1])];
  const powNegS = (lnb) => ctx.cexp(ctx.cmul(negS, ctx.creal(lnb)));
  let sum = [ctx.int(0), ctx.int(0)];
  for (let n = 0; n < N; n++) sum = ctx.cadd(sum, powNegS(ctx.ln(ctx.frac(BigInt(n) * v + u, v))));
  const X = ctx.frac(BigInt(N) * v + u, v);
  const lnX = ctx.ln(X);
  const Xs = powNegS(lnX); // X^(−s)
  // X^(1−s)/(s − 1), or its regular part
  const sm1 = [ctx.sub(s[0], ctx.int(1)), s[1]];
  if (atOne) sum = ctx.cadd(sum, ctx.creal(ctx.neg(lnX)));
  else {
    let top = ctx.cmul(ctx.creal(X), Xs);
    if (regular) top = ctx.csub(top, ctx.creal(ctx.int(1)));
    sum = ctx.cadd(sum, ctx.cdiv(top, sm1));
  }
  sum = ctx.cadd(sum, [ctx.mul(Xs[0], ctx.frac(1n, 2n)), ctx.mul(Xs[1], ctx.frac(1n, 2n))]);
  // Σ B_2j/(2j)! · (s)_(2j−1) · X^(−s−2j+1)
  const invX = ctx.inv(X);
  const invX2 = ctx.mul(invX, invX);
  // one running term T = (s)_(2j−1)·X^(−s−2j+1), kept moderate in size: fixed point cannot carry a huge
  // rising factorial times a tiny power separately
  let T = ctx.cmul(ctx.cmul(s, Xs), ctx.creal(invX));
  let fact = 2n; // (2j)!
  for (let j = 1; j <= M; j++) {
    const B = bernoulli(2 * j);
    const c = ctx.frac(B.numerator, B.denominator * fact);
    sum = ctx.cadd(sum, ctx.cmul(T, ctx.creal(c)));
    const s1 = [ctx.add(s[0], ctx.int(2 * j - 1)), s[1]];
    const s2 = [ctx.add(s[0], ctx.int(2 * j)), s[1]];
    T = ctx.cmul(ctx.cmul(ctx.cmul(T, s1), s2), ctx.creal(invX2));
    fact *= BigInt(2 * j + 1) * BigInt(2 * j + 2);
  }
  // remainder: |R| ≤ 4·|(s)_2M|/(2π)^2M · X^(−(σ+2M−1))/(σ+2M−1), and σ + 2M − 1 ≥ 1 by the choice of M
  // here T = (s)_(2M+1)·X^(−s−2M−1); use |(s)_2M| ≤ |(s)_(2M+1)| (as |s + 2M| ≥ 1) and (2π)² > 39
  const up = (a, b) => ((a * b) >> ctx.p) + 1n;
  let bound = ctx.cmag(T);
  const X2 = ctx.mag(ctx.mul(X, X));
  bound = up(bound, X2); // ·X² → X^(−σ−2M+1)
  // 1/|s + 2M| ≤ 1/(σ + 2M) ≤ 1 since σ + 2M ≥ 2
  bound = 4n * bound;
  for (let j = 0; j < M; j++) bound = bound / 39n + 1n;
  bound += 1n;
  return [
    [sum[0][0], sum[0][1] + bound],
    [sum[1][0], sum[1][1] + bound],
  ];
}
function chiBall(ctx, k, N) {
  const g = Number(gcd(BigInt(k), BigInt(N)));
  const [a, b] = [k / g, N / g];
  if (b === 1) return [ctx.int(1), ctx.int(0)];
  if (b === 2) return [ctx.int(-1), ctx.int(0)];
  if (b === 4) return [ctx.int(0), ctx.int(a === 1 ? 1 : -1)];
  const ang = ctx.mul(ctx.pi(), ctx.frac(2n * BigInt(a), BigInt(b)));
  const [c, s] = ctx.cosSin(ang);
  return [c, s];
}
function lBall(ctx, chi, sx) {
  const q = chi.modulus;
  if (q === 1) return hurwitzBall(ctx, sx, 1n, 1n, false);
  const atOne = sx.im.numerator === 0n && sx.re.numerator === sx.re.denominator;
  if (chi.isPrincipal && atOne) refuse("pole", "the L-function of a principal character has a pole at s = 1");
  // away from 1 the pole parts are kept; at 1 they cancel because Σ χ(a) = 0
  let sum = [ctx.int(0), ctx.int(0)];
  for (let a = 1; a <= q; a++) {
    const k = chi._k(a);
    if (k === null) continue;
    sum = ctx.cadd(sum, ctx.cmul(chiBall(ctx, k, chi.order), hurwitzBall(ctx, sx, BigInt(a), BigInt(q), atOne)));
  }
  // q^(−s)
  const s = [ctx.frac(sx.re.numerator, sx.re.denominator), ctx.frac(sx.im.numerator, sx.im.denominator)];
  const qs = ctx.cexp(ctx.cmul([ctx.neg(s[0]), ctx.neg(s[1])], ctx.creal(ctx.ln(ctx.int(q)))));
  return ctx.cmul(sum, qs);
}
function enclose(fn, bits, real) {
  bits = Number(bits);
  if (!(bits >= 8 && bits <= 4000)) throw new RangeError("bits must be between 8 and 4000");
  for (const extra of [48, 112, 240, 500]) {
    const ctx = new Ctx(bits + extra);
    const z = fn(ctx);
    const worst = z[0][1] > z[1][1] ? z[0][1] : z[1][1];
    if (worst >> BigInt(extra - 1) > 0n) continue;
    const sh = BigInt(extra - 2);
    const D = 1n << BigInt(bits + 2);
    const iv = ([m, r]) => new RationalInterval(new Rational((m - r) >> sh, D), new Rational(-((-(m + r)) >> sh), D));
    if (real) {
      if (z[1][0] - z[1][1] > 0n || z[1][0] + z[1][1] < 0n) throw new Error("a real value came out with a non-zero imaginary part");
      return iv(z[0]);
    }
    return { re: iv(z[0]), im: iv(z[1]) };
  }
  refuse("precision", "the enclosure did not reach the requested width; the argument may be too close to a pole or too large");
}

/**
 * A proven enclosure of ζ(s): a RationalInterval for real s, { re, im } for
 * s = [σ, t] or { re, im } in the plane. s is given exactly (whole numbers,
 * fraction text, Rationals); option bits (default 64) sets the width 2^−bits.
 */
export function zetaValue(s, { bits = 64 } = {}) {
  const sx = parseS(s);
  return enclose((ctx) => hurwitzBall(ctx, sx, 1n, 1n, false), bits, sx.im.numerator === 0n);
}
/** A proven enclosure of the Hurwitz zeta function ζ(s, a) = Σ (n + a)^−s for a rational a in (0, 1]. */
export function hurwitzZeta(s, a, { bits = 64 } = {}) {
  const sx = parseS(s);
  const A = toRat(a);
  if (A.numerator <= 0n || A.numerator > A.denominator) throw new RangeError("hurwitzZeta takes a rational a with 0 < a ≤ 1");
  return enclose((ctx) => hurwitzBall(ctx, sx, A.numerator, A.denominator, false), bits, sx.im.numerator === 0n);
}
/**
 * A proven enclosure of L(s, χ) = Σ χ(n)·n^−s, continued to the whole plane:
 * a RationalInterval when s is real and χ is a real character, { re, im }
 * otherwise. At s = 1 the character must not be principal.
 */
export function lValue(chi, s, { bits = 64 } = {}) {
  if (!(chi instanceof DirichletCharacter)) throw new TypeError("lValue takes a DirichletCharacter first");
  if (chi.modulus > 5000) refuse("too-large", "L-values are computed for moduli up to 5000");
  const sx = parseS(s);
  return enclose((ctx) => lBall(ctx, chi, sx), bits, sx.im.numerator === 0n && chi.isReal);
}
