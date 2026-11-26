// p-adic numbers, computed lazily: as many digits as are asked for.
//
// A p-adic number x is held as a rule that, given N, returns a rational r
// with v_p(x − r) ≥ N ("x to absolute precision N"), together with a lower
// bound on v_p(x). Rationals are exact; sums, products, quotients, square
// roots, Hensel roots, exp and log are new rules built on the old ones, and
// each asks its inputs for exactly the precision it needs. So nothing is
// lost to a fixed precision chosen in advance: ask for more digits and the
// whole computation is redone to that depth. Results are cached.
//
// What cannot be known lazily is reported, not guessed: the valuation of a
// number that may be zero is found by computing until a nonzero digit shows,
// and gives up (code "unknown-valuation") after a limit.

import { Rational } from "@ratmath/core";
import { big, inverse, mod, rat, refuse, val } from "./int.js";
import { sqrtMod } from "@ratmath/algebraic";

const VALUATION_SEARCH = 2000;

/** v_p of a Rational (Infinity for 0). */
function vq(r, p) {
  if (r.numerator === 0n) return Infinity;
  return val(r.numerator, p) - val(r.denominator, p);
}
/** The residue of a Rational with v_p ≥ 0 modulo p^k, as a BigInt in [0, p^k). */
function residue(r, p, k) {
  const pk = p ** BigInt(k);
  if (k <= 0) return 0n;
  return mod(r.numerator * inverse(mod(r.denominator, pk), pk), pk);
}

export class Padic {
  #rule;
  #cache = new Map();
  #val; // undefined: not yet known
  /**
   * Low level: p, a rule N ↦ Rational with v_p(x − rule(N)) ≥ N, and a lower
   * bound on the valuation. Use padic(x, p) to make one.
   */
  constructor(p, rule, vmin, exact = null) {
    this.p = p;
    this.#rule = rule;
    this.vmin = vmin;
    this.exact = exact; // a Rational when the number is rational and known exactly
    if (exact) this.#val = vq(exact, p);
    Object.freeze(this);
  }
  /** A rational r with v_p(x − r) ≥ N. */
  approx(N) {
    if (this.exact) return this.exact;
    N = Math.max(N, this.vmin);
    if (!this.#cache.has(N)) {
      // a cached higher precision serves as well
      let best = null;
      for (const [k, r] of this.#cache) if (k >= N && (best === null || k < best[0])) best = [k, r];
      this.#cache.set(N, best ? best[1] : this.#rule(N));
    }
    return this.#cache.get(N);
  }
  /** The valuation v_p(x); throws "unknown-valuation" if x agrees with 0 to the search limit. */
  valuation(limit = VALUATION_SEARCH) {
    if (this.#val !== undefined) return this.#val;
    for (let N = this.vmin + 8; ; N += Math.max(8, N - this.vmin)) {
      const r = this.approx(N);
      const v = vq(r, this.p);
      if (v < N) {
        this.#val = v;
        return v;
      }
      if (N - this.vmin > limit) refuse("unknown-valuation", `this ${this.p}-adic number is 0 to ${N} digits; it may be exactly 0, which a computation by digits cannot show`);
    }
  }
  /** |x|_p = p^(−v) as a Rational (0 for 0). */
  norm() {
    if (this.exact && this.exact.numerator === 0n) return new Rational(0n);
    const v = this.valuation();
    return v >= 0 ? new Rational(1n, this.p ** BigInt(v)) : new Rational(this.p ** BigInt(-v));
  }
  isZero() {
    return !!this.exact && this.exact.numerator === 0n;
  }
  /** x mod p^N as an integer in [0, p^N), for x with v_p(x) ≥ 0. */
  mod(N) {
    const k = Number(N);
    const r = this.approx(k);
    if (vq(r, this.p) < 0) refuse("bad-argument", "this number is not a p-adic integer (it has a negative valuation)");
    return residue(r, this.p, k);
  }
  /**
   * The digits d_v, d_(v+1), … of x = Σ d_i p^i, starting at the valuation v:
   * { valuation, digits } with `count` digits (each 0 ≤ d < p).
   */
  digits(count = 20) {
    const k = Number(count);
    if (this.isZero()) return { valuation: Infinity, digits: [] };
    const v = this.valuation();
    const r = this.approx(v + k);
    const u = r.multiply(v >= 0 ? new Rational(1n, this.p ** BigInt(v)) : new Rational(this.p ** BigInt(-v)));
    let m = residue(u, this.p, k);
    const out = [];
    for (let i = 0; i < k; i++) {
      out.push(m % this.p);
      m /= this.p;
    }
    return { valuation: v, digits: out };
  }
  /** The usual right-to-left writing: …d₂d₁d₀.d₋₁d₋₂ (digits separated by spaces when p > 10). */
  toString(count = 20) {
    if (this.isZero()) return "0";
    const { valuation: v, digits } = this.digits(count);
    const sep = this.p > 10n ? " " : "";
    const ds = digits.map(String);
    // position i of ds holds the digit of p^(v + i)
    let whole = [];
    let frac = [];
    for (let i = 0; i < ds.length; i++) (v + i >= 0 ? whole : frac).push(ds[i]);
    if (v > 0) whole = [...Array(v).fill("0"), ...whole];
    const left = whole.reverse().join(sep) || "0";
    const right = frac.reverse().join(sep);
    return `…${left}${right ? "." + right : ""}`;
  }
  /** Equal to absolute precision N: v_p(x − y) ≥ N. */
  equalsTo(y, N) {
    const d = this.sub(y).approx(N);
    return vq(d, this.p) >= N;
  }

  // ---- arithmetic ----
  #same(y) {
    const z = y instanceof Padic ? y : padic(y, this.p);
    if (z.p !== this.p) refuse("bad-argument", `a ${this.p}-adic number and a ${z.p}-adic number cannot be combined`);
    return z;
  }
  add(y) {
    const z = this.#same(y);
    if (this.exact && z.exact) return padic(this.exact.add(z.exact), this.p);
    return new Padic(this.p, (N) => this.approx(N).add(z.approx(N)), Math.min(this.vmin, z.vmin));
  }
  neg() {
    if (this.exact) return padic(this.exact.negate(), this.p);
    return new Padic(this.p, (N) => this.approx(N).negate(), this.vmin);
  }
  sub(y) {
    return this.add(this.#same(y).neg());
  }
  mul(y) {
    const z = this.#same(y);
    if (this.exact && z.exact) return padic(this.exact.multiply(z.exact), this.p);
    const vx = this.exact ? this.valuation() : this.vmin;
    const vy = z.exact ? z.valuation() : z.vmin;
    if (vx === Infinity || vy === Infinity) return padic(0, this.p);
    return new Padic(this.p, (N) => this.approx(Math.max(N - vy, vx)).multiply(z.approx(Math.max(N - vx, vy))), vx + vy);
  }
  /** 1/x: needs the valuation of x, which must be finite. */
  inv() {
    if (this.exact) {
      if (this.exact.numerator === 0n) refuse("bad-argument", "division by 0");
      return padic(this.exact.reciprocal(), this.p);
    }
    const v = this.valuation();
    const p = this.p;
    return new Padic(
      p,
      (N) => {
        // x = p^v·u; 1/x = p^(−v)·u⁻¹, and u⁻¹ is needed to precision N + v
        const k = N + v;
        if (k <= 0) return new Rational(0n);
        const r = this.approx(k + v);
        const u = r.multiply(v >= 0 ? new Rational(1n, p ** BigInt(v)) : new Rational(p ** BigInt(-v)));
        const ui = inverse(residue(u, p, k), p ** BigInt(k));
        return v >= 0 ? new Rational(ui, p ** BigInt(v)) : new Rational(ui * p ** BigInt(-v));
      },
      -v,
    );
  }
  div(y) {
    return this.mul(this.#same(y).inv());
  }
  /** x^k for a whole k (negative k divides). */
  pow(k) {
    let e = big(k, "k");
    if (e < 0n) return this.inv().pow(-e);
    let r = padic(1, this.p);
    let b = this;
    while (e > 0n) {
      if (e & 1n) r = r.mul(b);
      b = b.mul(b);
      e >>= 1n;
    }
    return r;
  }
  /** Is x a square in Q_p? (p odd: even valuation and the unit a residue mod p; p = 2: even valuation and the unit ≡ 1 mod 8.) */
  isSquare() {
    if (this.isZero()) return true;
    const v = this.valuation();
    if (v % 2) return false;
    const u = this.mul(padic(new Rational(1n), this.p).shift(-v));
    if (this.p === 2n) return u.mod(3) === 1n;
    return sqrtMod(u.mod(1), this.p) !== null;
  }
  /** x·p^k. */
  shift(k) {
    const K = Number(k);
    const f = K >= 0 ? new Rational(this.p ** BigInt(K)) : new Rational(1n, this.p ** BigInt(-K));
    return this.mul(padic(f, this.p));
  }
  /**
   * A square root, by Hensel's lemma (Newton's method): the one whose leading
   * digit is the smaller square root mod p (for p = 2, the one ≡ 1 mod 4).
   * Refuses a non-square, with the reason.
   */
  sqrt() {
    if (this.isZero()) return this;
    const p = this.p;
    const v = this.valuation();
    if (v % 2) refuse("not-a-square", `the valuation ${v} is odd`);
    const u = this.shift(-v);
    let r0;
    if (p === 2n) {
      if (u.mod(3) !== 1n) refuse("not-a-square", "a 2-adic unit is a square only when it is 1 mod 8");
      r0 = 1n;
    } else {
      r0 = sqrtMod(u.mod(1), p);
      if (r0 === null) refuse("not-a-square", `${u.mod(1)} is not a square mod ${p}`);
    }
    return henselRoot([u.neg(), padic(0, p), padic(1, p)], r0, p).shift(v / 2);
  }
  /** exp(x) = Σ xⁿ/n!, which converges when v(x) ≥ 1 (v(x) ≥ 2 for p = 2). */
  exp() {
    const p = this.p;
    const need = p === 2n ? 2 : 1;
    if (this.isZero()) return padic(1, p);
    const vx = this.vmin >= need ? this.vmin : this.valuation();
    if (vx < need) refuse("diverges", `exp converges only for v(x) ≥ ${need}; here v(x) = ${vx}`);
    const pm1 = Number(p - 1n);
    return new Padic(
      p,
      (N) => {
        // terms with n·v(x) − (n − 1)/(p − 1) ≥ N are below the precision
        let T = 1;
        while (T * vx - (T - 1) / pm1 < N) T++;
        const a = this.approx(N + 2);
        let term = new Rational(1n);
        let sum = new Rational(1n);
        for (let n = 1; n < T; n++) {
          term = term.multiply(a).divide(new Rational(BigInt(n)));
          sum = sum.add(term);
        }
        return sum;
      },
      0,
    );
  }
  /**
   * The p-adic logarithm. For x ≡ 1 (mod p) (mod 4 when p = 2) it is the series
   * log(1 + y) = y − y²/2 + y³/3 − …; it extends to every unit by
   * log u = log(u^(p−1))/(p−1) (log(u²)/2 for p = 2), and to every x ≠ 0 by
   * Iwasawa's convention log p = 0.
   */
  log() {
    const p = this.p;
    if (this.isZero()) refuse("bad-argument", "log 0 is undefined");
    const v = this.valuation();
    const u = v ? this.shift(-v) : this;
    const need = p === 2n ? 2 : 1;
    const y = u.sub(padic(1, p));
    if (y.isZero()) return padic(0, p);
    let vy;
    try {
      vy = y.valuation(200);
    } catch {
      return padic(0, p);
    }
    if (vy < need) {
      const k = p === 2n ? 2n : p - 1n;
      return u.pow(k).log().div(padic(k, p));
    }
    return new Padic(
      p,
      (N) => {
        // terms with n·v(y) − log_p(n) ≥ N are below the precision
        let T = 1;
        const lp = (n) => Math.floor(Math.log(n) / Math.log(Number(p)) + 1e-9);
        while (!(T * vy - lp(T) >= N && T * vy - lp(T) >= N + 1 - vy)) T++;
        T += 2;
        const a = y.approx(N + lp(T) + 2);
        let pw = new Rational(1n);
        let sum = new Rational(0n);
        for (let n = 1; n <= T; n++) {
          pw = pw.multiply(a);
          const t = pw.divide(new Rational(BigInt(n)));
          sum = n % 2 ? sum.add(t) : sum.subtract(t);
        }
        return sum;
      },
      vy,
    );
  }
}

/**
 * The root of a polynomial (p-adic or rational coefficients, constant first,
 * all p-adic integers) that starts with r0, by Hensel's lemma in its general
 * form: if v(f(r0)) > 2·v(f'(r0)), Newton's iteration from r0 converges to a
 * unique root. The root is a lazy p-adic number.
 */
export function henselRoot(coefficients, r0, p) {
  const P = big(p, "p");
  const cs = coefficients.map((c) => padic(c, P));
  if (cs.some((c) => !c.isZero() && c.vmin < 0 && c.valuation() < 0)) refuse("bad-argument", "the coefficients must be p-adic integers");
  const R0 = big(r0, "r0");
  const at = (M, x) => {
    // f(x) and f'(x) mod p^M with integer coefficients
    const pm = P ** BigInt(M);
    const ints = cs.map((c) => c.mod(M));
    let f = 0n;
    let d = 0n;
    for (let i = ints.length - 1; i >= 0; i--) {
      d = mod(d * x + f, pm);
      f = mod(f * x + ints[i], pm);
    }
    return { f, d, pm };
  };
  const M0 = 64;
  const z = at(M0, R0);
  const dv = val(z.d, P);
  if (dv >= M0) refuse("singular", "f'(r0) vanishes to the working precision; Hensel's lemma does not apply");
  if (val(z.f, P) <= 2 * dv) refuse("bad-argument", `Hensel's condition fails: v(f(r0)) = ${val(z.f, P)} is not more than 2·v(f'(r0)) = ${2 * dv}`);
  return new Padic(
    P,
    (N) => {
      const M = Math.max(N, 1) + 2 * dv + 2;
      let x = R0;
      for (let it = 0; it < 200; it++) {
        const { f, d, pm } = at(M, x);
        if (val(f, P) >= Math.max(N, 1) + 2 * dv || f === 0n) break;
        // x ← x − f/f', dividing out p^dv from both
        const pd = P ** BigInt(dv);
        const w = inverse(mod(d / pd, pm / pd), pm / pd);
        x = mod(x - (f / pd) * w, pm / pd);
      }
      return new Rational(x);
    },
    0,
  );
}

/** The Teichmüller representative of a (not divisible by p): the (p − 1)-th root of unity ≡ a (mod p), as lim a^(p^n). */
export function teichmuller(a, p) {
  const P = big(p, "p");
  const A = mod(big(a, "a"), P);
  if (A === 0n) return padic(0, P);
  return new Padic(P, (N) => {
    const k = Math.max(N, 1);
    const pk = P ** BigInt(k);
    let x = A;
    for (let i = 0; i < k; i++) x = powmodBig(x, P, pk);
    return new Rational(x);
  }, 0);
}
const powmodBig = (a, e, m) => {
  let r = 1n;
  a = mod(a, m);
  while (e > 0n) {
    if (e & 1n) r = (r * a) % m;
    a = (a * a) % m;
    e >>= 1n;
  }
  return r;
};

/**
 * The p-adic expansion of a rational, which is eventually periodic:
 * { valuation, preperiod, period } with digits read from p^valuation upward,
 * so x = Σ dᵢ p^(valuation + i) with the period repeating forever.
 */
export function expansion(x, p) {
  const P = big(p, "p");
  const r = rat(x, "x");
  if (r.numerator === 0n) return { valuation: Infinity, preperiod: [], period: [0n] };
  const v = vq(r, P);
  const u = r.multiply(v >= 0 ? new Rational(1n, P ** BigInt(v)) : new Rational(P ** BigInt(v < 0 ? -v : 0)));
  let n = u.numerator;
  const d = u.denominator;
  const dinv = inverse(mod(d, P), P);
  const seen = new Map();
  const digits = [];
  while (!seen.has(n)) {
    seen.set(n, digits.length);
    const dg = mod(n * dinv, P);
    digits.push(dg);
    n = (n - dg * d) / P;
  }
  const k = seen.get(n);
  return { valuation: v, preperiod: digits.slice(0, k), period: digits.slice(k) };
}

/** v_p(x) of a whole number or fraction, and |x|_p. */
export function padicValuation(x, p) {
  return vq(rat(x, "x"), big(p, "p"));
}
export function padicNorm(x, p) {
  const r = rat(x, "x");
  const P = big(p, "p");
  if (r.numerator === 0n) return new Rational(0n);
  const v = vq(r, P);
  return v >= 0 ? new Rational(1n, P ** BigInt(v)) : new Rational(P ** BigInt(-v));
}

/** A p-adic number from a whole number, a fraction ("3/7", "0.25"), a Rational, or another p-adic number. */
export function padic(x, p) {
  const P = big(p, "p");
  if (P < 2n) throw new RangeError("p must be a prime");
  if (x instanceof Padic) {
    if (x.p !== P) refuse("bad-argument", "a p-adic number for another prime");
    return x;
  }
  const r = rat(x, "x");
  return new Padic(P, () => r, r.numerator === 0n ? 0 : vq(r, P), r);
}
