// The ring of integers of a number field: an integral basis by the Round 2
// algorithm (Zassenhaus; Cohen, A Course in Computational Algebraic Number
// Theory, 6.1), the field discriminant, and the decomposition of rational
// primes into prime ideals by splitting the finite algebra O/pO, which works
// for every prime, including the ones dividing the index where
// Dedekind–Kummer does not apply (Buchmann–Lenstra; Cohen 6.2).
//
// Inside, the field is Q(θ) with θ = c·α a root of a monic integer polynomial
// T, an element is an integer vector over the integral basis ω₁ = 1, …, ωₙ
// (with a common denominator when it is not integral), and an ideal is the
// Hermite normal form of a Z-basis in those coordinates.

import { Rational } from "@ratmath/core";
import { NumberField, NumberFieldElement, countRealRoots, factorModP } from "@ratmath/algebraic";
import { NumberTheoryError, refuse, big, gcd, abs, val } from "./int.js";
import { factor } from "./primes.js";
import { hnfLower, spanMod, reduceMod, sameSpace, kernelMod, det, inverseInt, rowTimes, content, modp, unitVec, identity, isZeroVec } from "./linear.js";

const SUP = { 0: "⁰", 1: "¹", 2: "²", 3: "³", 4: "⁴", 5: "⁵", 6: "⁶", 7: "⁷", 8: "⁸", 9: "⁹", "-": "⁻" };
export const sup = (k) => String(k).replace(/[-\d]/g, (ch) => SUP[ch]);

const cache = new WeakMap();

/**
 * The ring of integers O_K of a number field: K is a NumberField from
 * @ratmath/algebraic or the coefficients (constant first) of an irreducible
 * polynomial. The result is cached on the field.
 */
export function ringOfIntegers(K, options = {}) {
  const F = K instanceof Order ? K.field : K instanceof NumberField ? K : new NumberField(K, options);
  if (K instanceof Order) return K;
  if (!cache.has(F)) cache.set(F, new Order(F));
  return cache.get(F);
}
/** The integral basis of K as field elements, ω₁ = 1 first. */
export const integralBasis = (K) => ringOfIntegers(K).basis;
/** The discriminant of the field K (of its ring of integers). */
export const fieldDiscriminant = (K) => ringOfIntegers(K).discriminant;
/** The prime ideals above p in K: [PrimeIdeal], with p·O = ∏ Pᵢ^eᵢ. */
export const primeDecomposition = (K, p) => ringOfIntegers(K).primesAbove(p);

export class Order {
  constructor(field) {
    this.field = field;
    const n = (this.n = field.n);
    this.degree = n;
    // θ = c·α with a monic integer minimal polynomial T
    const m = field.m; // pairs, monic, constant first
    let c = 1n;
    let L = 1n;
    for (const [, d] of m) L = (L * d) / gcd(L, d);
    for (const { prime: p } of L > 1n ? factor(L) : []) {
      let k = 0;
      for (let i = 1; i <= n; i++) k = Math.max(k, Math.ceil(val(m[n - i][1], p) / i));
      c *= p ** BigInt(k);
    }
    this.c = c;
    this.T = m.map(([a, b], j) => (a * c ** BigInt(n - j)) / b);
    this.polyDiscriminant = this.#discT();
    this.steps = [];
    this.#round2();
    const idx = this.dW ** BigInt(n) / this.W.reduce((s, r, i) => s * r[i], 1n);
    /** [O_K : Z[θ]], θ the integral generator used inside (θ = scale·α). */
    this.index = idx;
    this.scale = c;
    /** The field discriminant. */
    this.discriminant = this.polyDiscriminant / (idx * idx);
    const r1 = countRealRoots(this.T);
    /** [r₁, r₂]: the numbers of real embeddings and of pairs of complex ones. */
    this.signature = [r1, (n - r1) / 2];
    this._primes = new Map();
  }

  // ---- arithmetic in Z[θ] and in the order ----

  /** Product of two polynomials in θ, reduced by T. */
  _pmul(a, b) {
    const n = this.n;
    const r = new Array(2 * n - 1).fill(0n);
    for (let i = 0; i < n; i++) {
      if (a[i] === 0n) continue;
      for (let j = 0; j < n; j++) r[i + j] += a[i] * b[j];
    }
    for (let k = 2 * n - 2; k >= n; k--) {
      const t = r[k];
      if (t === 0n) continue;
      for (let j = 0; j < n; j++) r[k - n + j] -= t * this.T[j];
    }
    return r.slice(0, n);
  }
  #discT() {
    // det of the trace form of 1, θ, …: Newton sums of the roots of T
    const n = this.n;
    const T = this.T;
    const s = [BigInt(n)];
    for (let k = 1; k <= 2 * n - 2; k++) {
      // s_k + e1 s_{k−1} + … = 0 with T = xⁿ + T[n−1]xⁿ⁻¹ + …
      let t = 0n;
      for (let j = 1; j <= Math.min(k, n); j++) t += j === k ? BigInt(k) * T[n - j] : T[n - j] * s[k - j];
      s.push(-t);
    }
    return det(Array.from({ length: n }, (_, i) => Array.from({ length: n }, (_, j) => s[i + j])));
  }
  #setBasis(W, dW) {
    let g = dW;
    for (const row of W) g = gcd(g, content(row));
    this.W = W.map((r) => r.map((x) => x / g));
    this.dW = dW / g;
    const { inv, den } = inverseInt(this.W);
    this.Winv = inv;
    this.dInv = den;
    const n = this.n;
    // multiplication table: ωᵢωⱼ in the basis
    this.mt = this.W.map((wi) =>
      this.W.map((wj) => {
        const prod = this._pmul(wi, wj); // over dW²
        const v = rowTimes(prod, this.Winv); // over dW² · dInv / dW
        const d = this.dW * this.dInv;
        return v.map((x) => {
          if (x % d !== 0n) throw new Error("the basis does not span a ring");
          return x / d;
        });
      }),
    );
    void n;
  }
  /** Product of integer coordinate vectors over the integral basis. */
  _mul(a, b) {
    const n = this.n;
    const out = new Array(n).fill(0n);
    for (let i = 0; i < n; i++) {
      if (a[i] === 0n) continue;
      for (let j = 0; j < n; j++) {
        if (b[j] === 0n) continue;
        const t = a[i] * b[j];
        const row = this.mt[i][j];
        for (let k = 0; k < n; k++) out[k] += t * row[k];
      }
    }
    return out;
  }
  _mulMod(a, b, p) {
    return this._mul(a, b).map((x) => modp(x, p));
  }
  _powMod(a, e, p) {
    let r = unitVec(0, this.n);
    let b = a.map((x) => modp(x, p));
    while (e > 0n) {
      if (e & 1n) r = this._mulMod(r, b, p);
      b = this._mulMod(b, b, p);
      e >>= 1n;
    }
    return r;
  }
  /** Trace of an integer coordinate vector. */
  _trace(a) {
    if (!this._tr) this._tr = this.mt.map((rows) => rows.reduce((s, r, j) => s + r[j], 0n));
    return a.reduce((s, x, k) => s + x * this._tr[k], 0n);
  }
  /** Norm of an integer coordinate vector: the determinant of multiplication by it. */
  _norm(a) {
    return det(identity(this.n).map((e) => this._mul(e, a)));
  }

  // ---- Round 2 ----

  /** The radical of p·O in O/pO, as a list of vectors mod p. */
  _radical(p) {
    const n = this.n;
    if (p > BigInt(n)) {
      // the kernel of the trace form
      const M = identity(n).map((ei) => identity(n).map((ej) => modp(this._trace(this._mul(ei, ej)), p)));
      return kernelMod(M, p);
    }
    let q = p;
    while (q < BigInt(n)) q *= p;
    return kernelMod(identity(n).map((e) => this._powMod(e, q, p)), p);
  }
  #round2() {
    const n = this.n;
    this.#setBasis(identity(n), 1n);
    const D = this.polyDiscriminant;
    const fs = abs(D) > 1n ? factor(D) : [];
    for (const { prime: p, exponent } of fs) {
      if (exponent < 2) continue;
      const gains = [];
      for (;;) {
        const rad = this._radical(p);
        // a Z-basis of I_p = rad + pO, in order coordinates
        const G = hnfLower([...rad, ...identity(n).map((r) => r.map((x) => x * p))], n);
        const { inv: Ginv, den: gden } = inverseInt(G);
        // ωᵢ acting on I_p / pI_p
        const maps = identity(n).map((e) => {
          const flat = [];
          for (const g of G) {
            const v = rowTimes(this._mul(e, g), Ginv);
            for (const x of v) flat.push(modp(x / gden, p));
          }
          return flat;
        });
        const U = kernelMod(maps, p);
        if (!U.length) break;
        gains.push(U.length);
        const rows = [...this.W.map((r) => r.map((x) => x * p)), ...U.map((u) => rowTimes(u, this.W))];
        this.#setBasis(hnfLower(rows, n), this.dW * p);
      }
      const total = gains.reduce((s, x) => s + x, 0);
      this.steps.push({ p, exponent, indexExponent: total, rounds: gains.length, note: total ? `the order grows by index ${p}${total > 1 ? sup(total) : ""} in ${gains.length} round${gains.length > 1 ? "s" : ""}` : `already maximal at ${p}` });
    }
  }

  // ---- elements ----

  /** Coordinates of a field element over the integral basis as { v, d }: the element is Σ vᵢωᵢ / d. */
  _coords(x) {
    if (x && x.v && typeof x.d === "bigint") return x;
    const el = this.field.of(typeof x === "number" ? BigInt(x) : x);
    const n = this.n;
    // α^i = θ^i / c^i
    let d = 1n;
    const cs = Array.from({ length: n }, (_, i) => {
      const [a, b] = el.c[i] ?? [0n, 1n];
      return [a, b * this.c ** BigInt(i)];
    });
    for (const [, b] of cs) d = (d * b) / gcd(d, b);
    const theta = cs.map(([a, b]) => a * (d / b));
    const v = rowTimes(theta, this.Winv); // over d · dInv / dW
    return normal(v.map((t) => t * this.dW), d * this.dInv);
  }
  /** The field element with integer coordinates v over the integral basis, divided by d. */
  _element(v, d = 1n) {
    const theta = rowTimes(v, this.W); // over dW·d
    const den = this.dW * d;
    const pairs = theta.map((t, i) => {
      const a = t * this.c ** BigInt(i);
      const g = gcd(a, den) || 1n;
      return [a / g, den / g];
    });
    return this.field.of(pairs.some(([a]) => a !== 0n) ? pairs : 0n);
  }
  /** The integral basis as field elements; ω₁ = 1. */
  get basis() {
    return (this._basis ??= identity(this.n).map((e) => this._element(e)));
  }
  /** Coordinates of x over the integral basis, as Rationals: whole numbers exactly when x is an algebraic integer. */
  coordinates(x) {
    const { v, d } = this._coords(x);
    return v.map((t) => new Rational(t, d));
  }
  /** The element Σ cᵢωᵢ of the field. */
  element(coords) {
    if (coords.length !== this.n) throw new RangeError(`expected ${this.n} coordinates`);
    let d = 1n;
    const rs = coords.map((t) => (t instanceof Rational ? t : new Rational(typeof t === "number" ? BigInt(t) : t)));
    for (const r of rs) d = (d * r.denominator) / gcd(d, r.denominator);
    return this._element(rs.map((r) => r.numerator * (d / r.denominator)), d);
  }
  /** Is x an algebraic integer? */
  contains(x) {
    return this._coords(x).d === 1n;
  }
  /** The matrix of traces Tr(ωᵢωⱼ); its determinant is the discriminant. */
  get traceMatrix() {
    return identity(this.n).map((ei) => identity(this.n).map((ej) => this._trace(this._mul(ei, ej))));
  }

  // ---- primes ----

  /** The prime ideals above the rational prime p, with p·O = ∏ P^e. */
  primesAbove(p) {
    p = big(p, "p");
    if (this._primes.has(p)) return this._primes.get(p);
    if (p < 2n || factor(p).length !== 1 || factor(p)[0].exponent !== 1) refuse("bad-argument", `${p} is not a prime`);
    const n = this.n;
    const one = unitVec(0, n);
    const found = [];
    const split = (H) => {
      // H: subspace of O/pO containing the radical; O/H is a product of fields
      const free = [];
      for (let j = 0; j < n; j++) if (!H.pivots.includes(j)) free.push(j);
      const m = free.length;
      const toFree = (v) => free.map((j) => v[j]);
      const frob = free.map((j) => {
        const e = unitVec(j, n);
        const y = reduceMod(this._powMod(e, p, p).map((x, k) => x - e[k]), H);
        return toFree(y);
      });
      const ker = kernelMod(frob, p).map((k) => {
        const v = new Array(n).fill(0n);
        k.forEach((x, i) => (v[free[i]] = x));
        return v;
      });
      if (ker.length === 1) {
        found.push({ H, f: m });
        return;
      }
      // an element fixed by Frobenius that is not a constant mod H
      const oneR = reduceMod(one, H);
      const alpha = ker.find((k) => spanMod([k, oneR], p).rows.length === 2);
      // its minimal polynomial over F_p in O/H
      const pows = [oneR];
      let rel = null;
      while (!rel) {
        const next = reduceMod(this._mulMod(pows[pows.length - 1], alpha, p), H);
        pows.push(next);
        const k = kernelMod(pows, p);
        if (k.length) rel = k[0];
      }
      const roots = rootsOfSplit(rel, p);
      for (const r of roots) {
        const shifted = alpha.map((x, k) => modp(x - r * one[k], p));
        split(spanMod([...H.rows, ...identity(n).map((e) => this._mulMod(e, shifted, p))], p));
      }
    };
    split(spanMod(this._radical(p), p));
    const primes = found.map(({ H, f }) => new PrimeIdeal(this, p, H, f));
    primes.sort((a, b) => a.f - b.f || a.e - b.e || cmpVec(a.pi, b.pi));
    primes.forEach((P, i) => (P.index = i));
    const total = primes.reduce((s, P) => s + P.e * P.f, 0);
    if (total !== n) throw new Error(`Σ e·f = ${total} ≠ ${n} above ${p}`);
    this._primes.set(p, primes);
    return primes;
  }
  /** "p·O = P₁² · P₂", the decomposition of p in words, with e and f for each prime. */
  describePrime(p) {
    const Ps = this.primesAbove(p);
    const kind = Ps.length === 1 && Ps[0].e === 1 ? "inert" : Ps.every((P) => P.e === 1) ? (Ps.length === this.n ? "splits completely" : "splits") : Ps.length === 1 && Ps[0].f === 1 && this.n > 2 ? "totally ramified" : "ramified";
    return { p: big(p), kind, primes: Ps, text: `${p}·O = ${Ps.map((P) => `${P}${P.e > 1 ? sup(P.e) : ""}`).join(" · ")}` };
  }
  /** The primes that ramify: the prime divisors of the discriminant. */
  get ramifiedPrimes() {
    return abs(this.discriminant) > 1n ? factor(this.discriminant).map((x) => x.prime) : [];
  }

  // ---- ideals ----

  /** The (fractional) ideal generated by the given elements; a PrimeIdeal or Ideal is passed through. */
  ideal(...gens) {
    if (gens.length === 1 && gens[0] instanceof Ideal) return gens[0];
    if (gens.length === 1 && gens[0] instanceof PrimeIdeal) return gens[0].toIdeal();
    if (gens.length === 1 && Array.isArray(gens[0]) && gens[0].some((g) => g instanceof NumberFieldElement)) gens = gens[0];
    const cs = gens.map((g) => this._coords(g)).filter((c) => !isZeroVec(c.v));
    if (!cs.length) refuse("bad-argument", "the zero ideal is not supported");
    let d = 1n;
    for (const c of cs) d = (d * c.d) / gcd(d, c.d);
    const rows = [];
    for (const c of cs) {
      const v = c.v.map((x) => x * (d / c.d));
      for (const e of identity(this.n)) rows.push(this._mul(e, v));
    }
    return Ideal._make(this, hnfLower(rows, this.n), d);
  }
  /** The unit ideal O. */
  get unitIdeal() {
    return (this._unit ??= Ideal._make(this, identity(this.n), 1n));
  }

  toString() {
    return `the ring of integers of ${this.field}: discriminant ${this.discriminant}, basis ${this.basis.join(", ")}`;
  }
}

const cmpVec = (a, b) => {
  for (let i = 0; i < a.length; i++) if (a[i] !== b[i]) return a[i] < b[i] ? -1 : 1;
  return 0;
};
function normal(v, d) {
  if (d < 0n) {
    d = -d;
    v = v.map((x) => -x);
  }
  let g = d;
  for (const x of v) g = gcd(g, x);
  return g > 1n ? { v: v.map((x) => x / g), d: d / g } : { v, d };
}
/** Roots in F_p of a polynomial (coefficients constant first) that splits into distinct linear factors. */
function rootsOfSplit(f, p) {
  if (p < 200n) {
    const out = [];
    for (let r = 0n; r < p; r++) {
      let s = 0n;
      for (let i = f.length - 1; i >= 0; i--) s = (s * r + f[i]) % p;
      if (s === 0n) out.push(r);
    }
    return out;
  }
  return factorModP(f, p).map((g) => modp(-g[0], p));
}

export class PrimeIdeal {
  constructor(order, p, H, f) {
    this.order = order;
    /** The rational prime below. */
    this.p = p;
    /** The residue degree: the norm is p^f. */
    this.f = f;
    this.space = H; // P / pO as a subspace of O / pO
    const n = order.n;
    // a second generator: π with span(π·ωᵢ) = P mod p
    this.pi = this.#generator();
    // τ with τ·P ⊂ pO, τ ∉ pO: τ/p generates P⁻¹ over O
    if (!H.rows.length) this.tau = unitVec(0, n);
    else {
      const maps = identity(n).map((e) => H.rows.flatMap((b) => order._mulMod(e, b, p)));
      this.tau = kernelMod(maps, p)[0];
    }
    /** The ramification index. */
    this.e = this._valInt(unitVec(0, n).map((x) => x * p));
    /** The norm p^f. */
    this.norm = p ** BigInt(f);
  }
  #generator() {
    const { order: O, p, space: H } = this;
    const n = O.n;
    if (!H.rows.length) return new Array(n).fill(0n);
    const good = (v) => sameSpace(spanMod(identity(n).map((e) => O._mulMod(e, v, p)), p), H);
    // Dedekind–Kummer first: g(θ) for the minimal polynomial g of θ mod P, when θ is in the basis span
    if (O.index % p !== 0n && n > 1) {
      const theta = thetaCoords(O);
      const pows = [reduceMod(unitVec(0, n), H)];
      let rel = null;
      while (!rel) {
        pows.push(reduceMod(O._mulMod(pows[pows.length - 1], theta.v, p), H));
        const k = kernelMod(pows, p);
        if (k.length) rel = k[0];
      }
      // normalise to monic and evaluate at θ with centred coefficients
      if (rel[rel.length - 1] !== 0n) {
        const deg = rel.length - 1;
        const inv = modInv(rel[deg], p);
        const g = rel.map((x) => (x * inv) % p);
        let acc = new Array(n).fill(0n);
        for (let i = deg; i >= 0; i--) {
          acc = O._mul(acc, theta.v);
          acc[0] += g[i];
        }
        if (good(acc)) {
          this.kummer = g;
          return acc;
        }
      }
    }
    // otherwise search the subspace, smallest combinations first
    const B = H.rows;
    for (const b of B) if (good(b)) return centre(b, p);
    const lim = p < 5n ? Number(p) : 5;
    for (let bound = 2; bound <= lim; bound++) {
      const idx = new Array(B.length).fill(0);
      for (;;) {
        let k = 0;
        while (k < B.length && ++idx[k] === bound) idx[k++] = 0;
        if (k === B.length) break;
        if (!idx.includes(bound - 1)) continue;
        const v = new Array(n).fill(0n);
        idx.forEach((c, i) => {
          if (c) for (let j = 0; j < n; j++) v[j] += BigInt(c) * B[i][j];
        });
        if (good(v)) return centre(v, p);
      }
    }
    // random combinations for large p
    let seed = 12345n;
    for (let t = 0; t < 10000; t++) {
      const v = new Array(n).fill(0n);
      for (const b of B) {
        seed = (seed * 6364136223846793005n + 1442695040888963407n) % (1n << 64n);
        const c = (seed >> 20n) % p;
        for (let j = 0; j < n; j++) v[j] += c * b[j];
      }
      if (good(v)) return centre(v, p);
    }
    throw new Error("no second generator found");
  }
  /** v_P of an integer coordinate vector (not zero). */
  _valInt(v) {
    const { order: O, p } = this;
    if (isZeroVec(v)) return Infinity;
    const c = content(v);
    let k = 0;
    const vc = val(c, p);
    // the content contributes e·v_p(content); before e is known (the constructor) run the loop on everything
    let x = v;
    if (this.e !== undefined && vc > 0) {
      const q = p ** BigInt(vc);
      x = v.map((t) => t / q);
      k = vc * this.e;
    }
    for (;;) {
      const y = O._mul(x, this.tau);
      if (!y.every((t) => t % p === 0n)) return k;
      x = y.map((t) => t / p);
      k++;
    }
  }
  /** The valuation at this prime of a field element, an ideal or a prime ideal. */
  valuation(x) {
    if (x instanceof PrimeIdeal) return x === this || (x.order === this.order && x.p === this.p && sameSpace(x.space, this.space)) ? 1 : 0;
    if (x instanceof Ideal) {
      let best = Infinity;
      for (const row of x.H) best = Math.min(best, this._valInt(row));
      return best - this.e * val(x.d, this.p);
    }
    const { v, d } = this.order._coords(x);
    if (isZeroVec(v)) return Infinity;
    return this._valInt(v) - this.e * val(d, this.p);
  }
  /** The two generators: p and a second element (just p when the prime is inert). */
  get generators() {
    const O = this.order;
    return this.space.rows.length ? [O._element(unitVec(0, O.n).map((x) => x * this.p)), O._element(this.pi)] : [O._element(unitVec(0, O.n).map((x) => x * this.p))];
  }
  /** This prime as an Ideal, for products and sums. */
  toIdeal() {
    const O = this.order;
    return (this._ideal ??= Ideal._make(O, hnfLower([...this.space.rows, ...identity(O.n).map((r) => r.map((x) => x * this.p))], O.n), 1n));
  }
  /** P⁻¹ = O + (τ/p)·O. */
  _inverse() {
    const O = this.order;
    return (this._inv ??= Ideal._make(O, hnfLower([...identity(O.n).map((r) => r.map((x) => x * this.p)), ...identity(O.n).map((e) => O._mul(e, this.tau))], O.n), this.p));
  }
  equals(Q) {
    return Q instanceof PrimeIdeal && this.valuation(Q) === 1;
  }
  toString() {
    const gs = this.generators;
    return `(${gs.join(", ")})`;
  }
}
function thetaCoords(O) {
  // θ = c·α in order coordinates
  return O._coords(O.field.gen.mul(O.field.of(O.c)));
}
function modInv(a, p) {
  let [r0, r1, x0, x1] = [modp(a, p), p, 1n, 0n];
  while (r1) {
    const q = r0 / r1;
    [r0, r1, x0, x1] = [r1, r0 - q * r1, x1, x0 - q * x1];
  }
  return modp(x0, p);
}
const centre = (v, p) => v.map((x) => (2n * x > p ? x - p : x));

export class Ideal {
  static _make(order, H, d) {
    const I = new Ideal();
    let g = d;
    for (const row of H) g = gcd(g, content(row));
    I.order = order;
    I.H = g > 1n ? H.map((r) => r.map((x) => x / g)) : H;
    I.d = d / g;
    return I;
  }
  /** The norm: [O : I] for an integral ideal, extended multiplicatively. */
  get norm() {
    return new Rational(this.H.reduce((s, r, i) => s * r[i], 1n), this.d ** BigInt(this.order.n));
  }
  get isIntegral() {
    return this.d === 1n;
  }
  get isOne() {
    return this.d === 1n && this.H.every((r, i) => r[i] === 1n);
  }
  /** A Z-basis of the ideal, as field elements. */
  get basis() {
    return this.H.map((row) => this.order._element(row, this.d));
  }
  /** The smallest positive rational number in the ideal. */
  get minimum() {
    return new Rational(this.H[0][0], this.d);
  }
  mul(J) {
    const O = this.order;
    J = O.ideal(J);
    if (this.isOne) return J;
    if (J.isOne) return this;
    const rows = [];
    for (const a of this.H) for (const b of J.H) rows.push(O._mul(a, b));
    return Ideal._make(O, hnfLower(rows, O.n), this.d * J.d);
  }
  add(J) {
    const O = this.order;
    J = O.ideal(J);
    const rows = [...this.H.map((r) => r.map((x) => x * J.d)), ...J.H.map((r) => r.map((x) => x * this.d))];
    return Ideal._make(O, hnfLower(rows, O.n), this.d * J.d);
  }
  /** The inverse ideal, through the prime factorisation. */
  inverse() {
    let out = this.order.unitIdeal;
    for (const { prime, exponent } of this.factor()) out = out.mul(powIdeal(exponent > 0 ? prime._inverse() : prime.toIdeal(), Math.abs(exponent)));
    return out;
  }
  div(J) {
    return this.mul(this.order.ideal(J).inverse());
  }
  pow(k) {
    k = Number(k);
    return k >= 0 ? powIdeal(this, k) : powIdeal(this.inverse(), -k);
  }
  /** Does the ideal contain the element x? */
  contains(x) {
    const { v, d } = this.order._coords(x);
    // x = c·H/dI  ⇔  v·dI/d ∈ row lattice of H
    const { inv, den } = (this._inv ??= inverseInt(this.H));
    return rowTimes(v, inv).every((t) => (t * this.d) % (den * d) === 0n);
  }
  equals(J) {
    J = this.order.ideal(J);
    return this.d === J.d && this.H.every((r, i) => r.every((x, j) => x === J.H[i][j]));
  }
  /** The factorisation into prime ideals: [{ prime, exponent }], negative exponents for a fractional ideal. */
  factor() {
    if (this._factors) return this._factors;
    const O = this.order;
    const N = this.H.reduce((s, r, i) => s * r[i], 1n) * this.d;
    const out = [];
    for (const { prime: p } of N > 1n ? factor(N) : []) {
      for (const P of O.primesAbove(p)) {
        const v = P.valuation(this);
        if (v !== 0) out.push({ prime: P, exponent: v });
      }
    }
    return (this._factors = out);
  }
  /** Is this a prime ideal? */
  get isPrime() {
    const fs = this.factor();
    return fs.length === 1 && fs[0].exponent === 1;
  }
  toString() {
    if (this.isOne) return "(1)";
    const fs = this.factor();
    return fs.map(({ prime, exponent }) => `${prime}${exponent === 1 ? "" : sup(exponent)}`).join(" · ");
  }
}
function powIdeal(I, k) {
  let out = I.order.unitIdeal;
  let b = I;
  while (k > 0) {
    if (k & 1) out = out.mul(b);
    k >>= 1;
    if (k) b = b.mul(b);
  }
  return out;
}

export { NumberTheoryError };
