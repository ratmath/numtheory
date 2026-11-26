// Elliptic curves over Q and over the prime fields F_p.
//
// A curve is given by a Weierstrass equation
//     y² + a1·xy + a3·y = x³ + a2·x² + a4·x + a6,
// or by the short form y² = x³ + ax + b. Points are frozen objects { x, y },
// with Rational coordinates over Q and BigInt residues over F_p; the point at
// infinity, the zero of the group, is curve.infinity = { infinity: true }.
//
// Over Q: the group law, torsion (Nagell–Lutz, with Mazur's theorem as a
// check), minimal models, reduction types, a search for points of small
// height, and enclosures of the canonical height. Ranks are out of scope.
// Over F_p: the group law, point counts (by character sum, by Shanks's baby
// steps and giant steps with Mestre's twist, and by Schoof's algorithm),
// orders of points and, for small p, the group structure.

import { Rational } from "@ratmath/core";
import { isProbablePrime, sqrtMod } from "@ratmath/algebraic";
import { abs, big, gcd, isqrt, isSquare, lcm, mod, rat, refuse, val } from "./int.js";
import { factor, PROVEN_BELOW } from "./primes.js";
import { ecAdd, ecMul } from "./ecmod.js";
import { invMod } from "./polymod.js";
import { divisionPolynomials, schoofTrace } from "./schoof.js";

const INF = Object.freeze({ infinity: true });
const R = (n) => new Rational(BigInt(n));

const QQ = {
  zero: R(0),
  one: R(1),
  int: (n) => R(n),
  from: (x, what) => rat(x, what),
  add: (a, b) => a.add(b),
  sub: (a, b) => a.subtract(b),
  mul: (a, b) => a.multiply(b),
  neg: (a) => a.negate(),
  inv: (a) => a.reciprocal(),
  eq: (a, b) => a.equals(b),
  isZero: (a) => a.numerator === 0n,
};
const FF = (p) => ({
  zero: 0n,
  one: 1n % p,
  int: (n) => mod(BigInt(n), p),
  from: (x, what) => {
    if (typeof x === "bigint") return mod(x, p);
    if (typeof x === "number" || (typeof x === "string" && /^\s*[-+]?\d+\s*$/.test(x))) return mod(big(x, what), p);
    const q = rat(x, what);
    if (q.denominator % p === 0n) refuse("bad-argument", `${what} = ${q} has a denominator divisible by ${p}`);
    return mod(q.numerator * invMod(q.denominator, p), p);
  },
  add: (a, b) => (a + b) % p,
  sub: (a, b) => mod(a - b, p),
  mul: (a, b) => (a * b) % p,
  neg: (a) => (a === 0n ? 0n : p - a),
  inv: (a) => invMod(a, p),
  eq: (a, b) => a === b,
  isZero: (a) => a === 0n,
});

const onlyQ = (E, what) => {
  if (E.modulus !== null) refuse("bad-argument", `${what} is for curves over Q; this curve is over F_${E.modulus}`);
};
const onlyFp = (E, what) => {
  if (E.modulus === null) refuse("bad-argument", `${what} is for curves over F_p; reduce with curve.modulo(p) first`);
};

export class EllipticCurve {
  /** Use ellipticCurve(coefficients, options). */
  constructor(a, K, modulus) {
    this.K = K;
    this.modulus = modulus;
    [this.a1, this.a2, this.a3, this.a4, this.a6] = a;
    const { add, sub, mul, int } = K;
    const { a1, a2, a3, a4, a6 } = this;
    const sq = (x) => mul(x, x);
    this.b2 = add(sq(a1), mul(int(4), a2));
    this.b4 = add(mul(int(2), a4), mul(a1, a3));
    this.b6 = add(sq(a3), mul(int(4), a6));
    this.b8 = sub(add(add(mul(sq(a1), a6), mul(int(4), mul(a2, a6))), mul(a2, sq(a3))), add(mul(a1, mul(a3, a4)), sq(a4)));
    const { b2, b4, b6, b8 } = this;
    this.c4 = sub(sq(b2), mul(int(24), b4));
    this.c6 = sub(add(K.neg(mul(b2, sq(b2))), mul(int(36), mul(b2, b4))), mul(int(216), b6));
    // Δ = −b2²b8 − 8b4³ − 27b6² + 9b2b4b6
    this.discriminant = add(
      sub(sub(K.neg(mul(sq(b2), b8)), mul(int(8), mul(b4, sq(b4)))), mul(int(27), sq(b6))),
      mul(int(9), mul(b2, mul(b4, b6))),
    );
    if (K.isZero(this.discriminant)) refuse("singular", "the discriminant is 0: the cubic has a node or a cusp, so this is not an elliptic curve");
    this.j = mul(mul(this.c4, sq(this.c4)), K.inv(this.discriminant));
    this.infinity = INF;
    this._cache = {};
    Object.freeze(this);
  }
  get field() {
    return this.modulus === null ? "Q" : `F_${this.modulus}`;
  }
  get coefficients() {
    return [this.a1, this.a2, this.a3, this.a4, this.a6];
  }
  get isShort() {
    const z = this.K.isZero;
    return z(this.a1) && z(this.a2) && z(this.a3);
  }
  toString() {
    const t = (c, m) => {
      if (this.K.isZero(c)) return "";
      const s = c.toString();
      const neg = s.startsWith("-");
      const body = neg ? s.slice(1) : s;
      return ` ${neg ? "−" : "+"} ${body === "1" && m ? "" : body}${m}`;
    };
    const lhs = `y²${t(this.a1, "xy")}${t(this.a3, "y")}`;
    const rhs = `x³${t(this.a2, "x²")}${t(this.a4, "x")}${t(this.a6, "")}`;
    return `${lhs} = ${rhs}${this.modulus === null ? "" : ` over F_${this.modulus}`}`;
  }

  // ---- points and the group law ----

  _pt(P) {
    if (P === null || P === undefined || P.infinity) return INF;
    if (Array.isArray(P)) return this.point(P[0], P[1]);
    if (Object.isFrozen(P) && this._isElt(P.x) && this._isElt(P.y)) return P;
    return this.point(P.x, P.y);
  }
  _isElt(v) {
    return this.modulus === null ? v instanceof Rational : typeof v === "bigint" && v >= 0n && v < this.modulus;
  }
  /** The point (x, y), checked to lie on the curve. */
  point(x, y) {
    const P = Object.freeze({ x: this.K.from(x, "x"), y: this.K.from(y, "y") });
    if (!this._on(P)) refuse("not-on-curve", `(${P.x}, ${P.y}) is not on ${this}`);
    return P;
  }
  _on(P) {
    const { add, mul } = this.K;
    const { x, y } = P;
    const lhs = add(mul(y, y), add(mul(this.a1, mul(x, y)), mul(this.a3, y)));
    const rhs = add(mul(add(mul(add(x, this.a2), x), this.a4), x), this.a6);
    return this.K.eq(lhs, rhs);
  }
  /** True when P (a point, or [x, y]) lies on the curve. */
  has(P) {
    if (P === null || P === undefined || P.infinity) return true;
    const [x, y] = Array.isArray(P) ? P : [P.x, P.y];
    return this._on({ x: this.K.from(x, "x"), y: this.K.from(y, "y") });
  }
  equal(P, Q) {
    P = this._pt(P);
    Q = this._pt(Q);
    if (P.infinity || Q.infinity) return Boolean(P.infinity && Q.infinity);
    return this.K.eq(P.x, Q.x) && this.K.eq(P.y, Q.y);
  }
  /** −P = (x, −y − a1·x − a3). */
  neg(P) {
    P = this._pt(P);
    if (P.infinity) return INF;
    const K = this.K;
    return Object.freeze({ x: P.x, y: K.sub(K.neg(P.y), K.add(K.mul(this.a1, P.x), this.a3)) });
  }
  /** P + Q by the chord and tangent rule. */
  add(P, Q) {
    P = this._pt(P);
    Q = this._pt(Q);
    if (P.infinity) return Q;
    if (Q.infinity) return P;
    const K = this.K;
    const { add, sub, mul, int } = K;
    const { a1, a2, a3, a4, a6 } = this;
    const [x1, y1, x2, y2] = [P.x, P.y, Q.x, Q.y];
    let lam;
    let nu;
    if (K.eq(x1, x2)) {
      if (K.isZero(add(add(y1, y2), add(mul(a1, x1), a3)))) return INF;
      const di = K.inv(add(add(mul(int(2), y1), mul(a1, x1)), a3));
      lam = mul(sub(add(add(mul(int(3), mul(x1, x1)), mul(int(2), mul(a2, x1))), a4), mul(a1, y1)), di);
      nu = mul(sub(add(add(K.neg(mul(x1, mul(x1, x1))), mul(a4, x1)), mul(int(2), a6)), mul(a3, y1)), di);
    } else {
      const di = K.inv(sub(x2, x1));
      lam = mul(sub(y2, y1), di);
      nu = mul(sub(mul(y1, x2), mul(y2, x1)), di);
    }
    const x3 = sub(sub(sub(add(mul(lam, lam), mul(a1, lam)), a2), x1), x2);
    const y3 = sub(sub(K.neg(mul(add(lam, a1), x3)), nu), a3);
    return Object.freeze({ x: x3, y: y3 });
  }
  sub(P, Q) {
    return this.add(P, this.neg(Q));
  }
  double(P) {
    return this.add(P, P);
  }
  /** k·P for any whole k, by doubling and adding. */
  mul(k, P) {
    let n = big(k, "k");
    P = this._pt(P);
    if (n < 0n) {
      n = -n;
      P = this.neg(P);
    }
    let Rr = INF;
    if (n === 0n || P.infinity) return Rr;
    const bits = n.toString(2);
    for (let i = 0; i < bits.length; i++) {
      Rr = this.add(Rr, Rr);
      if (bits[i] === "1") Rr = this.add(Rr, P);
    }
    return Rr;
  }
  /**
   * The order of P. Over F_p a BigInt (found from the group order). Over Q a
   * BigInt from 1 to 12 for a torsion point and Infinity otherwise (Mazur:
   * no rational point has finite order above 12).
   */
  order(P) {
    P = this._pt(P);
    if (P.infinity) return 1n;
    if (this.modulus === null) {
      let Q = P;
      for (let k = 1n; k <= 12n; k++) {
        if (Q.infinity) return k;
        Q = this.add(Q, P);
      }
      return Infinity;
    }
    let n = this.count();
    for (const { prime } of factor(n)) while (n % prime === 0n && this.mul(n / prime, P).infinity) n /= prime;
    return n;
  }

  // ---- models over Q ----

  /** True when every coefficient is a whole number. */
  get isIntegral() {
    onlyQ(this, "isIntegral");
    return this.coefficients.every((c) => c.denominator === 1n);
  }
  /**
   * A short model Y² = X³ + AX + B with whole A and B, and the maps both
   * ways: { curve, A, B, to(P), from(P) }. X = s²(36x + 3b2) and
   * Y = 108·s³(2y + a1·x + a3) for a rational s chosen to clear denominators.
   */
  shortModel() {
    onlyQ(this, "shortModel");
    if (this._cache.short) return this._cache.short;
    let A = this.c4.multiply(R(-27));
    let B = this.c6.multiply(R(-54));
    const d = lcm(A.denominator, B.denominator);
    let [An, Bn] = [A.numerator * (d ** 4n / A.denominator), B.numerator * (d ** 6n / B.denominator)];
    let s = new Rational(d, 1n);
    for (const p of [2n, 3n]) {
      while (An % p ** 4n === 0n && Bn % p ** 6n === 0n) {
        An /= p ** 4n;
        Bn /= p ** 6n;
        s = s.divide(R(p));
      }
    }
    const curve = new EllipticCurve([R(0), R(0), R(0), R(An), R(Bn)], QQ, null);
    const s2 = s.multiply(s);
    const s3 = s2.multiply(s);
    const to = (P) => {
      P = this._pt(P);
      if (P.infinity) return INF;
      const X = s2.multiply(P.x.multiply(R(36)).add(this.b2.multiply(R(3))));
      const Y = s3.multiply(R(108)).multiply(P.y.multiply(R(2)).add(this.a1.multiply(P.x)).add(this.a3));
      return curve.point(X, Y);
    };
    const from = (P) => {
      P = curve._pt(P);
      if (P.infinity) return INF;
      const x = P.x.divide(s2).subtract(this.b2.multiply(R(3))).divide(R(36));
      const y = P.y.divide(s3.multiply(R(108))).subtract(this.a1.multiply(x)).subtract(this.a3).divide(R(2));
      return this.point(x, y);
    };
    return (this._cache.short = Object.freeze({ curve, A: An, B: Bn, to, from }));
  }

  /**
   * A global minimal model: whole coefficients, a1 and a3 in {0, 1}, a2 in
   * {−1, 0, 1}, and a discriminant that no change of variables over Q makes
   * smaller. Returns { curve, u, r, s, t, to(P), from(P) }, where the old
   * coordinates are x = u²x′ + r, y = u³y′ + u²s·x′ + t.
   *
   * At each prime the largest u is found by trying it: the candidate (c4, c6)
   * either are the invariants of a model with whole coefficients, which is
   * then written down and checked, or they are not.
   */
  minimalModel() {
    onlyQ(this, "minimalModel");
    if (this._cache.minimal) return this._cache.minimal;
    // an integral model first: a_i → w^i·a_i
    const w = this.coefficients.reduce((m, c) => lcm(m, c.denominator), 1n);
    const c4 = (this.c4.numerator * w ** 4n) / this.c4.denominator;
    const c6 = (this.c6.numerator * w ** 6n) / this.c6.denominator;
    const g = gcd(c4, c6);
    const cap = (p) => {
      const v4 = c4 === 0n ? Infinity : Math.floor(val(c4, p) / 4);
      const v6 = c6 === 0n ? Infinity : Math.floor(val(c6, p) / 6);
      return Math.min(v4, v6);
    };
    let u = 1n;
    if (g > 1n) for (const { prime } of factor(g)) if (prime > 3n) u *= prime ** BigInt(cap(prime));
    const attempt = (uu) => reducedFrom(c4 / uu ** 4n, c6 / uu ** 6n);
    for (const p of [2n, 3n]) {
      let d = cap(p);
      while (d > 0 && !attempt(u * p ** BigInt(d))) d--;
      u *= p ** BigInt(d);
    }
    const model = attempt(u);
    if (!model) throw new Error("minimalModel: no integral model found; this should not happen");
    const curve = new EllipticCurve(model.map(R), QQ, null);
    // the change of variables from this curve to the minimal one
    const uq = new Rational(u, w);
    const s = uq.multiply(curve.a1).subtract(this.a1).divide(R(2));
    const r = uq.pow(2).multiply(curve.a2).subtract(this.a2).add(s.multiply(this.a1)).add(s.multiply(s)).divide(R(3));
    const t = uq.pow(3).multiply(curve.a3).subtract(this.a3).subtract(r.multiply(this.a1)).divide(R(2));
    const [u2, u3] = [uq.pow(2), uq.pow(3)];
    const to = (P) => {
      P = this._pt(P);
      if (P.infinity) return INF;
      const x = P.x.subtract(r).divide(u2);
      return curve.point(x, P.y.subtract(t).subtract(u2.multiply(s).multiply(x)).divide(u3));
    };
    const from = (P) => {
      P = curve._pt(P);
      if (P.infinity) return INF;
      return this.point(u2.multiply(P.x).add(r), u3.multiply(P.y).add(u2.multiply(s).multiply(P.x)).add(t));
    };
    return (this._cache.minimal = Object.freeze({ curve, u: uq, r, s, t, to, from }));
  }
  /** The primes of bad reduction: those dividing the minimal discriminant. */
  badPrimes() {
    onlyQ(this, "badPrimes");
    return factor(this.minimalModel().curve.discriminant.numerator).map((f) => f.prime);
  }
  /**
   * How the curve reduces at the prime p:
   * { prime, type: "good" | "multiplicative" | "additive", split, discriminantValuation, conductorExponent }.
   * `split` is set for multiplicative reduction (a node): true when the two
   * tangents at the node are defined over F_p. The conductor exponent is 0,
   * 1, or 2 for additive reduction at p ≥ 5; for additive reduction at 2 or 3
   * it needs Tate's algorithm, which is not here, and is null.
   */
  reduction(p) {
    onlyQ(this, "reduction");
    p = big(p, "p");
    if (p < 2n || !isProbablePrime(p)) refuse("bad-argument", `${p} is not a prime`);
    const E = this.minimalModel().curve;
    const v = val(E.discriminant.numerator, p);
    if (v === 0) return { prime: p, type: "good", discriminantValuation: 0, conductorExponent: 0 };
    const c4 = E.c4.numerator;
    if (mod(c4, p) !== 0n) {
      const split = p <= 3n ? nodeIsSplit(E, p) : legendreSymbol(-E.c6.numerator, p) === 1;
      return { prime: p, type: "multiplicative", split, discriminantValuation: v, conductorExponent: 1 };
    }
    return { prime: p, type: "additive", discriminantValuation: v, conductorExponent: p >= 5n ? 2 : null };
  }
  /** The conductor, when no prime has additive reduction at 2 or 3 (those exponents need Tate's algorithm). */
  conductor() {
    onlyQ(this, "conductor");
    let N = 1n;
    for (const p of this.badPrimes()) {
      const e = this.reduction(p).conductorExponent;
      if (e === null) refuse("too-hard", `the reduction at ${p} is additive; its conductor exponent needs Tate's algorithm, which is not implemented`);
      N *= p ** BigInt(e);
    }
    return N;
  }
  /** This equation read modulo the prime p, as a curve over F_p (p must not divide a denominator or the discriminant). */
  modulo(p) {
    onlyQ(this, "modulo");
    p = big(p, "p");
    return ellipticCurve(this.coefficients, { mod: p });
  }

  // ---- torsion over Q ----

  /**
   * The torsion subgroup: { order, structure, points, generators }.
   * structure is [n] for Z/n or [2, m] for Z/2 × Z/m. Found by Nagell–Lutz on
   * a short integral model: a torsion point there has whole coordinates and
   * y = 0 or y² dividing 4A³ + 27B². Each candidate is kept only if one of
   * its first 12 multiples is the zero (Mazur's bound).
   */
  torsion() {
    onlyQ(this, "torsion");
    if (this._cache.torsion) return this._cache.torsion;
    const { curve: S, A, B, from } = this.shortModel();
    const D = abs(4n * A ** 3n + 27n * B ** 2n);
    // every y ≥ 0 with y = 0 or y² | D
    let ys = [1n];
    for (const { prime, exponent } of factor(D)) {
      const next = [];
      for (const y of ys) for (let e = 0, q = 1n; 2 * e <= exponent; e++, q *= prime) next.push(y * q);
      ys = next;
    }
    ys.push(0n);
    const found = [INF];
    for (const y of ys) {
      for (const x of integerCubicRoots(A, B - y * y)) {
        for (const yy of y === 0n ? [0n] : [y, -y]) {
          const P = S.point(x, yy);
          if (S.order(P) !== Infinity) found.push(P);
        }
      }
    }
    const order = found.length;
    const ords = found.map((P) => S.order(P));
    const twos = ords.filter((o) => o === 2n).length;
    let structure;
    let gens;
    if (twos === 3) {
      const m = BigInt(order / 2);
      structure = [2, order / 2];
      const g1 = found[ords.indexOf(m)];
      // a point of order 2 outside the subgroup generated by g1
      const inside = S.mul(m / 2n, g1);
      const g2 = found.find((P, i) => ords[i] === 2n && !S.equal(P, inside));
      gens = [g2, g1];
    } else {
      structure = [order];
      gens = order === 1 ? [] : [found[ords.indexOf(BigInt(order))]];
    }
    const mazur = twos === 3 ? [4, 8, 12, 16].includes(order) : order <= 12 && order !== 11;
    if (!mazur || gens.includes(undefined)) throw new Error(`torsion: a group of order ${order} contradicts Mazur's theorem; this should not happen`);
    const sortKey = (P) => (P.infinity ? -1 : 0);
    const points = found.map(from).sort((P, Q) => sortKey(P) - sortKey(Q) || P.x.compareTo(Q.x) || P.y.compareTo(Q.y));
    return (this._cache.torsion = Object.freeze({ order, structure, points, generators: gens.map(from) }));
  }

  // ---- heights over Q ----

  /** The naive height H(P) = max(|numerator|, |denominator|) of the x-coordinate (1 for the zero). */
  naiveHeight(P) {
    onlyQ(this, "naiveHeight");
    P = this._pt(P);
    if (P.infinity) return 1n;
    const a = abs(P.x.numerator);
    return a > P.x.denominator ? a : P.x.denominator;
  }
  /**
   * An enclosure { lower, upper } (Rationals) of the canonical height
   *     ĥ(P) = lim log H(x(2ⁿP)) / 4ⁿ,
   * taken on the short integral model (the value does not depend on the
   * model). ĥ is 0 exactly for torsion points, is a quadratic form, and is
   * twice the height in Silverman's normalisation. Option `doublings`
   * (default 6) is n: each further doubling divides the width by 4 and
   * multiplies the size of the numbers by 4.
   *
   * The bounds are proved for the curve at hand: for x = p/q the numerator
   * and denominator of x(2Q) are quartic forms Φ, Ψ in (p, q), and identities
   * F·Φ + G·Ψ = c·p⁷ and F′·Φ + G′·Ψ = c′·q⁷ (found by solving a linear
   * system exactly) bound how much can cancel. That gives constants with
   * lo ≤ log H(2Q) − 4·log H(Q) ≤ hi for every Q, hence
   * lo/3 ≤ ĥ(Q) − log H(Q) ≤ hi/3.
   */
  height(P, { doublings = 6 } = {}) {
    onlyQ(this, "height");
    if (!Number.isInteger(doublings) || doublings < 0 || doublings > 12) throw new RangeError("doublings must be a whole number from 0 to 12");
    P = this._pt(P);
    if (this.order(P) !== Infinity) return { lower: R(0), upper: R(0), exact: true };
    const { curve: S, to } = this.shortModel();
    const { lo, hi } = S._heightConstants();
    let Q = to(P);
    for (let i = 0; i < doublings; i++) Q = S.double(Q);
    const [l0, l1] = lnBounds(S.naiveHeight(Q));
    const scale = R(4n ** BigInt(doublings));
    const lower = l0.add(lo.divide(R(3))).divide(scale);
    const upper = l1.add(hi.divide(R(3))).divide(scale);
    return { lower: lower.numerator < 0n ? R(0) : lower, upper, exact: false };
  }
  /** Rational bounds { lo, hi } with lo ≤ log H(2Q) − 4 log H(Q) ≤ hi on this short integral curve. */
  _heightConstants() {
    if (this._cache.hc) return this._cache.hc;
    const A = this.a4.numerator;
    const B = this.a6.numerator;
    // x(2Q) = Φ/Ψ, as forms in (p, q) listed by the power of p
    const Phi = [A * A, -8n * B, -2n * A, 0n, 1n];
    const Psi = [4n * B, 4n * A, 0n, 4n, 0n];
    const norm = (v) => v.reduce((s, c) => s + abs(c), 0n);
    const kappa = norm(Phi) > norm(Psi) ? norm(Phi) : norm(Psi);
    // solve F·Φ + G·Ψ = p^7 and = q^7 for cubic forms F, G (8 unknowns, 8 equations)
    const solveFor = (target) => {
      const M = Array.from({ length: 8 }, (_, k) => {
        const row = [];
        for (let i = 0; i < 4; i++) row.push(R(Phi[k - i] ?? 0n));
        for (let i = 0; i < 4; i++) row.push(R(Psi[k - i] ?? 0n));
        row.push(R(k === target ? 1 : 0));
        return row;
      });
      const sol = solveSquare(M);
      const c = sol.reduce((m, v) => lcm(m, v.denominator), 1n);
      const K = sol.reduce((s, v) => s + abs((v.numerator * c) / v.denominator), 0n);
      return { c, K };
    };
    const s7 = solveFor(7);
    const s0 = solveFor(0);
    // max(|Φ|, |Ψ|) ≥ H⁴·min(c/K), and gcd(Φ, Ψ) divides lcm(c, c′)
    const ratio = [new Rational(s7.c, s7.K), new Rational(s0.c, s0.K)].reduce((a, b) => (a.lessThan(b) ? a : b));
    const rho = ratio.divide(R(lcm(s7.c, s0.c)));
    const lo = lnBounds(rho.numerator)[0].subtract(lnBounds(rho.denominator)[1]);
    const hi = lnBounds(kappa)[1];
    return (this._cache.hc = { lo, hi });
  }
  /**
   * The rational points whose x-coordinate on the short integral model has
   * naive height at most `bound` (default 100), found by search; the zero is
   * not listed. A search, not a description of the group.
   */
  rationalPoints(bound = 100) {
    onlyQ(this, "rationalPoints");
    const N = big(bound, "bound");
    if (N < 1n || N > 100000n) refuse("too-large", "the bound must be from 1 to 100000");
    const { A, B, from } = this.shortModel();
    const out = [];
    for (let d = 1n; d * d <= N; d++) {
      const [d2, d4, d6] = [d * d, d ** 4n, d ** 6n];
      for (let a = -N; a <= N; a++) {
        if (d > 1n && gcd(a, d) !== 1n) continue;
        const rhs = a * a * a + A * a * d4 + B * d6;
        if (rhs < 0n || !isSquare(rhs)) continue;
        const s = isqrt(rhs);
        for (const y of s === 0n ? [0n] : [s, -s]) out.push(from({ x: new Rational(a, d2), y: new Rational(y, d2 * d) }));
      }
    }
    return out.sort((P, Q) => P.x.compareTo(Q.x) || P.y.compareTo(Q.y));
  }

  // ---- curves over F_p ----

  /** The short form y² = x³ + Ax + B of a curve over F_p, p > 3: [A, B]. */
  _shortFp() {
    const p = this.modulus;
    return [mod(-27n * this.c4, p), mod(-54n * this.c6, p)];
  }
  /**
   * The number of points over F_p, counting the point at infinity. Option
   * method: "naive" (a character sum, p up to about 4 million), "bsgs"
   * (Shanks and Mestre, about p^(1/4) steps), "schoof" (polynomial in log p),
   * or "auto". Every method gives a proven count. With method "schoof", an
   * array passed as option log receives { l, t } for each small prime l: the
   * trace modulo l.
   */
  count({ method = "auto", log } = {}) {
    onlyFp(this, "count");
    const p = this.modulus;
    if (method === "auto" && this._cache.count !== undefined) return this._cache.count;
    let N;
    if (p <= 3n) N = BigInt(this.points().length) + 1n;
    else {
      const [A, B] = this._shortFp();
      const m = method === "auto" ? (p < 1n << 22n ? "naive" : p < 1n << 72n ? "bsgs" : "schoof") : method;
      if (m === "naive") {
        if (p >= 1n << 22n) refuse("too-large", "the naive count is limited to p below 2^22; use method \"bsgs\" or \"schoof\"");
        N = countNaive(A, B, Number(p));
      } else if (m === "bsgs") {
        if (p < 1000n) N = countNaive(A, B, Number(p));
        else if (p >= 1n << 90n) refuse("too-large", "baby steps and giant steps is limited to p below 2^90; use method \"schoof\"");
        else N = countBsgs(A, B, p);
      } else if (m === "schoof") N = p + 1n - schoofTrace(A, B, p, { log });
      else throw new RangeError(`unknown method ${method}`);
    }
    return (this._cache.count = N);
  }
  /**
   * The n-th division polynomial of a short curve over F_p (p > 3), as
   * coefficients with the constant term first. Its roots are the
   * x-coordinates of the points P ≠ 0 with n·P = 0 (for even n, those that
   * are not of order 2: the full ψ_n is 2y times this polynomial).
   */
  divisionPolynomial(n) {
    onlyFp(this, "divisionPolynomial");
    if (!this.isShort || this.modulus <= 3n) refuse("bad-argument", "division polynomials are given for short curves y² = x³ + ax + b over F_p, p > 3");
    if (!Number.isInteger(n) || n < 1 || n > 200) throw new RangeError("n must be a whole number from 1 to 200");
    return divisionPolynomials(this.a4, this.a6, this.modulus, Math.max(n, 4))[n];
  }
  /** The trace of Frobenius t = p + 1 − #E(F_p); |t| ≤ 2√p by Hasse's theorem. */
  trace() {
    onlyFp(this, "trace");
    return this.modulus + 1n - this.count();
  }
  /** The points with the given x-coordinate (none, one or two). */
  liftX(x) {
    onlyFp(this, "liftX");
    const p = this.modulus;
    const K = this.K;
    x = K.from(x, "x");
    const rhs = K.add(K.mul(K.add(K.mul(K.add(x, this.a2), x), this.a4), x), this.a6);
    const lin = K.add(K.mul(this.a1, x), this.a3);
    if (p === 2n) {
      return [0n, 1n].filter((y) => (y * y + lin * y) % 2n === rhs).map((y) => Object.freeze({ x, y }));
    }
    // (2y + lin)² = 4·rhs + lin²
    const d = mod(4n * rhs + lin * lin, p);
    const s = sqrtMod(d, p);
    if (s === null) return [];
    const half = invMod(2n, p);
    const ys = s === 0n ? [mod(-lin * half, p)] : [mod((s - lin) * half, p), mod((-s - lin) * half, p)];
    return ys.sort((a, b) => (a < b ? -1 : 1)).map((y) => Object.freeze({ x, y }));
  }
  /** Every affine point, for p up to 10⁶, in order of x then y. */
  points() {
    onlyFp(this, "points");
    if (this.modulus > 1000000n) refuse("too-large", "listing the points is limited to p ≤ 10⁶");
    const out = [];
    for (let x = 0n; x < this.modulus; x++) out.push(...this.liftX(x));
    return out;
  }
  /** A point of the curve chosen by walking x upward from `seed` (reproducible, not uniformly random). */
  randomPoint(seed = 1n) {
    onlyFp(this, "randomPoint");
    for (let x = mod(big(seed, "seed"), this.modulus), i = 0n; i < this.modulus; i++, x = (x + 1n) % this.modulus) {
      const ps = this.liftX(x);
      if (ps.length) return ps[0];
    }
    return INF;
  }
  /**
   * The group E(F_p) as Z/n1 × Z/n2 with n1 | n2: { order, structure: [n1, n2] }
   * (structure [n] when cyclic). Proven: cyclic at once when no d > 1 has
   * d² | #E and d | p − 1; otherwise by the orders of all points, for p ≤ 5000.
   */
  structure() {
    onlyFp(this, "structure");
    const p = this.modulus;
    const N = this.count();
    let canSplit = false;
    for (const { prime, exponent } of factor(N)) if (exponent >= 2 && (p - 1n) % prime === 0n) canSplit = true;
    if (!canSplit) return { order: N, structure: [N] };
    if (p > 5000n) refuse("too-large", "the group may not be cyclic, and deciding that is limited to p ≤ 5000 here");
    let e = 1n;
    for (const P of this.points()) e = lcm(e, this.order(P));
    return { order: N, structure: e === N ? [N] : [N / e, e] };
  }
}

/**
 * An elliptic curve from its coefficients: [a, b] for y² = x³ + ax + b, or
 * [a1, a2, a3, a4, a6] for the long Weierstrass form. Over Q by default
 * (whole numbers, fraction text or Rationals); option mod: p gives the curve
 * over F_p. Singular equations are refused.
 */
export function ellipticCurve(coefficients, { mod: p } = {}) {
  if (!Array.isArray(coefficients) || (coefficients.length !== 2 && coefficients.length !== 5)) {
    throw new TypeError("give [a, b] for y² = x³ + ax + b, or [a1, a2, a3, a4, a6]");
  }
  const cs = coefficients.length === 2 ? [0n, 0n, 0n, ...coefficients] : coefficients;
  if (p === undefined || p === null) return new EllipticCurve(cs.map((c, i) => rat(c, `coefficient ${i + 1}`)), QQ, null);
  const q = big(p, "mod");
  if (q < 2n || !isProbablePrime(q)) refuse("bad-argument", `${q} is not a prime; curves are over prime fields F_p`);
  if (q >= PROVEN_BELOW && !factorProven(q)) refuse("unproven", `${q} is only a probable prime; no certificate was found`);
  const K = FF(q);
  return new EllipticCurve(cs.map((c, i) => K.from(c, `coefficient ${i + 1}`)), K, q);
}
const proven = new Set();
function factorProven(q) {
  if (proven.has(q)) return true;
  const f = factor(q);
  const ok = f.length === 1 && f[0].proven;
  if (ok) proven.add(q);
  return ok;
}

// ---- helpers ----

/** The Legendre symbol (a/p) for an odd prime p, by Euler's criterion on the binary method. */
function legendreSymbol(a, p) {
  a = mod(a, p);
  if (a === 0n) return 0;
  let r = 1n;
  let b = a;
  let e = (p - 1n) / 2n;
  while (e > 0n) {
    if (e & 1n) r = (r * b) % p;
    b = (b * b) % p;
    e >>= 1n;
  }
  return r === 1n ? 1 : -1;
}

/** For a curve with a node mod p: are the tangents at the node defined over F_p? By direct search (small p). */
export function nodeIsSplit(E, p) {
  const [a1, a2, a3, a4] = [E.a1, E.a2, E.a3, E.a4].map((c) => mod(c.numerator, p));
  const a6 = mod(E.a6.numerator, p);
  for (let x = 0n; x < p; x++) {
    for (let y = 0n; y < p; y++) {
      const F = mod(y * y + a1 * x * y + a3 * y - x * x * x - a2 * x * x - a4 * x - a6, p);
      const Fx = mod(a1 * y - 3n * x * x - 2n * a2 * x - a4, p);
      const Fy = mod(2n * y + a1 * x + a3, p);
      if (F !== 0n || Fx !== 0n || Fy !== 0n) continue;
      // tangent cone at the singular point: T² + a1·T − (3x + a2) = 0 for the slope T
      for (let T = 0n; T < p; T++) if (mod(T * T + a1 * T - 3n * x - a2, p) === 0n) return true;
      return false;
    }
  }
  throw new Error("nodeIsSplit: no singular point");
}

/** The reduced integral model with invariants (c4, c6), or null when there is none. */
function reducedFrom(c4, c6) {
  const b2 = { 0: 0n, 1: 1n, 4: 4n, 5: 5n, 8: -4n, 9: -3n }[Number(mod(-c6, 12n))];
  if (b2 === undefined) return null;
  if (mod(b2 * b2 - c4, 24n) !== 0n) return null;
  const b4 = (b2 * b2 - c4) / 24n;
  const n6 = -b2 * b2 * b2 + 36n * b2 * b4 - c6;
  if (mod(n6, 216n) !== 0n) return null;
  const b6 = n6 / 216n;
  const a1 = mod(b2, 2n);
  const a2 = (b2 - a1) / 4n;
  const a3 = mod(b6, 2n);
  if (mod(b6 - a3, 4n) !== 0n || mod(b4 - a1 * a3, 2n) !== 0n) return null;
  const a6 = (b6 - a3) / 4n;
  const a4 = (b4 - a1 * a3) / 2n;
  // check: these coefficients do have the invariants asked for
  const B2 = a1 * a1 + 4n * a2;
  const B4 = 2n * a4 + a1 * a3;
  const B6 = a3 * a3 + 4n * a6;
  if (B2 * B2 - 24n * B4 !== c4 || -B2 * B2 * B2 + 36n * B2 * B4 - 216n * B6 !== c6) return null;
  if (c4 ** 3n === c6 ** 2n) return null;
  return [a1, a2, a3, a4, a6];
}

/** The whole roots of x³ + Ax + C. */
export function integerCubicRoots(A, C) {
  const f = (x) => x * x * x + A * x + C;
  const Rb = 2n + (abs(A) > abs(C) ? abs(A) : abs(C));
  const out = [];
  const search = (lo, hi, increasing) => {
    if (lo > hi) return;
    while (lo < hi) {
      const mid = (lo + hi) >> 1n;
      const v = f(mid);
      if ((v < 0n) === increasing) lo = mid + 1n;
      else hi = mid;
    }
    if (f(lo) === 0n) out.push(lo);
  };
  if (A >= 0n) search(-Rb, Rb, true);
  else {
    const c = isqrt(-A / 3n);
    search(-Rb, -c - 1n, true);
    // on −c … c the cubic does not increase; search its negative
    let [lo, hi] = [-c, c];
    while (lo < hi) {
      const mid = (lo + hi) >> 1n;
      if (f(mid) > 0n) lo = mid + 1n;
      else hi = mid;
    }
    if (f(lo) === 0n) out.push(lo);
    search(c + 1n, Rb, true);
  }
  return [...new Set(out)].sort((a, b) => (a < b ? -1 : 1));
}

/** Solve a square system given as an augmented matrix of Rationals. */
function solveSquare(M) {
  const n = M.length;
  for (let c = 0; c < n; c++) {
    const piv = M.findIndex((row, r) => r >= c && row[c].numerator !== 0n);
    if (piv < 0) throw new Error("height: the linear system is singular; this should not happen");
    [M[c], M[piv]] = [M[piv], M[c]];
    const inv = M[c][c].reciprocal();
    M[c] = M[c].map((v) => v.multiply(inv));
    for (let r = 0; r < n; r++) {
      if (r === c || M[r][c].numerator === 0n) continue;
      const f = M[r][c];
      M[r] = M[r].map((v, k) => v.subtract(f.multiply(M[c][k])));
    }
  }
  return M.map((row) => row[n]);
}

/** Rational bounds [lo, hi] on the natural logarithm of a whole number m ≥ 1. */
export function lnBounds(m) {
  if (m < 1n) throw new RangeError("lnBounds needs m ≥ 1");
  const PREC = 96n;
  const one = 1n << PREC;
  const TERMS = 48;
  const SLACK = BigInt(TERMS * TERMS + 4 * TERMS + 8);
  // bounds on ln f for f = fnum/2^PREC in [1, 2], from ln f = 2·atanh((f − 1)/(f + 1))
  const lnFix = (fLo, fHi) => {
    const series = (z) => {
      let sum = 0n;
      let term = z;
      const z2 = (z * z) >> PREC;
      for (let i = 0; i < TERMS; i++) {
        sum += term / BigInt(2 * i + 1);
        term = (term * z2) >> PREC;
      }
      return sum;
    };
    const zLo = ((fLo - one) << PREC) / (fLo + one);
    const zHi = ((fHi - one) << PREC) / (fHi + one) + 1n;
    // truncation only lowers the sum; the tail is below (1/3)^(2·TERMS)
    return [2n * series(zLo), 2n * (series(zHi) + SLACK)];
  };
  const [ln2Lo, ln2Hi] = lnFix(2n * one, 2n * one);
  const k = BigInt(m.toString(2).length - 1);
  // f = m / 2^k in [1, 2): keep PREC bits of it
  const fLo = k <= PREC ? m << (PREC - k) : m >> (k - PREC);
  const fHi = k <= PREC ? fLo : fLo + 1n;
  const [a, b] = fLo === one && fHi === one ? [0n, 0n] : lnFix(fLo, fHi > 2n * one ? 2n * one : fHi);
  return [new Rational(k * ln2Lo + a, one), new Rational(k * ln2Hi + b, one)];
}

/** #E(F_p) for y² = x³ + Ax + B by the character sum p + 1 + Σ χ(x³ + Ax + B), p < 2^22. */
function countNaive(A, B, p) {
  const a = Number(A);
  const b = Number(B);
  const chi = new Int8Array(p).fill(-1);
  chi[0] = 0;
  for (let i = 1; i <= (p - 1) / 2; i++) chi[(i * i) % p] = 1;
  let s = 0;
  for (let x = 0; x < p; x++) s += chi[(((((x * x) % p) * x) % p) + ((a * x) % p) + b) % p];
  return BigInt(p + 1 + s);
}

/**
 * #E(F_p) by Shanks and Mestre. The order of a point P divides #E, and #E
 * lies in the Hasse interval [p + 1 − 2√p, p + 1 + 2√p]. Baby steps and giant
 * steps find a multiple of the order of P in that interval, from which the
 * exact order follows. Orders of several points are combined by lcm, on the
 * curve and on its quadratic twist (whose count is 2p + 2 − #E); as soon as
 * one of the two lcms has a single multiple in its interval, the count is
 * proved. Mestre showed this must happen for p > 229.
 */
function countBsgs(A, B, p) {
  const T = isqrt(4n * p);
  let g = 2n;
  while (legendreSymbol(g, p) !== -1) g++;
  const curves = [
    { A, B, L: 1n },
    { A: (A * g * g) % p, B: (((B * g * g) % p) * g) % p, L: 1n },
  ];
  const unique = (L) => {
    const lo = p + 1n - T;
    const first = ((lo + L - 1n) / L) * L;
    return first + L > p + 1n + T ? first : null;
  };
  let x = 0n;
  for (let tries = 0; tries < 400; tries++) {
    for (let side = 0; side < 2; side++) {
      const C = curves[side];
      let P = null;
      for (; P === null; x++) {
        const rhs = mod(x * x * x + C.A * x + C.B, p);
        if (rhs === 0n) continue; // points of order 2 add little
        const y = sqrtMod(rhs, p);
        if (y !== null) P = [x, y];
      }
      C.L = lcm(C.L, pointOrderBsgs(P, C.A, p, T));
      const N = unique(C.L);
      if (N !== null) return side === 0 ? N : 2n * p + 2n - N;
    }
  }
  refuse("too-hard", "baby steps and giant steps did not settle the count; use method \"schoof\"");
}
/** The exact order of P on y² = x³ + Ax + B over F_p, from a multiple of it in the Hasse interval. */
function pointOrderBsgs(P, A, p, T) {
  const s = isqrt(T) + 1n;
  // baby steps: x(jP) for j = 1 … s
  const baby = new Map();
  let J = P;
  for (let j = 1n; j <= s; j++) {
    if (J === null) return reduceOrder(j, P, A, p);
    if (!baby.has(J[0])) baby.set(J[0], j);
    J = ecAdd(J, P, A, p);
  }
  const G = ecMul(2n * s + 1n, P, A, p);
  const I = T / (2n * s + 1n) + 1n;
  let Rr = ecAdd(ecMul(p + 1n, P, A, p), ecMul(I, G === null ? null : [G[0], mod(-G[1], p)], A, p), A, p);
  for (let i = -I; i <= I; i++) {
    const base = p + 1n + i * (2n * s + 1n);
    if (Rr === null) return reduceOrder(base, P, A, p);
    const j = baby.get(Rr[0]);
    if (j !== undefined) {
      for (const M of [base - j, base + j]) if (M > 0n && ecMul(M, P, A, p) === null) return reduceOrder(M, P, A, p);
    }
    Rr = ecAdd(Rr, G, A, p);
  }
  throw new Error("pointOrderBsgs: no multiple of the order in the Hasse interval; this should not happen");
}
function reduceOrder(M, P, A, p) {
  for (const { prime } of factor(M)) while (M % prime === 0n && ecMul(M / prime, P, A, p) === null) M /= prime;
  return M;
}
