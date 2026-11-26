// Phase 4: elliptic curves, the elliptic curve method, and elliptic curve
// primality certificates. Each algorithm is checked against a slower one that
// shares no code with it: brute-force point lists, counts by three methods,
// reduction modulo primes, and the independent certificate verifier.

import { describe, expect, test } from "bun:test";
import { Rational } from "@ratmath/core";
import * as N from "../index.js";
import { nodeIsSplit, lnBounds, integerCubicRoots } from "../src/elliptic.js";
import { schoofTrace } from "../src/schoof.js";
import { cornacchia4 } from "../src/ecpp.js";

const E = N.ellipticCurve;
let seed = 12345n;
const rand = (n) => {
  seed = (seed * 6364136223846793005n + 1442695040888963407n) % (1n << 64n);
  return (seed >> 11n) % n;
};
const gcd = (a, b) => {
  while (b) [a, b] = [b, a % b];
  return a < 0n ? -a : a;
};
const Q = (s) => new Rational(s);

describe("the group law", () => {
  test("multiples of (0, 0) on y² + y = x³ − x", () => {
    const C = E([0, 0, 1, -1, 0]);
    expect(String(C)).toBe("y² + y = x³ − x");
    expect(C.discriminant.toString()).toBe("37");
    expect(C.j.toString()).toBe("110592/37");
    const P = C.point(0, 0);
    const xs = [1, 2, 3, 4, 5, 6, 7].map((k) => C.mul(k, P)).map((R) => `${R.x},${R.y}`);
    expect(xs).toEqual(["0,0", "1,0", "-1,-1", "2,-3", "1/4,-5/8", "6,14", "-5/9,8/27"]);
    expect(C.mul(-3, P)).toEqual(C.neg(C.mul(3, P)));
    expect(C.add(C.mul(5, P), C.mul(-5, P)).infinity).toBe(true);
    expect(C.order(P)).toBe(Infinity);
    expect(() => C.point(1, 1)).toThrow(/not on/);
    expect(C.has([1, 0])).toBe(true);
    expect(C.has([1, 1])).toBe(false);
  });
  test("singular cubics and bad input are refused", () => {
    expect(() => E([0, 0])).toThrow(/node or a cusp/);
    expect(() => E([-3, 2])).toThrow(/node or a cusp/);
    expect(() => E([1, 2, 3])).toThrow(TypeError);
    expect(() => E([1.5, 2])).toThrow(TypeError);
    expect(() => E([1, 1], { mod: 15 })).toThrow(/not a prime/);
    expect(() => E([0, 1], { mod: 3 })).toThrow(/node or a cusp/);
  });
  test("over F_p the points form a group whose order is the count (all Weierstrass forms, including p = 2 and 3)", () => {
    for (const p of [2n, 3n, 5n, 7n, 11n, 13n, 31n]) {
      for (let trial = 0; trial < 12; trial++) {
        let C;
        try {
          C = E([rand(p), rand(p), rand(p), rand(p), rand(p)], { mod: p });
        } catch {
          continue;
        }
        const pts = [C.infinity, ...C.points()];
        const Np = C.count();
        expect(BigInt(pts.length)).toBe(Np);
        // Hasse
        expect((C.trace()) ** 2n <= 4n * p).toBe(true);
        for (const P of pts) {
          expect(C.has(P)).toBe(true);
          expect(C.mul(Np, P).infinity).toBe(true);
          expect(Np % C.order(P)).toBe(0n);
          expect(C.add(P, C.neg(P)).infinity).toBe(true);
        }
        for (let i = 0; i < 40; i++) {
          const [P, R, S] = [pts[Number(rand(Np))], pts[Number(rand(Np))], pts[Number(rand(Np))]];
          expect(C.equal(C.add(P, R), C.add(R, P))).toBe(true);
          expect(C.equal(C.add(C.add(P, R), S), C.add(P, C.add(R, S)))).toBe(true);
        }
        // the structure Z/n1 × Z/n2 has n1 | n2, n1 | p − 1 and the right exponent
        const st = C.structure();
        const exps = pts.map((P) => C.order(P)).reduce((a, b) => (a * b) / gcd(a, b), 1n);
        expect(st.structure[st.structure.length - 1]).toBe(exps);
        expect(st.structure.reduce((a, b) => a * b, 1n)).toBe(Np);
      }
    }
  });
  test("associativity over Q on a curve of rank 2", () => {
    const C = E([0, 1, 1, -2, 0]);
    const [P, R] = [C.point(0, 0), C.point(-1, 1)];
    const lhs = C.add(C.add(C.mul(3, P), R), C.mul(2, R));
    expect(C.equal(lhs, C.add(C.mul(3, P), C.mul(3, R)))).toBe(true);
    expect(C.equal(C.sub(C.add(P, R), R), P)).toBe(true);
  });
});

describe("counting points", () => {
  test("Schoof agrees with the character sum on every curve over small fields", () => {
    for (const p of [5n, 7n, 11n, 13n, 17n, 19n, 23n]) {
      for (let A = 0n; A < p; A++) {
        for (let B = 0n; B < p; B++) {
          if ((4n * A ** 3n + 27n * B ** 2n) % p === 0n) continue;
          const C = E([A, B], { mod: p });
          const naive = BigInt(C.points().length) + 1n;
          expect(C.count()).toBe(naive);
          expect(p + 1n - schoofTrace(A, B, p)).toBe(naive);
        }
      }
    }
  });
  test("naive, baby-step giant-step and Schoof agree on larger fields", () => {
    for (const p of [1009n, 10007n, 65537n, 1000003n]) {
      for (let i = 0; i < 4; i++) {
        const C = E([rand(p), rand(p) + 1n], { mod: p });
        const n = C.count({ method: "naive" });
        expect(C.count({ method: "bsgs" })).toBe(n);
        expect(C.count({ method: "schoof" })).toBe(n);
      }
    }
    // long Weierstrass forms are counted through their short form
    const C = E([1, 2, 3, 4, 5], { mod: 10007 });
    expect(C.count()).toBe(BigInt(C.points().length) + 1n);
    const big = E([2, 3], { mod: 1000000007 });
    expect(big.count()).toBe(1000004178n);
    const log = [];
    expect(big.count({ method: "schoof", log })).toBe(1000004178n);
    for (const { l, t } of log) expect((((big.trace() - BigInt(t)) % BigInt(l)) + BigInt(l)) % BigInt(l)).toBe(0n);
    const P = big.randomPoint(7);
    expect(big.mul(big.count(), P).infinity).toBe(true);
    expect(big.mul(big.order(P), P).infinity).toBe(true);
  });
  test("a 61-bit field by baby steps and giant steps", () => {
    const C = E([2, 3], { mod: 2n ** 61n - 1n });
    const n = C.count();
    expect(n).toBe(2305843011631544440n);
    expect(C.mul(n, C.randomPoint(5)).infinity).toBe(true);
  });
  test("division polynomials vanish exactly on the points killed by n", () => {
    const p = 103n;
    const C = E([5, 7], { mod: p });
    const pts = C.points();
    for (const n of [2, 3, 4, 5, 6, 7, 9]) {
      const f = C.divisionPolynomial(n);
      expect(f.length - 1).toBe(n % 2 ? (n * n - 1) / 2 : (n * n - 4) / 2);
      for (const P of pts) {
        const value = f.reduceRight((acc, c) => (acc * P.x + c) % p, 0n);
        const killedOdd = C.mul(n, P).infinity === true;
        const killed = killedOdd && C.mul(2, P).infinity !== true;
        expect(value === 0n).toBe(n % 2 ? killedOdd : killed);
      }
    }
  });
});

describe("curves over Q", () => {
  const known = [
    [[0, 1], [6]],
    [[-1, 0], [2, 2]],
    [[0, -2], [1]],
    [[4, 0], [4]],
    [[0, -1, 1, 0, 0], [5]],
    [[1, 0, 1, 4, -6], [6]],
    [[1, 1, 1, -10, -10], [2, 4]],
    [[-43, 166], [7]],
    [[1, -1, 1, -14, 29], [9]],
    [[1, 0, 1, -19, 26], [2, 6]],
    [[1, 0, 0, -1070, 7812], [2, 8]],
    [[0, 0, 1, -1, 0], [1]],
  ];
  test("torsion subgroups of known curves", () => {
    for (const [cs, structure] of known) {
      const C = E(cs);
      const T = C.torsion();
      expect(T.structure).toEqual(structure);
      expect(T.order).toBe(structure.reduce((a, b) => a * b, 1));
      expect(T.points.length).toBe(T.order);
    }
  });
  test("torsion is a group, is generated by its generators, and embeds in E(F_p)", () => {
    const curves = known.map(([cs]) => cs);
    for (let i = 0; i < 40; i++) curves.push([Number(rand(9n)) - 4, Number(rand(21n)) - 10]);
    for (const a of [0, 1]) for (const b of [-1, 0, 1]) for (let c = -3; c <= 3; c++) curves.push([a, b, 1, c, 0], [a, b, 0, c, 1]);
    let nontrivial = 0;
    for (const cs of curves) {
      let C;
      try {
        C = E(cs);
      } catch {
        continue;
      }
      const T = C.torsion();
      if (T.order > 1) nontrivial++;
      const key = (P) => (P.infinity ? "O" : `${P.x},${P.y}`);
      const set = new Set(T.points.map(key));
      for (const P of T.points) for (const R of T.points) expect(set.has(key(C.add(P, R)))).toBe(true);
      // the generators give every point
      const span = new Set();
      const [g1, g2] = T.generators;
      for (let i = 0; i < 16; i++) for (let j = 0; j < 16; j++) span.add(key(C.add(C.mul(i, g1 ?? C.infinity), C.mul(j, g2 ?? C.infinity))));
      expect(span.size).toBe(T.order);
      // reduction is injective on torsion at odd primes of good reduction
      let g = 0n;
      const disc = C.discriminant;
      for (const p of N.primes(120)) {
        if (p < 3n || disc.numerator % p === 0n || disc.denominator % p === 0n) continue;
        g = gcd(g, C.modulo(p).count());
      }
      expect(g % BigInt(T.order)).toBe(0n);
    }
    expect(nontrivial).toBeGreaterThan(20);
  });
  test("a torsion point found on a model with fractions", () => {
    const C = E(["1/4", 0]); // y² = x³ + x/4, isomorphic to y² = x³ + 4x
    expect(C.torsion().structure).toEqual([4]);
    expect(C.isIntegral).toBe(false);
    expect(C.shortModel().curve.isIntegral).toBe(true);
  });
  test("whole roots of cubics", () => {
    for (let i = 0; i < 300; i++) {
      const [r, s, t] = [rand(41n) - 20n, rand(41n) - 20n, 0n];
      // (x − r)(x − s)(x + r + s) has no x² term
      const u = -r - s;
      const A = r * s + r * u + s * u;
      const C = -r * s * u + t;
      expect(integerCubicRoots(A, C)).toEqual([...new Set([r, s, u])].sort((a, b) => (a < b ? -1 : 1)));
      const shifted = integerCubicRoots(A, C + 1n);
      for (const x of shifted) expect(x ** 3n + A * x + C + 1n).toBe(0n);
      for (let x = -45n; x <= 45n; x++) if (x ** 3n + A * x + C + 1n === 0n) expect(shifted.includes(x)).toBe(true);
    }
  });
  test("minimal models: a change of variables is undone exactly", () => {
    expect(String(E([-432, 8208]).minimalModel().curve)).toBe("y² + y = x³ − x²");
    let done = 0;
    for (let i = 0; i < 150; i++) {
      let C;
      try {
        C = E([rand(2n), rand(3n) - 1n, rand(2n), rand(13n) - 6n, rand(13n) - 6n]);
      } catch {
        continue;
      }
      const M = C.minimalModel();
      // a reduced model with small coefficients is usually minimal already; start from the minimal one
      const m = M.curve;
      const [u, r, s, t] = [rand(6n) + 1n, rand(7n) - 3n, rand(5n) - 2n, rand(7n) - 3n];
      // x = u²x′ + r, y = u³y′ + u²s·x′ + t applied backwards: the curve with a_i′·u^i = …
      const [a1, a2, a3, a4, a6] = m.coefficients.map((c) => c.numerator);
      const b1 = a1 + 2n * s;
      const b2 = a2 - s * a1 + 3n * r - s * s;
      const b3 = a3 + r * a1 + 2n * t;
      const b4 = a4 - s * a3 + 2n * r * a2 - (t + r * s) * a1 + 3n * r * r - 2n * s * t;
      const b6 = a6 + r * a4 + r * r * a2 + r * r * r - t * a3 - t * t - r * t * a1;
      const big = E([b1 * u, b2 * u ** 2n, b3 * u ** 3n, b4 * u ** 4n, b6 * u ** 6n]);
      expect(big.j.equals(m.j)).toBe(true);
      const back = big.minimalModel();
      expect(back.curve.coefficients.map(String)).toEqual(m.coefficients.map(String));
      expect(big.discriminant.divide(back.curve.discriminant).equals(back.u.pow(12))).toBe(true);
      // the maps carry points both ways
      for (const P of big.rationalPoints(30).slice(0, 4)) {
        const Pm = back.to(P);
        expect(back.curve.has(Pm)).toBe(true);
        expect(big.equal(back.from(Pm), P)).toBe(true);
      }
      done++;
    }
    expect(done).toBeGreaterThan(100);
    // fractions in, whole numbers out
    const F = E(["1/4", "1/8"]);
    expect(F.minimalModel().curve.isIntegral).toBe(true);
    expect(F.minimalModel().curve.j.equals(F.j)).toBe(true);
  });
  test("reduction types, against a direct look at the node", () => {
    let seen = 0;
    for (let a2 = -2; a2 <= 2; a2++) {
      for (let a4 = -4; a4 <= 4; a4++) {
        for (let a6 = -4; a6 <= 4; a6++) {
          for (const [a1, a3] of [[0, 0], [1, 0], [0, 1], [1, 1]]) {
            let C;
            try {
              C = E([a1, a2, a3, a4, a6]);
            } catch {
              continue;
            }
            const M = C.minimalModel().curve;
            for (const p of C.badPrimes()) {
              const r = C.reduction(p);
              expect(r.discriminantValuation).toBeGreaterThan(0);
              if (r.type === "multiplicative" && p < 60n) {
                expect(r.split).toBe(nodeIsSplit(M, p));
                seen++;
              }
              // additive exactly when the reduced cubic has a cusp: c4 ≡ 0
              expect(r.type === "additive").toBe(M.c4.numerator % p === 0n);
            }
            expect(C.reduction(1000003n).type === "good").toBe(M.discriminant.numerator % 1000003n !== 0n);
          }
        }
      }
    }
    expect(seen).toBeGreaterThan(500);
  });
  test("conductors of the first curves in the tables", () => {
    expect(E([0, -1, 1, 0, 0]).conductor()).toBe(11n);
    expect(E([1, 0, 1, 4, -6]).conductor()).toBe(14n);
    expect(E([1, 1, 1, -10, -10]).conductor()).toBe(15n);
    expect(E([0, 0, 1, -1, 0]).conductor()).toBe(37n);
    expect(E([0, 0, 1, -1, 0]).reduction(37)).toMatchObject({ type: "multiplicative", split: false });
    expect(E([0, -1, 1, 0, 0]).reduction(11)).toMatchObject({ type: "multiplicative", split: true });
    expect(E([-1, 0]).reduction(2).type).toBe("additive");
    expect(() => E([-1, 0]).conductor()).toThrow(/Tate/);
    expect(E([0, 0, 0, 0, 5 ** 2]).reduction(5)).toMatchObject({ type: "additive", conductorExponent: 2 });
  });
  test("points of small height", () => {
    const C = E([0, 0, 1, -1, 0]);
    const pts = C.rationalPoints(30);
    for (const P of pts) expect(C.has(P)).toBe(true);
    const xs = new Set(pts.map((P) => P.x.toString()));
    for (const x of ["0", "1", "-1", "2", "1/4", "6"]) expect(xs.has(x)).toBe(true);
  });
});

describe("heights", () => {
  test("logarithm bounds enclose the logarithm", () => {
    for (const m of [1n, 2n, 3n, 10n, 12345678901234567890n, 2n ** 300n + 1n, 2n ** 96n, 2n ** 97n - 1n]) {
      const [lo, hi] = lnBounds(m);
      expect(lo.lessThanOrEqual(hi)).toBe(true);
      const f = m < 2n ** 900n ? Math.log(Number(m)) : null;
      expect(lo.toNumber()).toBeLessThanOrEqual(f + 1e-12);
      expect(hi.toNumber()).toBeGreaterThanOrEqual(f - 1e-12);
      expect(hi.subtract(lo).toNumber()).toBeLessThan(1e-20);
    }
    expect(lnBounds(1n).map(String)).toEqual(["0", "0"]);
  });
  test("the canonical height of the generator of 37a", () => {
    const C = E([0, 0, 1, -1, 0]);
    const P = C.point(0, 0);
    const known = 0.0511114082399688; // the regulator of 37a1 (Cremona's tables)
    let width = Infinity;
    for (const doublings of [3, 5, 7]) {
      const h = C.height(P, { doublings });
      expect(h.lower.toNumber()).toBeLessThanOrEqual(known);
      expect(h.upper.toNumber()).toBeGreaterThanOrEqual(known);
      const w = h.upper.subtract(h.lower).toNumber();
      expect(w).toBeLessThan(width / 8);
      width = w;
    }
    expect(width).toBeLessThan(1e-3);
    // ĥ(kP) = k²·ĥ(P)
    const h3 = C.height(C.mul(3, P), { doublings: 5 });
    expect(h3.lower.toNumber()).toBeLessThanOrEqual(9 * known);
    expect(h3.upper.toNumber()).toBeGreaterThanOrEqual(9 * known);
  });
  test("torsion has height 0, and the parallelogram law holds within the enclosures", () => {
    const T = E([0, 1]);
    expect(T.height(T.point(2, 3))).toMatchObject({ exact: true });
    expect(T.height(T.point(2, 3)).upper.toString()).toBe("0");
    const C = E([0, 1, 1, -2, 0]);
    const [P, R] = [C.point(0, 0), C.point(-1, 1)];
    const h = (X) => C.height(X, { doublings: 6 });
    const [hp, hr, hs, hd] = [h(P), h(R), h(C.add(P, R)), h(C.sub(P, R))];
    // ĥ(P + R) + ĥ(P − R) = 2ĥ(P) + 2ĥ(R)
    const lhsLo = hs.lower.add(hd.lower).toNumber();
    const lhsHi = hs.upper.add(hd.upper).toNumber();
    const rhsLo = 2 * (hp.lower.toNumber() + hr.lower.toNumber());
    const rhsHi = 2 * (hp.upper.toNumber() + hr.upper.toNumber());
    expect(lhsLo).toBeLessThanOrEqual(rhsHi);
    expect(rhsLo).toBeLessThanOrEqual(lhsHi);
    expect(hp.lower.toNumber()).toBeGreaterThan(0);
  });
  test("the doubling bounds hold on every point tried", () => {
    let tried = 0;
    for (const cs of [[0, 0, 1, -1, 0], [0, 1, 1, -2, 0], [-2, 2], [1, 0, 0, -45, 81], [3, -5], [-7, 6]]) {
      const C = E(cs);
      const S = C.shortModel().curve;
      const { lo, hi } = S._heightConstants();
      for (const P0 of S.rationalPoints(200)) {
        let P = P0;
        for (let k = 0; k < 4; k++) {
          const D = S.double(P);
          if (D.infinity) break;
          const [a0, a1] = lnBounds(S.naiveHeight(P));
          const [b0, b1] = lnBounds(S.naiveHeight(D));
          // lo ≤ ln H(2P) − 4 ln H(P) ≤ hi
          expect(b1.subtract(a0.multiply(Q("4"))).greaterThanOrEqual(lo)).toBe(true);
          expect(b0.subtract(a1.multiply(Q("4"))).lessThanOrEqual(hi)).toBe(true);
          P = D;
          tried++;
        }
      }
    }
    expect(tried).toBeGreaterThan(40);
  });
});

describe("the elliptic curve method", () => {
  test("F₇ = 2^128 + 1 (Morrison and Brillhart, 1970, by continued fractions)", () => {
    const F7 = 2n ** 128n + 1n;
    const r = N.ecm(F7, { B1: 2000, curves: 60, details: true });
    expect(r.factor).toBe(59649589127497217n);
    // the same curve again, by its σ
    expect(N.ecmCurve(F7, r.sigma, { B1: 2000 }).factor).toBe(r.factor);
    expect(N.formatFactors(N.factor(F7))).toBe("59649589127497217 · 5704689200685129054721");
    expect(N.factor(F7).every((f) => f.proven)).toBe(true);
  });
  test("factors are found by size, and are always proper factors", () => {
    for (const [p, q] of [[1000003n, 2n ** 89n - 1n], [1000000007n, 2n ** 107n - 1n], [100000000003n, 2n ** 127n - 1n], [10n ** 15n + 37n, 10n ** 18n + 9n]]) {
      const n = p * q;
      const d = N.ecm(n, { curves: 100 });
      expect(d === p || d === q).toBe(true);
    }
    expect(N.ecm(2n ** 89n - 1n, { curves: 5 })).toBe(null);
    expect(N.ecm(91n * 10n)).toBe(2n);
    expect(N.ecmCurve(1000003n * 1000033n, 6, { B1: 2, B2: 2 }).factor).toBe(null);
    expect(() => N.ecmCurve(15n, 3)).toThrow(RangeError);
  });
  test("a curve finds p exactly when its order mod p is smooth", () => {
    // for the Suyama curve with parameter σ, compute the group order mod p independently and compare
    const p = 1000003n;
    const q = 2n ** 61n - 1n;
    const [B1, B2] = [50, 400];
    let hits = 0;
    for (let s = 6n; s < 60n; s++) {
      const u = (s * s - 5n) % p;
      const v = (4n * s) % p;
      const inv = (a) => N.modPow(a, p - 2n, p);
      const den = (16n * u ** 3n * v) % p;
      if (den === 0n) continue;
      const A = ((((v - u + p) ** 3n * (3n * u + v)) % p) * inv(den) * 4n - 2n + p) % p;
      let C;
      try {
        C = E([0n, A, 0n, 1n, 0n], { mod: p });
      } catch {
        continue;
      }
      // the x-only ladder works on the curve or on its twist, whichever holds the starting point
      const x0 = (u ** 3n * inv(v ** 3n % p)) % p;
      const order = C.liftX(x0).length ? C.count() : 2n * p + 2n - C.count();
      const fs = N.factor(order);
      // found in stage 1 or 2 when the order is B1-powersmooth apart from at most one prime ≤ B2
      const bigOnes = fs.filter((f) => f.prime ** BigInt(f.exponent) > BigInt(B1));
      const sure = bigOnes.length === 0 || (bigOnes.length === 1 && bigOnes[0].exponent === 1 && bigOnes[0].prime <= BigInt(B2));
      const r = N.ecmCurve(p * q, s, { B1, B2 });
      if (sure) {
        expect(r.factor).toBe(p);
        hits++;
      }
      if (r.factor !== null) expect(r.factor === p || r.factor === q).toBe(true);
    }
    expect(hits).toBeGreaterThan(3);
  });
});

describe("elliptic curve primality certificates", () => {
  const kinds = (c, acc = {}) => {
    acc[c.kind] = (acc[c.kind] ?? 0) + 1;
    if (c.certificate) kinds(c.certificate, acc);
    for (const f of c.factors ?? []) kinds(f.certificate, acc);
    return acc;
  };
  const firstElliptic = (c) => {
    if (c.kind === "elliptic") return c;
    for (const f of c.factors ?? []) {
      const e = firstElliptic(f.certificate);
      if (e) return e;
    }
    return c.certificate ? firstElliptic(c.certificate) : null;
  };
  test("class polynomials", () => {
    expect(N.classPolynomial(-3)).toEqual([0n, 1n]);
    expect(N.classPolynomial(-4)).toEqual([-1728n, 1n]);
    expect(N.classPolynomial(-163)).toEqual([262537412640768000n, 1n]);
    expect(N.classPolynomial(-15)).toEqual([-121287375n, 191025n, 1n]);
    expect(N.classPolynomial(-23)).toEqual([12771880859375n, -5151296875n, 3491750n, 1n]);
    for (let D = -3; D > -400; D--) {
      if ((((D % 4) + 4) % 4) > 1) continue;
      const H = N.classPolynomial(D);
      expect(H.length - 1).toBe(Number(N.classNumber(D)));
    }
    // H_D splits into linear factors modulo a prime p with 4p = t² + |D|v², and its roots give curves with p + 1 ± t points
    const D = -47;
    const H = N.classPolynomial(D);
    let checked = 0;
    for (const p of N.primesBetween(1000, 6000)) {
      const tv = cornacchia4(p, D);
      const roots = [];
      for (let x = 0n; x < p; x++) if (N.polyMod(H, x, p) === 0n) roots.push(x);
      if (!tv) continue;
      expect(roots.length).toBe(5);
      const j = roots[0];
      const c = (j * N.modInverse(1728n - j + p, p)) % p;
      const t = E([3n * c, 2n * c], { mod: p }).trace();
      expect(t === tv[0] || t === -tv[0]).toBe(true);
      checked++;
    }
    expect(checked).toBeGreaterThan(5);
  });
  test("primes that n − 1 cannot certify", () => {
    for (const [bits, offset] of [[100, 12345n], [128, 0n], [160, 999n], [200, 5n], [256, 77n], [384, 3n]]) {
      let n = ((1n << BigInt(bits)) + offset * 7919n) | 1n;
      while (!N.isProbablePrime(n)) n += 2n;
      // a small budget for n − 1, so that the elliptic step does the work
      const r = N.provePrime(n, { budget: 256 });
      expect(r.proven).toBe(true);
      expect(N.verifyCertificate(r.certificate)).toEqual({ valid: true });
      expect(N.provePrime(n, { budget: 256, elliptic: false }).proven === false || !kinds(r.certificate).elliptic).toBe(true);
    }
    const p = 2n ** 255n - 19n;
    const r = N.provePrime(p, { budget: 1 << 12 });
    expect(r.proven).toBe(true);
    expect(kinds(r.certificate).elliptic).toBeGreaterThan(0);
    expect(N.verifyCertificate(r.certificate).valid).toBe(true);
    expect(N.isPrime(10n ** 60n + 7n)).toBe(true);
    expect(N.nextPrime(10n ** 40n)).toBe(10n ** 40n + 121n);
  });
  test("a tampered elliptic certificate is rejected", () => {
    const r = N.provePrime(2n ** 255n - 19n, { budget: 1 << 12 });
    const e = firstElliptic(r.certificate);
    expect(e).not.toBe(null);
    expect(N.verifyCertificate(e).valid).toBe(true);
    const bad = (change) => N.verifyCertificate({ ...e, ...change }).valid;
    expect(bad({ x: (e.x + 1n) % e.n })).toBe(false);
    expect(bad({ y: (e.n - e.y) % e.n })).toBe(true); // −P is as good as P
    expect(bad({ a: (e.a + 1n) % e.n })).toBe(false);
    expect(bad({ m: e.m + e.q })).toBe(false);
    expect(bad({ q: e.q + 2n })).toBe(false);
    expect(bad({ certificate: { n: e.q, kind: "small" } })).toBe(false);
    expect(bad({ certificate: undefined })).toBe(false);
    expect(bad({ n: e.n + 2n })).toBe(false);
    // a small prime q: the point has that order, but q is not above (n^(1/4) + 1)²
    expect(N.verifyCertificate({ n: 1000003n * 1000033n, kind: "elliptic", a: 1n, b: 0n, x: 0n, y: 0n, m: 2n, q: 2n, certificate: { n: 2n, kind: "small" } })).toMatchObject({ valid: false });
  });
  test("no certificate for a composite, even with points of large prime order on both factors", () => {
    // n = p·r; on a curve mod n take a point and m a multiple of its order mod p and mod r. Whatever prime q | m is
    // offered, the verifier must refuse: by size, or because the arithmetic meets a divisor that is not a unit.
    let offered = 0;
    for (const [p, r] of [[10007n, 10009n], [65537n, 65539n], [1009n, 1000003n], [100003n, 100019n]]) {
      const n = p * r;
      for (let a = 1n; a <= 6n; a++) {
        const b = a + 2n;
        let Cp, Cr;
        try {
          Cp = E([a, b], { mod: p });
          Cr = E([a, b], { mod: r });
        } catch {
          continue;
        }
        const m = Cp.count() * Cr.count();
        for (let x = 1n; x < 40n; x++) {
          const y2 = (x ** 3n + a * x + b) % n;
          const yp = N.sqrtMod(y2 % p, p);
          const yr = N.sqrtMod(y2 % r, r);
          if (yp === null || yr === null) continue;
          const y = N.crt([yp, yr], [p, r]).x ?? N.crt([yp, yr], [p, r]);
          for (const { prime: q } of N.factor(m)) {
            const cert = { n, kind: "elliptic", a, b, x, y: BigInt(y), m, q, certificate: N.provePrime(q).certificate };
            expect(N.verifyCertificate(cert).valid).toBe(false);
            offered++;
          }
          break;
        }
      }
    }
    expect(offered).toBeGreaterThan(40);
  });
});
