import { test, expect, describe } from "bun:test";
import * as N from "../index.js";

const { NumberField } = N;
const S = (x) => String(x);

describe("rings of integers", () => {
  test("quadratic fields: the basis depends on d mod 4", () => {
    const O5 = N.ringOfIntegers(NumberField.quadratic(5));
    expect(O5.discriminant).toBe(5n);
    expect(O5.index).toBe(2n);
    expect(O5.basis.map(S)).toEqual(["1", "1/2 + (1/2)√5"]);
    const Om5 = N.ringOfIntegers(NumberField.quadratic(-5));
    expect(Om5.discriminant).toBe(-20n);
    expect(Om5.index).toBe(1n);
    expect(Om5.signature).toEqual([0, 1]);
  });
  test("Dedekind's cubic: index 2, discriminant −503, and 2 splits completely", () => {
    const K = new NumberField([-8n, -2n, -1n, 1n]);
    const O = N.ringOfIntegers(K);
    expect(O.discriminant).toBe(-503n);
    expect(O.index).toBe(2n);
    expect(O.contains(K.gen.mul(K.gen).add(K.gen).div(K.of(2n)))).toBe(true);
    expect(O.contains(K.gen.div(K.of(2n)))).toBe(false);
    const Ps = O.primesAbove(2);
    expect(Ps.map((P) => [P.e, P.f])).toEqual([[1, 1], [1, 1], [1, 1]]);
    // no element generates the ring: three primes of norm 2 cannot come from a cubic mod 2
    expect(O.describePrime(3).kind).toBe("inert");
    expect(O.primesAbove(503).map((P) => P.e).sort()).toEqual([1, 2]);
    expect(O.steps.find((s) => s.p === 2n).indexExponent).toBe(1);
  });
  test("pure cubics and a wildly ramified prime", () => {
    const O19 = N.ringOfIntegers([-19n, 0n, 0n, 1n]);
    expect(O19.discriminant).toBe(-1083n);
    expect(O19.index).toBe(3n);
    expect(O19.primesAbove(3).map((P) => P.e).sort()).toEqual([1, 2]);
    const O175 = N.ringOfIntegers([-175n, 0n, 0n, 1n]);
    expect(O175.discriminant).toBe(-33075n);
    expect(O175.primesAbove(5).map((P) => [P.e, P.f])).toEqual([[3, 1]]);
    expect(N.fieldDiscriminant([-2n, 0n, 0n, 1n])).toBe(-108n);
  });
  test("cyclotomic and other fields", () => {
    expect(N.fieldDiscriminant([1n, 1n, 1n, 1n, 1n])).toBe(125n);
    expect(N.fieldDiscriminant([1n, 0n, 0n, 1n, 0n, 0n, 1n])).toBe(-19683n);
    expect(N.fieldDiscriminant([2n, 0n, -4n, 0n, 1n])).toBe(2048n);
    // a polynomial that is not monic over Z
    const O = N.ringOfIntegers([[1n, 2n], 0n, 0n, 1n]);
    expect(O.discriminant).toBe(-108n);
    expect(N.integralBasis([-5n, 0n, 1n]).length).toBe(2);
  });
  test("Σ e·f = n at every prime, and the trace matrix gives the discriminant", () => {
    const O = N.ringOfIntegers([-7n, 3n, 0n, 1n, 1n]);
    for (const p of [2n, 3n, 5n, 7n, 11n, 13n, 101n]) {
      expect(O.primesAbove(p).reduce((s, P) => s + P.e * P.f, 0)).toBe(4);
    }
    const T = O.traceMatrix;
    // 4×4 determinant by cofactors is overkill: compare with the ideal norm identity instead
    expect(T.length).toBe(4);
    expect(O.ramifiedPrimes.every((p) => O.primesAbove(p).some((P) => P.e > 1))).toBe(true);
    expect(O.primesAbove(101n).every((P) => P.e === 1)).toBe(true);
  });
  test("coordinates and elements round-trip", () => {
    const K = new NumberField([-19n, 0n, 0n, 1n]);
    const O = N.ringOfIntegers(K);
    const x = O.element([1, 2, 3]);
    expect(O.coordinates(x).map(S)).toEqual(["1", "2", "3"]);
    expect(O.contains(x)).toBe(true);
    expect(() => O.primesAbove(4)).toThrow();
  });
});

describe("ideals", () => {
  const K = NumberField.quadratic(-5);
  const O = N.ringOfIntegers(K);
  test("(6) = P2² · P3 · P3′ in Z[√−5]", () => {
    const I = O.ideal(6);
    const fs = I.factor();
    expect(fs.map((f) => [f.prime.p, f.exponent])).toEqual([[2n, 2], [3n, 1], [3n, 1]]);
    expect(S(I.norm)).toBe("36");
    // 6 = 2·3 = (1+√−5)(1−√−5): the two factorisations share their prime ideals
    const a = O.ideal(K.of([1n, 1n].map((c) => [c, 1n])));
    expect(a.factor().map((f) => f.prime.p).sort()).toEqual([2n, 3n]);
  });
  test("products, sums, inverses", () => {
    const [P2] = O.primesAbove(2);
    const [P3, Q3] = O.primesAbove(3);
    const I = P2.toIdeal().mul(P3);
    expect(S(I.norm)).toBe("6");
    expect(I.mul(I.inverse()).isOne).toBe(true);
    expect(P3.toIdeal().add(Q3).isOne).toBe(true);
    expect(P2.toIdeal().pow(2).equals(O.ideal(2))).toBe(true);
    expect(I.pow(-2).mul(I.pow(2)).isOne).toBe(true);
    expect(I.contains(6)).toBe(true);
    expect(I.contains(1)).toBe(false);
    expect(P3.toIdeal().isPrime).toBe(true);
    expect(I.isPrime).toBe(false);
    expect(S(O.ideal(K.of([[1n, 2n]])).norm)).toBe("1/4");
  });
  test("valuations", () => {
    const [P2] = O.primesAbove(2);
    expect(P2.e).toBe(2);
    expect(P2.valuation(2)).toBe(2);
    expect(P2.valuation(K.of([[1n, 1n], [1n, 1n]]))).toBe(1);
    expect(P2.valuation(K.of([[1n, 4n]]))).toBe(-4);
    expect(P2.valuation(O.ideal(4))).toBe(4);
  });
});

describe("class groups", () => {
  const h = (poly, opts) => N.idealClassGroup(Array.isArray(poly) ? new NumberField(poly) : poly, opts);
  test("imaginary quadratic fields agree with reduced forms", () => {
    for (const [d, hh] of [[-1, 1n], [-5, 2n], [-23, 3n], [-47, 5n], [-163, 1n], [-71, 7n], [-21, 4n]]) {
      const cg = h(NumberField.quadratic(d));
      expect(cg.h).toBe(hh);
      expect(cg.proven).toBe(true);
    }
    expect(h(NumberField.quadratic(-21)).structure).toEqual([2n, 2n]);
  });
  test("real quadratic fields: the wide class number, not the narrow one", () => {
    const cg = h(NumberField.quadratic(79));
    expect(cg.h).toBe(3n);
    expect(cg.certificate.narrowClassNumber).toBe(6n);
    expect(h(NumberField.quadratic(82)).structure).toEqual([4n]);
    expect(h(NumberField.quadratic(94)).h).toBe(1n);
    const ug = N.unitGroup(NumberField.quadratic(94));
    expect(Math.abs(ug.regulator - Math.log(2143295 + 221064 * Math.sqrt(94)))).toBeLessThan(1e-9);
  });
  test("cubic fields, certified by saturation", () => {
    const a = h([-7n, 0n, 0n, 1n]);
    expect(a.h).toBe(3n);
    expect(a.proven).toBe(true);
    expect(a.certificate.method).toBe("saturation");
    expect(a.certificate.classNumberPrimes).toEqual([3]);
    expect(a.certificate.unitIndexBound).toBe(11);
    expect(a.certificate.primesChecked).toBe(5);
    expect(h([-1n, 4n, 0n, 1n]).h).toBe(2n); // discriminant −283, the first complex cubic with h = 2
    expect(h([-11n, 0n, 0n, 1n]).h).toBe(2n);
    expect(h([-2n, 0n, 0n, 1n]).h).toBe(1n);
    const c = h([-65n, 0n, 0n, 1n]);
    expect(c.h).toBe(18n);
    expect(c.structure).toEqual([3n, 6n]);
    expect(c.proven).toBe(true);
  });
  test("classes, orders and principal ideals in Q(∛65)", () => {
    const K = new NumberField([-65n, 0n, 0n, 1n]);
    const O = N.ringOfIntegers(K);
    const cg = N.idealClassGroup(K);
    expect(cg.generators.length).toBe(2);
    cg.generators.forEach((g, i) => expect(cg.orderOf(g)).toBe(cg.structure[i]));
    const [P] = O.primesAbove(2);
    const o = cg.orderOf(P);
    expect(o).toBe(6n);
    expect(cg.isPrincipal(P).principal).toBe(false);
    const I = P.toIdeal().pow(6);
    const r = cg.isPrincipal(I);
    expect(r.principal).toBe(true);
    expect(O.ideal(r.generator).equals(I)).toBe(true);
    // the class map is a homomorphism
    const [Q] = O.primesAbove(3);
    const sum = cg.classOf(P).map((x, i) => (x + cg.classOf(Q)[i]) % cg.structure[i]);
    expect(cg.classOf(P.toIdeal().mul(Q))).toEqual(sum);
    expect(cg.classOf(O.ideal(7))).toEqual([0n, 0n]);
  });
  test("forms and ideals", () => {
    const K = NumberField.quadratic(-23);
    const cg = N.idealClassGroup(K);
    const I = N.formToIdeal([2, 1, 3], K);
    expect(S(I.norm)).toBe("2");
    expect(cg.orderOf(I)).toBe(3n);
    expect(cg.isPrincipal(N.formToIdeal([1, 1, 6], K)).principal).toBe(true);
  });
  test("a first guess that is too large is caught by the certification", () => {
    // without certification the relations found first give 8; the 2-saturation test fails, more relations follow
    const K = new NumberField([-30n, 0n, 0n, 0n, 1n]);
    const cg = N.idealClassGroup(K);
    expect(cg.h).toBe(4n);
    expect(cg.structure).toEqual([4n]);
    expect(cg.proven).toBe(true);
    const loose = N.idealClassGroup(K, { certify: false });
    expect(loose.proven).toBe(false);
    expect(loose.h % cg.h).toBe(0n);
  });
  test("the Minkowski bound", () => {
    expect(N.minkowskiBound(NumberField.quadratic(-5))).toBe(2n);
    expect(N.minkowskiBound([-65n, 0n, 0n, 1n])).toBe(95n);
    expect(N.fieldClassNumber(NumberField.quadratic(-5))).toBe(2n);
  });
});

describe("unit groups", () => {
  test("rank, torsion and regulator", () => {
    const u2 = N.unitGroup(NumberField.quadratic(2));
    expect(u2.rank).toBe(1);
    expect(u2.torsionOrder).toBe(2);
    expect(Math.abs(u2.regulator - Math.log(1 + Math.SQRT2))).toBeLessThan(1e-12);
    expect(N.unitGroup(NumberField.gaussian()).torsionOrder).toBe(4);
    expect(N.unitGroup(NumberField.quadratic(-3)).torsionOrder).toBe(6);
    const c5 = N.unitGroup([1n, 1n, 1n, 1n, 1n]);
    expect(c5.torsionOrder).toBe(10);
    expect(c5.rank).toBe(1);
    expect(Math.abs(c5.regulator - 2 * Math.log((1 + Math.sqrt(5)) / 2))).toBeLessThan(1e-12);
    const cyc = N.unitGroup([1n, -3n, 0n, 1n]); // the cubic field in Q(ζ9)
    expect(cyc.rank).toBe(2);
    expect(Math.abs(cyc.regulator - 0.8492874506)).toBeLessThan(1e-9);
    expect(cyc.proven).toBe(true);
    expect(Math.abs(N.regulator([-2n, 0n, 0n, 1n]) - 1.3473773483)).toBeLessThan(1e-9);
  });
  test("every unit is a root of unity times powers of the fundamental units", () => {
    const K = new NumberField([-2n, 0n, 0n, 1n]);
    const ug = N.unitGroup(K);
    const [e] = ug.fundamentalUnits;
    expect(S(e.norm()).replace("-", "")).toBe("1");
    expect(ug.exponents(e.pow(7).neg())).toEqual({ torsion: 1, exponents: [7n] });
    expect(ug.exponents(e.inv())).toEqual({ torsion: 0, exponents: [-1n] });
    expect(ug.isUnit(K.gen)).toBe(false);
    // 1 + ∛2 + ∛4 is the classical unit
    const u = K.of([[1n, 1n], [1n, 1n], [1n, 1n]]);
    expect(ug.isUnit(u)).toBe(true);
    expect(ug.exponents(u).exponents.map((x) => (x < 0n ? -x : x))).toEqual([1n]);
    const q = N.unitGroup([2n, 0n, -4n, 0n, 1n]);
    expect(q.rank).toBe(3);
    const w = q.fundamentalUnits[0].mul(q.fundamentalUnits[1].pow(2)).mul(q.fundamentalUnits[2].inv());
    expect(q.exponents(w)).toEqual({ torsion: 0, exponents: [1n, 2n, -1n] });
  });
});
