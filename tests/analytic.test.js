import { test, expect, describe } from "bun:test";
import * as N from "../index.js";

const S = (x) => String(x);
const num = (q) => Number(q.numerator) / Number(q.denominator);
const mid = (iv) => (num(iv.low) + num(iv.high)) / 2;
const width = (iv) => num(iv.high.subtract(iv.low));

describe("counting primes", () => {
  test("π(x) by Meissel's formula", () => {
    expect(N.primePi(10n ** 8n)).toBe(5761455n);
    expect(N.primePi(10n ** 9n)).toBe(50847534n);
    expect(N.primePi(10n ** 10n)).toBe(455052511n);
    expect(N.primePi(10n ** 11n)).toBe(4118054813n);
    expect(N.primePi(2n ** 32n)).toBe(203280221n);
    expect(() => N.primePi(10n ** 13n)).toThrow();
  });
  test("the formula agrees with the sieve on small and awkward x", () => {
    for (const x of [1, 2, 3, 4, 7, 8, 9, 26, 27, 28, 63, 64, 65, 100, 124, 125, 126, 1000, 12345, 99991, 1000003]) {
      expect(N.meissel(x).pi).toBe(N.primePi(x));
    }
    const m = N.meissel(1000);
    expect([m.y, m.a, m.phi, m.P2]).toEqual([10n, 4n, 228n, 63n]);
    expect(m.phi + m.a - 1n - m.P2).toBe(168n);
  });
  test("Legendre's φ", () => {
    expect(N.legendrePhi(100, 3)).toBe(26n);
    expect(N.legendrePhi(1000, 4)).toBe(228n);
    expect(N.legendrePhi(30, 3)).toBe(8n);
    expect(N.legendrePhi(50, 0)).toBe(50n);
    // π(x) = φ(x, a) + a − 1 with a = π(√x)
    expect(N.legendrePhi(10000, 25) + 25n - 1n).toBe(1229n);
  });
});

describe("Dirichlet characters", () => {
  test("the group mod 15", () => {
    const G = N.dirichletGroup(15);
    expect(G.order).toBe(8);
    expect(G.structure).toEqual([2, 4]);
    expect(G.primitive.map((c) => c.number)).toEqual([2, 8, 14]);
    expect(G.real.map((c) => c.number)).toEqual([1, 4, 11, 14]);
    expect(G.orthogonality(2, 2)).toEqual({ re: 8, im: 0 });
    expect(G.orthogonality(2, 7)).toEqual({ re: 0, im: 0 });
  });
  test("Conrey's numbering", () => {
    const c = N.dirichletCharacter(7, 3);
    expect(c.order).toBe(6);
    expect([1, 2, 3, 4, 5, 6].map((n) => S(c.exponent(n)))).toEqual(["0", "1/3", "1/6", "2/3", "5/6", "1/2"]);
    expect(c.exponent(7)).toBe(null);
    expect(c.isEven).toBe(false);
    // symmetric: χ_q(m, n) = χ_q(n, m)
    for (const q of [9, 16, 20, 35]) {
      const units = N.dirichletGroup(q).characters.map((c) => c.number);
      for (const m of units) for (const n of units) expect(S(N.dirichletCharacter(q, m).exponent(n))).toBe(S(N.dirichletCharacter(q, n).exponent(m)));
    }
    expect(() => N.dirichletCharacter(15, 3)).toThrow();
  });
  test("multiplicativity, products and conjugates", () => {
    const c = N.dirichletCharacter(20, 3);
    for (const [a, b] of [[3, 7], [9, 11], [13, 17]]) {
      const e = c.exponent(a).add(c.exponent(b));
      const want = c.exponent(a * b);
      expect(S(e.subtract(want).denominator)).toBe("1");
    }
    expect(S(N.dirichletCharacter(3, 2).mul(N.dirichletCharacter(5, 2)))).toBe("χ_15(2, ·)");
    expect(c.mul(c.conjugate()).isPrincipal).toBe(true);
    expect(c.pow(c.order).isPrincipal).toBe(true);
  });
  test("conductors and inducing characters", () => {
    expect(N.dirichletCharacter(15, 11).conductor).toBe(3);
    expect(S(N.dirichletCharacter(15, 11).primitive())).toBe("χ_3(2, ·)");
    expect(N.dirichletCharacter(16, 7).conductor).toBe(8);
    expect(N.dirichletCharacter(16, 15).conductor).toBe(4);
    expect(N.dirichletCharacter(16, 3).conductor).toBe(16);
    for (const q of [8, 12, 16, 45]) {
      for (const c of N.dirichletGroup(q).characters) {
        const p = c.primitive();
        expect(p.isPrimitive).toBe(true);
        for (let n = 1; n < q; n++) if (c.exponent(n) !== null) expect(S(p.exponent(n))).toBe(S(c.exponent(n)));
      }
    }
  });
  test("Gauss sums of primitive characters have absolute value √q", () => {
    for (const c of N.dirichletGroup(13).primitive) {
      const g = c.gaussSum();
      expect(Math.abs(g.re * g.re + g.im * g.im - 13)).toBeLessThan(1e-9);
    }
  });
  test("the Kronecker symbol as a character", () => {
    const k = N.kroneckerCharacter(-23);
    expect(k.modulus).toBe(23);
    expect(k.isPrimitive).toBe(true);
    expect(k.isEven).toBe(false);
    for (let n = 0; n < 23; n++) expect(k.at(n)).toBe(BigInt(N.kronecker(-23, n)));
    expect(N.kroneckerCharacter(5).isEven).toBe(true);
    expect(N.kroneckerCharacter(12).conductor).toBe(12);
  });
});

describe("exact values", () => {
  test("Bernoulli numbers and polynomials", () => {
    expect([0, 1, 2, 4, 12].map((n) => S(N.bernoulli(n)))).toEqual(["1", "-1/2", "1/6", "-1/30", "-691/2730"]);
    expect(S(N.bernoulli(7))).toBe("0");
    expect(S(N.bernoulliPolynomial(2, "1/2"))).toBe("-1/12");
    expect(S(N.bernoulliPolynomial(3, 1))).toBe("0");
  });
  test("ζ at integers", () => {
    expect(S(N.zetaExact(0))).toBe("-1/2");
    expect(S(N.zetaExact(-1))).toBe("-1/12");
    expect(S(N.zetaExact(-3))).toBe("1/120");
    expect(S(N.zetaExact(-2))).toBe("0");
    expect(S(N.zetaExact(2).coefficient)).toBe("1/6");
    expect(S(N.zetaExact(4).coefficient)).toBe("1/90");
    expect(S(N.zetaExact(12).coefficient)).toBe("691/638512875");
    expect(() => N.zetaExact(1)).toThrow();
    expect(() => N.zetaExact(3)).toThrow();
  });
  test("L(0, χ_D) is the class number for D < −4", () => {
    for (const [D, h] of [[-7, 1], [-23, 3], [-47, 5], [-163, 1], [-20, 2], [-84, 4]]) {
      expect(S(N.lValueExact(N.kroneckerCharacter(D), 0))).toBe(String(h));
      expect(N.classNumber(D)).toBe(BigInt(h));
    }
    expect(S(N.lValueExact(N.kroneckerCharacter(5), -1))).toBe("-2/5");
    expect(S(N.lValueExact(N.kroneckerCharacter(-4), 0))).toBe("1/2");
    expect(S(N.lValueExact(N.dirichletCharacter(1, 1), -1))).toBe("-1/12");
  });
  test("a complex character: the value lives in a cyclotomic field", () => {
    const c = N.dirichletCharacter(5, 2); // order 4
    const v = N.lValueExact(c, 0);
    // L(0, χ) = −(1/5)Σ a·χ(a) = −(1 + 2i − 3i − 4)/5 = (3 + i)/5 with i = ζ4
    expect(S(v)).toBe("3/5 + (1/5)ζ4");
    const enc = N.lValue(c, 0);
    expect(Math.abs(mid(enc.re) - 0.6)).toBeLessThan(1e-15);
    expect(Math.abs(mid(enc.im) - 0.2)).toBeLessThan(1e-15);
  });
});

describe("enclosures", () => {
  test("ζ on the real line", () => {
    const z2 = N.zetaValue(2);
    expect(width(z2)).toBeLessThan(2 ** -63);
    expect(Math.abs(mid(z2) - Math.PI ** 2 / 6)).toBeLessThan(1e-15);
    expect(Math.abs(mid(N.zetaValue(-3)) - 1 / 120)).toBeLessThan(1e-16);
    expect(Math.abs(mid(N.zetaValue("1/2")) + 1.4603545088095868)).toBeLessThan(1e-14);
    expect(Math.abs(mid(N.zetaValue(0)) + 0.5)).toBeLessThan(1e-16);
    expect(() => N.zetaValue(1)).toThrow();
    // the exact value lies inside the enclosure
    const z = N.zetaValue(-5, { bits: 100 });
    const exact = N.zetaExact(-5);
    expect(z.low.lessThanOrEqual(exact) && exact.lessThanOrEqual(z.high)).toBe(true);
    expect(width(N.zetaValue(3, { bits: 200 }))).toBeLessThan(2 ** -199);
  });
  test("ζ in the plane, near its first zero", () => {
    const z = N.zetaValue(["1/2", 14]);
    expect(Math.abs(mid(z.re) - 0.022241142609993589)).toBeLessThan(1e-14);
    expect(Math.abs(mid(z.im) + 0.10325812326645005)).toBeLessThan(1e-14);
    const far = N.zetaValue([2, 100], { bits: 40 });
    expect(Math.abs(mid(far.re) - 1.1899)).toBeLessThan(0.3);
  });
  test("L(1, χ): Leibniz and the class number formula", () => {
    expect(Math.abs(mid(N.lValue(N.kroneckerCharacter(-4), 1)) - Math.PI / 4)).toBeLessThan(1e-15);
    expect(Math.abs(mid(N.lValue(N.kroneckerCharacter(-23), 1)) - (3 * Math.PI) / Math.sqrt(23))).toBeLessThan(1e-14);
    // real quadratic: L(1, χ_8) = ln(1 + √2)/√2
    expect(Math.abs(mid(N.lValue(N.kroneckerCharacter(8), 1)) - Math.log(1 + Math.SQRT2) / Math.SQRT2)).toBeLessThan(1e-14);
    expect(() => N.lValue(N.dirichletCharacter(5, 1), 1)).toThrow();
  });
  test("L(s, χ) agrees with exact values and with an Euler product", () => {
    const k = N.kroneckerCharacter(-23);
    const l0 = N.lValue(k, 0);
    expect(l0.low.lessThanOrEqual(N.lValueExact(k, 0)) && N.lValueExact(k, 0).lessThanOrEqual(l0.high)).toBe(true);
    const c = N.dirichletCharacter(7, 3);
    const v = N.lValue(c, [3, 2], { bits: 40 });
    // partial Dirichlet series
    let re = 0;
    let im = 0;
    for (let n = 1; n < 20000; n++) {
      const x = c.value(n);
      const mag = n ** -3;
      const ang = -2 * Math.log(n);
      re += mag * (x.re * Math.cos(ang) - x.im * Math.sin(ang));
      im += mag * (x.re * Math.sin(ang) + x.im * Math.cos(ang));
    }
    expect(Math.abs(mid(v.re) - re)).toBeLessThan(1e-7);
    expect(Math.abs(mid(v.im) - im)).toBeLessThan(1e-7);
  });
  test("Hurwitz zeta", () => {
    // ζ(2, 1/4) = π² + 8G
    expect(Math.abs(mid(N.hurwitzZeta(2, "1/4")) - (Math.PI ** 2 + 8 * 0.915965594177219))).toBeLessThan(1e-12);
    expect(Math.abs(mid(N.hurwitzZeta(2, 1)) - Math.PI ** 2 / 6)).toBeLessThan(1e-15);
    expect(() => N.hurwitzZeta(2, "3/2")).toThrow();
  });
});
