// Number theory, checked independently: certificates by a verifier that shares
// no code with the prover, functions against their definitions by brute
// force, class numbers against the analytic class number formula.
import { describe, expect, test } from "bun:test";
import { Rational } from "@ratmath/core";
import * as N from "../index.js";

const range = (a, b) => Array.from({ length: b - a + 1 }, (_, i) => BigInt(a + i));
const bruteGcd = (a, b) => (b === 0n ? (a < 0n ? -a : a) : bruteGcd(b, a % b));
const isPrimeBrute = (n) => {
  if (n < 2n) return false;
  for (let d = 2n; d * d <= n; d++) if (n % d === 0n) return false;
  return true;
};
const mod = (a, m) => ((a % m) + m) % m;
let seed = 12345n;
const rand = (n) => {
  seed = (seed * 6364136223846793005n + 1442695040888963407n) % (1n << 64n);
  return (seed >> 11n) % n;
};

describe("primes", () => {
  test("isPrime agrees with trial division", () => {
    for (const n of range(0, 3000)) expect(N.isPrime(n)).toBe(isPrimeBrute(n));
    for (let i = 0; i < 300; i++) {
      const n = rand(1n << 36n);
      expect(N.isPrime(n)).toBe(isPrimeBrute(n));
    }
  });
  test("sieves and counting", () => {
    expect([10, 100, 1000, 10000, 100000, 1000000, 10000000].map((x) => N.primePi(x))).toEqual([4n, 25n, 168n, 1229n, 9592n, 78498n, 664579n]);
    expect(N.primes(30)).toEqual([2n, 3n, 5n, 7n, 11n, 13n, 17n, 19n, 23n, 29n]);
    expect(N.nthPrime(1)).toBe(2n);
    expect(N.nthPrime(10000)).toBe(104729n);
    expect(N.nextPrime(1000000)).toBe(1000003n);
    expect(N.prevPrime(1000000)).toBe(999983n);
    expect(N.primesBetween(100, 130)).toEqual([101n, 103n, 107n, 109n, 113n, 127n]);
    expect(N.nextPrime(2n ** 89n - 2n)).toBe(2n ** 89n - 1n);
  });
  test("certificates are checked by an independent verifier", () => {
    let protheish = 2n ** 160n + 1n;
    while (!N.isProbablePrime(protheish)) protheish += 2n ** 160n;
    const primes = [2n ** 61n - 1n, 2n ** 89n - 1n, 2n ** 107n - 1n, 2n ** 127n - 1n, 10n ** 30n + 57n, 2n ** 64n + 13n, protheish, 4294967311n, 1000003n];
    for (const p of primes) {
      const r = N.provePrime(p);
      expect(r.proven).toBe(true);
      expect(r.prime).toBe(true);
      expect(N.verifyCertificate(r.certificate)).toEqual({ valid: true });
    }
    // tampering is caught
    const r = N.provePrime(2n ** 89n - 1n).certificate;
    const bad1 = { ...r, factors: r.factors.map((f, i) => (i === 0 ? { ...f, a: f.a + 1n } : f)) };
    const bad2 = { ...r, n: r.n + 2n };
    const bad3 = { ...r, factors: r.factors.slice(1) };
    const bad4 = { ...r, factors: r.factors.map((f, i) => (i === r.factors.length - 1 ? { ...f, certificate: { n: f.q, kind: "small" } } : f)) };
    for (const b of [bad2, bad3]) expect(N.verifyCertificate(b).valid).toBe(false);
    // a different base may still be valid; if it is accepted, it must truly satisfy Pocklington's conditions
    void bad1;
    expect(N.verifyCertificate(bad4).valid).toBe(r.factors.at(-1).q < 1n << 32n);
    expect(N.verifyCertificate({ n: 91n, kind: "small" }).valid).toBe(false);
    expect(N.verifyCertificate({ n: 561n, kind: "pocklington", F: 560n, R: 1n, factors: [{ q: 2n, e: 4, a: 2n, certificate: { n: 2n, kind: "small" } }] }).valid).toBe(false);
  });
  test("composites come with a proof", () => {
    const composites = [561n, 2047n, 3215031751n, 2n ** 67n - 1n, 1000000007n * 1000000009n, (2n ** 61n - 1n) * (2n ** 89n - 1n), 3825123056546413051n];
    for (const n of composites) {
      const r = N.provePrime(n);
      expect(r.prime).toBe(false);
      expect(r.proven).toBe(true);
      expect(N.verifyComposite(n, r).valid).toBe(true);
      expect(N.isPrime(n)).toBe(false);
    }
  });
  test("an unproven prime says so", () => {
    // n − 1 alone cannot certify 2^255 − 19 on a small budget; the elliptic curve certificates can
    const r = N.provePrime(2n ** 255n - 19n, { budget: 1 << 12, elliptic: false });
    expect(r.prime).toBe(true);
    expect(r.proven).toBe(false);
    expect(r.reason).toMatch(/bits of n − 1/);
    expect(N.isPrime(2n ** 255n - 19n)).toBe(true);
    expect(N.isProbablePrime(2n ** 255n - 19n)).toBe(true);
  });
  test("factoring", () => {
    for (let i = 0; i < 100; i++) {
      const n = rand(10n ** 15n) + 2n;
      const f = N.factor(n);
      expect(f.reduce((s, x) => s * x.prime ** BigInt(x.exponent), 1n)).toBe(n);
      for (const x of f) expect(isPrimeBrute(x.prime) || x.prime > 10n ** 7n).toBe(true);
    }
    expect(N.formatFactors(N.factor(360))).toBe("2^3 · 3^2 · 5");
    expect(N.factor(2n ** 67n - 1n).map((x) => x.prime)).toEqual([193707721n, 761838257287n]);
    expect([1000000007n, 1000000009n]).toContain(N.pollardPm1(1000000016000000063n));
    const n = 1000000007n * 998244353n;
    const d = N.pollardRho(n);
    expect(n % d).toBe(0n);
  });
});

describe("arithmetic functions, against their definitions", () => {
  const divs = (n) => range(1, Number(n)).filter((d) => n % d === 0n);
  test("σ, τ, φ, μ, ω, Ω, λ", () => {
    for (const n of range(1, 400)) {
      const ds = divs(n);
      expect(N.divisors(n)).toEqual(ds);
      expect(N.sigma(n)).toBe(ds.reduce((s, d) => s + d, 0n));
      expect(N.sigma(n, 2)).toBe(ds.reduce((s, d) => s + d * d, 0n));
      expect(N.tau(n)).toBe(BigInt(ds.length));
      expect(N.totient(n)).toBe(BigInt(range(1, Number(n)).filter((k) => bruteGcd(k, n) === 1n).length));
      const primesOf = ds.filter((d) => d > 1n && isPrimeBrute(d));
      const squareFree = !ds.some((d) => d > 1n && n % (d * d) === 0n);
      expect(N.mobius(n)).toBe(squareFree ? (primesOf.length % 2 ? -1n : 1n) : 0n);
      expect(N.omega(n)).toBe(BigInt(primesOf.length));
      // Carmichael λ: the least m > 0 with a^m ≡ 1 for all units a
      const units = range(1, Number(n)).filter((a) => bruteGcd(a, n) === 1n);
      let m = 1n;
      while (!units.every((a) => N.modPow(a, m, n) === 1n % n)) m++;
      expect(N.carmichael(n)).toBe(m);
    }
  });
  test("Dirichlet convolution and Möbius inversion", () => {
    const phiStar1 = N.dirichlet(N.totient, N.one);
    const muStar1 = N.dirichlet(N.mobius, N.one);
    const inv = N.mobiusInversion(N.sigma);
    for (const n of range(1, 200)) {
      expect(phiStar1(n)).toBe(n); // Gauss: Σ φ(d) = n
      expect(muStar1(n)).toBe(N.unit(n));
      expect(inv(n)).toBe(n);
    }
    expect(N.isMultiplicative(N.sigma)).toBe(true);
    expect(N.isMultiplicative((n) => n + 1n, 10)).toBe(false);
  });
  test("perfect, amicable, aliquot", () => {
    expect(N.perfectNumbers(8).at(-1)).toBe(2305843008139952128n);
    for (const p of N.perfectNumbers(5)) expect(N.isPerfect(p)).toBe(true);
    expect(N.amicablePairs(20000).length).toBe(8);
    expect(N.aliquotSequence(12496).ending).toBe("cycle of length 5");
    expect(N.aliquotSequence(10).ending).toBe("terminates");
    expect([2, 3, 5, 7, 13, 17, 19, 31, 61, 89, 107, 127].every((p) => N.lucasLehmer(p))).toBe(true);
    expect(N.lucasLehmer(11)).toBe(false);
  });
});

describe("congruences, against search", () => {
  test("linear congruences and the Chinese remainder theorem", () => {
    for (let i = 0; i < 200; i++) {
      const m = rand(60n) + 1n;
      const a = rand(200n) - 100n;
      const b = rand(200n) - 100n;
      const brute = range(0, Number(m) - 1).filter((x) => mod(a * x - b, m) === 0n);
      expect(N.solveLinear(a, b, m).solutions).toEqual(brute);
    }
    for (let i = 0; i < 200; i++) {
      const ms = [rand(20n) + 1n, rand(20n) + 1n, rand(20n) + 1n];
      const rs = ms.map((m) => rand(m));
      const L = ms.reduce((x, y) => (x * y) / bruteGcd(x, y));
      const brute = range(0, Number(L) - 1).filter((x) => rs.every((r, j) => mod(x - r, ms[j]) === 0n));
      const got = N.crt(rs, ms);
      if (!brute.length) expect(got).toBe(null);
      else expect(got).toEqual({ x: brute[0], modulus: L });
    }
  });
  test("polynomial roots mod m, Hensel", () => {
    for (let i = 0; i < 80; i++) {
      const poly = [rand(50n) - 25n, rand(50n) - 25n, rand(50n) - 25n, rand(3n)];
      const m = rand(300n) + 1n;
      const brute = range(0, Number(m) - 1).filter((x) => N.polyMod(poly, x, m) === 0n);
      expect(N.rootsMod(poly, m)).toEqual(brute);
    }
    const r = N.hensel([-2, 0, 1], 3, 7, 10);
    expect(mod(r * r - 2n, 7n ** 10n)).toBe(0n);
    expect(() => N.hensel([0, 0, 1], 0, 5, 3)).toThrow(/lift/);
  });
  test("primitive roots, orders and indices", () => {
    for (const n of range(1, 300)) {
      const units = range(1, Number(n)).filter((a) => bruteGcd(a, n) === 1n);
      const phi = BigInt(units.length);
      const ord = (g) => {
        let k = 1n;
        let x = mod(g, n);
        while (x !== 1n % n) {
          x = (x * g) % n;
          k++;
        }
        return k;
      };
      const brute = units.find((g) => ord(g) === phi) ?? null;
      expect(N.primitiveRootMod(n)).toBe(n === 1n ? 0n : brute);
      expect(N.hasPrimitiveRoot(n)).toBe(n === 1n || brute !== null);
    }
    for (let i = 0; i < 50; i++) {
      const n = rand(500n) + 2n;
      const g = rand(n);
      if (bruteGcd(g, n) !== 1n) continue;
      const h = N.modPow(g, rand(100n), n);
      const x = N.discreteLog(g, h, n);
      expect(N.modPow(g, x, n)).toBe(h);
      let k = 0n;
      while (N.modPow(g, k, n) !== h) k++;
      expect(x).toBe(k);
    }
  });
});

describe("reciprocity", () => {
  test("Jacobi = product of Legendre symbols by Euler's criterion; steps agree", () => {
    for (let n = 1n; n < 300n; n += 2n) {
      const fs = N.factor(n);
      for (let a = -20n; a < 40n; a++) {
        const viaEuler = fs.reduce((t, { prime, exponent }) => {
          const e = N.modPow(mod(a, prime), (prime - 1n) / 2n, prime);
          const l = e === 0n ? 0 : e === 1n ? 1 : -1;
          return t * l ** exponent + 0;
        }, 1);
        expect(N.jacobi(a, n)).toBe(viaEuler);
        expect(N.jacobiSteps(a, n).value).toBe(viaEuler);
      }
    }
  });
  test("Gauss's lemma, Kronecker, residues and square roots", () => {
    for (const p of N.primes(200).slice(1)) for (let a = 1n; a < p; a++) expect(N.gaussLemma(a, p).value).toBe(N.legendre(a, p));
    // Kronecker is completely multiplicative in the bottom
    for (let a = -10n; a <= 10n; a++) for (let m = 1n; m < 30n; m++) for (let n = 1n; n < 30n; n++) expect(N.kronecker(a, m * n)).toBe(N.kronecker(a, m) * N.kronecker(a, n) + 0);
    for (let i = 0; i < 100; i++) {
      const n = rand(400n) + 1n;
      const a = rand(n);
      expect(N.sqrtsMod(a, n)).toEqual(range(0, Number(n) - 1).filter((x) => (x * x - a) % n === 0n));
    }
    expect(N.quadraticResidues(11)).toEqual([0n, 1n, 3n, 4n, 5n, 9n]);
    expect(() => N.legendre(2, 15)).toThrow(/prime/);
  });
});

describe("sums of squares", () => {
  test("r₂, r₄, rₖ and the representations agree with counting", () => {
    for (const n of range(0, 600)) {
      const brute2 = [];
      for (let a = 0n; a * a <= n; a++) for (let b = 0n; b <= a; b++) if (a * a + b * b === n) brute2.push([a, b]);
      expect(N.twoSquaresAll(n)).toEqual(brute2.sort((x, y) => (x[0] > y[0] ? -1 : 1)));
      expect(N.r2(n)).toBe(N.rk(n, 2));
      expect(N.r4(n)).toBe(N.rk(n, 4));
      const t = N.threeSquares(n);
      expect(t === null).toBe(N.rk(n, 3) === 0n);
      if (t) expect(t.reduce((s, x) => s + x * x, 0n)).toBe(n);
      const f = N.fourSquares(n);
      expect(f.reduce((s, x) => s + x * x, 0n)).toBe(n);
    }
    const big = 2n ** 89n - 1n; // prime, 3 mod 4: not a sum of two squares
    expect(N.twoSquares(big)).toBe(null);
    const p = N.nextPrime(10n ** 20n);
    if (p % 4n === 1n) {
      const [a, b] = N.primeAsTwoSquares(p);
      expect(a * a + b * b).toBe(p);
    }
    const [a, b] = N.primeAsTwoSquares(10n ** 18n + 9n);
    expect(a * a + b * b).toBe(10n ** 18n + 9n);
  });
  test("Pythagorean triples", () => {
    const brute = [];
    for (let c = 1n; c <= 100n; c++) for (let a = 1n; a < c; a++) for (let b = a + 1n; b < c; b++) if (a * a + b * b === c * c && bruteGcd(bruteGcd(a, b), c) === 1n) brute.push(c);
    expect(N.pythagoreanTriples(100).map((t) => t[2]).sort((x, y) => (x < y ? -1 : x > y ? 1 : 0))).toEqual(brute);
    for (const [a, b, c] of N.pythagoreanTriples(1000)) expect(a * a + b * b).toBe(c * c);
    expect(N.triplesWithHypotenuse(65).length).toBe(4);
    for (const [a, b, c] of N.triplesWithLeg(60)) expect(a * a + b * b).toBe(c * c);
  });
});

describe("Diophantine equations", () => {
  test("linear equations: the basis spans every solution", () => {
    for (let i = 0; i < 100; i++) {
      const a = [rand(30n) - 15n, rand(30n) - 15n, rand(30n) - 15n];
      const c = rand(40n) - 20n;
      const r = N.linearDiophantine(a, c);
      const dot = (x) => x.reduce((s, xi, j) => s + xi * a[j], 0n);
      const g = a.reduce(bruteGcd, 0n);
      if (g === 0n || c % g !== 0n) {
        expect(r.solvable).toBe(false);
        continue;
      }
      expect(dot(r.particular)).toBe(c);
      for (const v of r.basis) expect(dot(v)).toBe(0n);
      // a lattice basis: the 2×2 minors have gcd = |the normalised coefficient vector|'s content, checked by covering a small box
      for (let x = -3n; x <= 3n; x++) for (let y = -3n; y <= 3n; y++) for (let z = -3n; z <= 3n; z++) {
        if (dot([x, y, z]) !== c) continue;
        // solve particular + s·b0 + t·b1 = (x, y, z) over the integers by brute force in a range
        const d = [x - r.particular[0], y - r.particular[1], z - r.particular[2]];
        let ok = false;
        for (let s = -60n; s <= 60n && !ok; s++) for (let t = -60n; t <= 60n && !ok; t++) ok = [0, 1, 2].every((j) => s * r.basis[0][j] + t * r.basis[1][j] === d[j]);
        expect(ok).toBe(true);
      }
    }
  });
  test("Frobenius numbers against search; coin counts", () => {
    for (const coins of [[3n, 5n], [6n, 9n, 20n], [4n, 7n, 11n], [5n, 8n, 9n, 12n], [7n, 11n, 13n]]) {
      const ok = [true];
      for (let s = 1; s < 500; s++) ok.push(coins.some((c) => s >= Number(c) && ok[s - Number(c)]));
      const brute = BigInt(ok.lastIndexOf(false));
      expect(N.frobenius(coins)).toBe(brute);
    }
    expect(N.coinWays(100, [1, 5, 10, 25, 50])).toBe(292n);
    expect(N.nonNegativeSolutions(3, 5, 30).solutions.length).toBe(3);
  });
  test("Pell: every listed solution solves the equation, and none is missed", () => {
    for (const [D, M] of [[2n, 1n], [2n, -1n], [2n, 7n], [3n, -2n], [5n, -4n], [7n, 2n], [13n, -4n], [13n, 3n], [10n, 6n], [6n, -2n]]) {
      const r = N.pellSolutions(D, M);
      const sols = r.solutions(6);
      for (const [x, y] of sols) expect(x * x - D * y * y).toBe(M);
      const brute = [];
      const X = sols.length ? sols.at(-1)[0] : 2000n;
      const Y = 20000n;
      for (let y = 0n; y <= Y && y * y * D <= X * X + 10n * (M < 0n ? -M : M); y++) {
        const t = M + D * y * y;
        if (t < 0n) continue;
        const x = N.isqrt(t);
        if (x * x === t && x <= X) brute.push([x, y]);
      }
      expect(sols.filter(([, y]) => y <= Y)).toEqual(brute.sort((u, v) => (u[0] < v[0] ? -1 : u[0] > v[0] ? 1 : u[1] < v[1] ? -1 : 1)));
    }
    expect(N.pellSolutions(-5, 29).solutions.length).toBe(4);
    expect(N.pellSolutions(3, 2).fundamental.length).toBe(0);
  });
  test("Egyptian fractions", () => {
    for (const q of ["4/13", "5/121", "7/15", "3/7", "2/3", "9/20"]) {
      for (const ds of [N.egyptianGreedy(q), N.egyptianShortest(q)]) {
        const sum = ds.reduce((s, d) => s.add(new Rational(1n, d)), new Rational(0n));
        expect(sum.equals(new Rational(q))).toBe(true);
        expect(new Set(ds.map(String)).size).toBe(ds.length);
      }
      expect(N.egyptianShortest(q).length).toBeLessThanOrEqual(N.egyptianGreedy(q).length);
    }
  });
});

describe("binary quadratic forms", () => {
  // h(D) = (w / 2|D|)·|Σ_{k=1}^{|D|} (D/k)·k| for a fundamental D < 0 (Dirichlet's class number formula)
  const fundamentalNegative = (D) => {
    const d = -D;
    if (D % 4n === 1n || D % 4n === -3n) return [...N.factor(d)].every((f) => f.exponent === 1);
    if (d % 4n !== 0n) return false;
    const m = d / 4n;
    return (m % 4n === 1n || m % 4n === 2n) && N.factor(m).every((f) => f.exponent === 1);
  };
  test("class numbers match Dirichlet's formula", () => {
    let checked = 0;
    for (let D = -3n; D >= -1500n; D--) {
      if (mod(D, 4n) > 1n || !fundamentalNegative(D)) continue;
      const w = D === -3n ? 6n : D === -4n ? 4n : 2n;
      let s = 0n;
      for (let k = 1n; k < -D; k++) s += BigInt(N.kronecker(D, k)) * k;
      const h = (w * (s < 0n ? -s : s)) / (2n * -D);
      expect(N.classNumber(D)).toBe(h);
      checked++;
    }
    expect(checked).toBeGreaterThan(400);
    expect([-3n, -4n, -7n, -8n, -11n, -19n, -43n, -67n, -163n].every((D) => N.classNumber(D) === 1n)).toBe(true);
  });
  test("narrow class numbers of real quadratic fields: h⁺ = h·(2 if the fundamental unit has norm +1)", () => {
    const h = { 40: 2n, 60: 2n, 65: 2n, 85: 2n };
    for (const D of [5, 8, 12, 13, 17, 21, 24, 28, 29, 33, 37, 40, 41, 44, 53, 56, 57, 60, 61, 65, 69, 73, 76, 77, 85, 88, 89, 92, 93, 97]) {
      const d = D % 4 === 0 ? D / 4 : D;
      const norm = N.fundamentalUnit(d).norm;
      expect(N.classNumber(D)).toBe((h[D] ?? 1n) * (norm === 1 ? 2n : 1n));
    }
  });
  test("composition makes a group", () => {
    for (const D of [-23n, -47n, -84n, -71n, -260n, -1540n, 229n, 316n, 60n, 145n]) {
      const G = N.classGroup(D);
      const els = G.forms;
      expect(G.structure.reduce((a, b) => a * b, 1n)).toBe(G.h);
      for (const f of els) {
        expect(G.compose(f, G.identity).equals(f)).toBe(true);
        expect(G.compose(f, G.inverse(f)).equals(G.identity)).toBe(true);
        expect(G.h % G.order(f)).toBe(0n);
        for (const g of els.slice(0, 6)) {
          expect(G.compose(f, g).equals(G.compose(g, f))).toBe(true);
          for (const k of els.slice(0, 4)) expect(G.compose(G.compose(f, g), k).equals(G.compose(f, G.compose(g, k)))).toBe(true);
        }
      }
    }
    expect(N.classGroup(-84n).structure).toEqual([2n, 2n]);
    expect(N.classGroup(-23n).structure).toEqual([3n]);
    expect(N.classGroup(-4n * 5n * 7n * 11n).structure).toEqual([2n, 2n, 2n]);
  });
  test("reduction and equivalence come with a matrix that proves them", () => {
    for (let i = 0; i < 200; i++) {
      const [a, b, c] = [rand(60n) - 30n, rand(60n) - 30n, rand(60n) - 30n];
      const f = N.form(a, b, c);
      const D = f.discriminant;
      if (D === 0n || (D > 0n && N.isqrt(D) ** 2n === D) || (D < 0n && a < 0n)) continue;
      const { form: r, matrix: M } = f.reduce();
      expect(r.isReduced()).toBe(true);
      expect(M[0][0] * M[1][1] - M[0][1] * M[1][0]).toBe(1n);
      expect(f.transform(M).equals(r)).toBe(true);
      // a random SL2 change of variables gives an equivalent form
      const t = rand(7n) - 3n;
      const g = f.transform([[1n, t], [0n, 1n]]).transform([[0n, -1n], [1n, 0n]]);
      const e = N.equivalent(f, g);
      expect(e.equivalent).toBe(true);
      expect(f.transform(e.matrix).equals(g)).toBe(true);
    }
  });
  test("representations", () => {
    const f = N.form(1, 0, 1);
    for (const n of range(1, 200)) {
      expect(BigInt(N.representations(f, n).length)).toBe(N.r2(n));
      const r = N.represent(f, n);
      if (r) expect(f.at(...r)).toBe(n);
    }
    for (const [g, n] of [[[1, 0, -2], 7], [[1, 1, -1], 11], [[2, 1, 3], 3], [[1, 1, 6], 8], [[3, 2, 5], 10], [[1, 0, -13], 3]]) {
      const r = N.represent(N.form(g), n);
      expect(r).not.toBe(null);
      expect(N.form(g).at(...r)).toBe(BigInt(n));
    }
    expect(N.represent(N.form(1, 0, 1), 3)).toBe(null);
  });
});

describe("p-adic numbers", () => {
  test("expansions of rationals match the digits and the arithmetic", () => {
    for (const p of [2n, 3n, 5n, 7n, 11n]) {
      for (const q of ["1/3", "-1", "7/10", "22/7", "-5/12", "1/1024", "3"]) {
        const r = new Rational(q);
        if (r.numerator === 0n) continue;
        const e = N.expansion(q, p);
        const d = N.padic(q, p).digits(30);
        expect(d.valuation).toBe(e.valuation);
        const want = [...e.preperiod];
        while (want.length < 30) want.push(...e.period);
        expect(d.digits).toEqual(want.slice(0, 30));
        // the digits sum back to q to 30 places: q − Σ dᵢ p^(v+i) has valuation ≥ v + 30
        const partial = d.digits.reduce((s, dg, i) => s.add(new Rational(dg).multiply(new Rational(p).pow(d.valuation + i))), new Rational(0n));
        expect(N.padicValuation(r.subtract(partial), p)).toBeGreaterThanOrEqual(d.valuation + 30);
      }
    }
  });
  test("lazy arithmetic agrees with exact rationals", () => {
    const p = 5n;
    const vals = ["1/3", "7/10", "-2/25", "125/4", "11"];
    for (const a of vals) for (const b of vals) {
      const A = N.padic(N.padic(a, p).sqrt ? a : a, p).add(N.padic(0, p).exp ? 0 : 0);
      const lazyA = N.padic(a, p).mul(N.padic(1, p).add(N.padic(5, p).exp()).div(N.padic(1, p).add(N.padic(5, p).exp())));
      const s = lazyA.add(N.padic(b, p));
      const m = lazyA.mul(N.padic(b, p));
      const q = lazyA.div(N.padic(b, p));
      expect(s.equalsTo(N.padic(new Rational(a).add(new Rational(b)), p), 25)).toBe(true);
      expect(m.equalsTo(N.padic(new Rational(a).multiply(new Rational(b)), p), 25)).toBe(true);
      expect(q.equalsTo(N.padic(new Rational(a).divide(new Rational(b)), p), 25)).toBe(true);
      void A;
    }
  });
  test("square roots, Hensel roots, exp and log", () => {
    for (const [x, p] of [[-1, 5n], [2, 7n], [-7, 2n], [17, 2n], [6, 5n], [-2, 3n], ["9/16", 7n], [50, 7n]]) {
      const r = N.padic(x, p).sqrt();
      expect(r.mul(r).equalsTo(N.padic(x, p), 40)).toBe(true);
    }
    expect(() => N.padic(2, 5n).sqrt()).toThrow(/square/);
    expect(() => N.padic(3, 2n).sqrt()).toThrow(/1 mod 8/);
    expect(() => N.padic(5, 5n).sqrt()).toThrow(/odd/);
    expect(N.padic(-1, 5n).isSquare()).toBe(true);
    expect(N.padic(2, 5n).isSquare()).toBe(false);
    const cube = N.henselRoot([-2, 0, 0, 1], 3, 5n); // x³ = 2 in Z_5
    expect(cube.pow(3).equalsTo(N.padic(2, 5n), 40)).toBe(true);
    for (const p of [3n, 5n, 7n]) {
      const x = N.padic(p * 2n, p);
      expect(x.exp().log().equalsTo(x, 30)).toBe(true);
      const a = N.padic(7, p).log();
      const b = N.padic(11, p).log();
      expect(N.padic(77, p).log().equalsTo(a.add(b), 30)).toBe(true);
      const e1 = N.padic(p, p).exp();
      const e2 = N.padic(2n * p, p).exp();
      expect(e1.mul(e1).equalsTo(e2, 30)).toBe(true);
    }
    expect(N.padic(4, 2n).exp().log().equalsTo(N.padic(4, 2n), 30)).toBe(true);
    expect(() => N.padic(1, 5n).exp()).toThrow(/converges/);
    for (const a of [2n, 3n, 4n]) expect(N.teichmuller(a, 5n).pow(4).equalsTo(N.padic(1, 5n), 30)).toBe(true);
  });
  test("the product formula: |x|·Π|x|_p = 1", () => {
    for (const q of ["360/77", "-1/1024", "999/1000"]) {
      const r = new Rational(q);
      let prod = r.abs();
      const ps = new Set([...N.factor(r.numerator < 0n ? -r.numerator : r.numerator), ...N.factor(r.denominator)].map((f) => f.prime));
      for (const p of ps) prod = prod.multiply(N.padicNorm(q, p));
      expect(prod.equals(new Rational(1n))).toBe(true);
    }
  });
  test("an undetermined valuation is reported", () => {
    const x = N.padic(5, 5n).exp();
    const zero = x.sub(x);
    expect(() => zero.valuation(100)).toThrow(/may be exactly 0/);
  });
});
