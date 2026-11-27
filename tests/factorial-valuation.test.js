import { describe, expect, test } from "bun:test";
import { factorialValuation, valuation } from "../index.js";

const fact = (n) => { let f = 1n; for (let k = 2n; k <= BigInt(n); k++) f *= k; return f; };
const digitSum = (n, p) => { let s = 0n; for (let m = BigInt(n); m > 0n; m /= BigInt(p)) s += m % BigInt(p); return s; };

describe("factorialValuation", () => {
  test("Legendre: the largest power of 2 in 100! is 97", () => {
    expect(factorialValuation(100, 2)).toBe(97);
    expect(factorialValuation(100, 5)).toBe(24);
    expect(factorialValuation(0, 2)).toBe(0);
    expect(factorialValuation(1, 3)).toBe(0);
  });

  test("agrees with the valuation of n! itself (independent check)", () => {
    for (let n = 0; n <= 120; n += 7) {
      const f = fact(n);
      for (const b of [2, 3, 4, 6, 7, 10, 12, 16, 36]) {
        // the largest k with b^k | n!
        let k = 0, r = f;
        while (r % BigInt(b) === 0n) { r /= BigInt(b); k++; }
        expect(factorialValuation(n, b)).toBe(k);
      }
      for (const p of [2, 3, 5, 7, 11]) expect(factorialValuation(n, p)).toBe(valuation(f, p));
    }
  });

  test("agrees with (n − s_p(n))/(p − 1) for huge n", () => {
    const n = 10n ** 30n + 12345n;
    for (const p of [2n, 3n, 5n, 97n]) expect(BigInt(factorialValuation(n, p))).toBe((n - digitSum(n, p)) / (p - 1n));
  });

  test("refuses bad input", () => {
    expect(() => factorialValuation(-1, 2)).toThrow();
    expect(() => factorialValuation(5, 1)).toThrow();
    expect(() => factorialValuation(2.5, 2)).toThrow();
  });
});
