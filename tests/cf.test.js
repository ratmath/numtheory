import { describe, expect, test } from "bun:test";
// Moved here from @ratmath/algebraic with the continued fractions themselves.
import { continuedFraction, convergents, fromContinuedFraction, periodicContinuedFraction, pell, fundamentalUnit, fromQuadratic } from "../index.js";
import { sqrt, nthRoot, fromRational, format, roots } from "@ratmath/algebraic";

const q = (x) => fromRational(x);
const nums = (a) => a.map(String);

describe("continued fractions", () => {
  test("rationals and algebraic numbers", () => {
    expect(nums(continuedFraction("355/113"))).toEqual(["3", "7", "16"]);
    expect(nums(continuedFraction(sqrt(q(2)), 5))).toEqual(["1", "2", "2", "2", "2"]);
    expect(nums(continuedFraction(nthRoot(q(2), 3), 12))).toEqual(["1", "3", "1", "5", "1", "1", "4", "1", "1", "8", "1", "14"]);
    expect(nums(continuedFraction(fromQuadratic(0, -1, 3, 1), 4))).toEqual(["-2", "3", "1", "2"]);
    expect(convergents([1, 2, 2, 2]).map(String)).toEqual(["1", "3/2", "7/5", "17/12"]);
  });
  test("quadratic irrationals are periodic", () => {
    const p = periodicContinuedFraction(sqrt(q(7)));
    expect(nums(p.preperiod)).toEqual(["2"]);
    expect(nums(p.period)).toEqual(["1", "1", "1", "4"]);
    expect(nums(periodicContinuedFraction(roots([-1, -1, 1])[1].value).period)).toEqual(["1"]);
    expect(format(fromContinuedFraction([1], [2]))).toBe("√2");
    const x = fromQuadratic(3, 1, 5, 7);
    const { preperiod, period } = periodicContinuedFraction(x);
    expect(format(fromContinuedFraction(preperiod, period))).toBe("(3 + √5)/7");
  });
});

describe("Pell equations and units", () => {
  test("pell", () => {
    expect(pell(61)).toEqual({ x: 1766319049n, y: 226153980n, negative: { x: 29718n, y: 3805n } });
    expect(pell(94).negative).toBe(null);
    const { x, y } = pell(94);
    expect(x * x - 94n * y * y).toBe(1n);
  });
  test("fundamental units", () => {
    expect(String(fundamentalUnit(2).unit)).toBe("1 + √2");
    expect(String(fundamentalUnit(5).unit)).toBe("1/2 + (1/2)√5");
    expect(fundamentalUnit(5).norm).toBe(-1);
    expect(String(fundamentalUnit(94).unit)).toBe("2143295 + 221064√94");
    expect(fundamentalUnit(21).norm).toBe(1);
  });
});

