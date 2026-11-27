// Small integer helpers shared by the package: reading whole numbers,
// gcd, integer roots, and the error type.

import { Rational } from "@ratmath/core";

export class NumberTheoryError extends Error {
  constructor(code, message) {
    super(message);
    this.name = "NumberTheoryError";
    this.code = code;
  }
}
export const refuse = (code, message) => {
  throw new NumberTheoryError(code, message);
};

/** A whole number as a BigInt: a JavaScript integer, a BigInt, decimal text, or a Rational with denominator 1. */
export function big(x, what = "argument") {
  if (typeof x === "bigint") return x;
  if (typeof x === "number") {
    if (!Number.isInteger(x)) throw new TypeError(`${what} must be a whole number, not ${x}`);
    if (!Number.isSafeInteger(x)) throw new RangeError(`${what} = ${x} is beyond 2^53 as a JavaScript number; pass a BigInt`);
    return BigInt(x);
  }
  if (typeof x === "string" && /^\s*[-+]?\d+\s*$/.test(x)) return BigInt(x.trim());
  if (x instanceof Rational || (x && typeof x.numerator === "bigint" && typeof x.denominator === "bigint")) {
    if (x.denominator === 1n) return x.numerator;
  }
  if (x && typeof x === "object" && typeof x.value === "bigint") return x.value; // Core's Integer
  throw new TypeError(`${what} must be a whole number, not ${typeof x === "string" ? JSON.stringify(x) : String(x)}`);
}
/** A whole number ≥ lo. */
export function atLeast(x, lo, what = "argument") {
  const n = big(x, what);
  if (n < lo) throw new RangeError(`${what} must be at least ${lo}, not ${n}`);
  return n;
}
/** An exact Rational from a whole number, a fraction text ("3/7", "0.25"), a BigInt or a Rational. */
export function rat(x, what = "argument") {
  if (x instanceof Rational) return x;
  if (typeof x === "bigint") return new Rational(x);
  if (typeof x === "number") return new Rational(big(x, what));
  if (typeof x === "string" && x.trim() !== "") return new Rational(x.trim());
  if (x && typeof x.numerator === "bigint") return new Rational(x.numerator, x.denominator);
  throw new TypeError(`${what} must be a whole number or a fraction, not ${String(x)}`);
}

export const abs = (a) => (a < 0n ? -a : a);
export const mod = (a, m) => ((a % m) + m) % m;
export function gcd(a, b) {
  a = abs(a);
  b = abs(b);
  while (b) [a, b] = [b, a % b];
  return a;
}
export const lcm = (a, b) => (a === 0n || b === 0n ? 0n : abs((a / gcd(a, b)) * b));
/** Floor division, rounding towards −∞. */
export function fdiv(a, b) {
  const q = a / b;
  return a % b !== 0n && (a < 0n) !== (b < 0n) ? q - 1n : q;
}
/** ⌊√n⌋ for n ≥ 0. */
export function isqrt(n) {
  if (n < 0n) throw new RangeError("isqrt of a negative number");
  if (n < 2n) return n;
  // a float start when n fits in a double, a power of two above the root otherwise
  let x = n >> 1000n === 0n ? BigInt(Math.floor(Math.sqrt(Number(n)))) : 1n << BigInt((n.toString(2).length + 1) >> 1);
  if (x === 0n) x = 1n;
  // Newton from a float start; correct the last step either way
  for (;;) {
    const y = (x + n / x) >> 1n;
    if (y >= x - 1n && y <= x + 1n) {
      x = y;
      break;
    }
    x = y;
  }
  while (x * x > n) x -= 1n;
  while ((x + 1n) * (x + 1n) <= n) x += 1n;
  return x;
}
/** ⌊n^(1/k)⌋ for n ≥ 0. */
export function iroot(n, k) {
  k = BigInt(k);
  if (n < 2n) return n;
  if (k === 1n) return n;
  let x = 1n << (BigInt(n.toString(2).length) / k + 1n);
  for (;;) {
    const y = ((k - 1n) * x + n / x ** (k - 1n)) / k;
    if (y >= x) break;
    x = y;
  }
  while (x ** k > n) x -= 1n;
  while ((x + 1n) ** k <= n) x += 1n;
  return x;
}
export const isSquare = (n) => n >= 0n && isqrt(n) ** 2n === n;
/** v_p(n): how many times p divides n ≠ 0. */
export function val(n, p) {
  if (n === 0n) return Infinity;
  let k = 0;
  n = abs(n);
  while (n % p === 0n) {
    n /= p;
    k++;
  }
  return k;
}
/** aᵉ mod m for e ≥ 0, m ≥ 1. */
export function powmod(a, e, m) {
  let r = 1n % m;
  a = mod(a, m);
  while (e > 0n) {
    if (e & 1n) r = (r * a) % m;
    a = (a * a) % m;
    e >>= 1n;
  }
  return r;
}
/** { g, x, y } with a·x + b·y = g = gcd(a, b) ≥ 0. */
export function egcd(a, b) {
  let [r0, r1, x0, x1, y0, y1] = [a, b, 1n, 0n, 0n, 1n];
  while (r1) {
    const q = r0 / r1;
    [r0, r1, x0, x1, y0, y1] = [r1, r0 - q * r1, x1, x0 - q * x1, y1, y0 - q * y1];
  }
  return r0 < 0n ? { g: -r0, x: -x0, y: -y0 } : { g: r0, x: x0, y: y0 };
}
export function inverse(a, m) {
  const { g, x } = egcd(mod(a, m), m);
  if (g !== 1n) refuse("no-inverse", `${a} has no inverse mod ${m}: they share the factor ${g}`);
  return mod(x, m);
}

// The public forms: any whole number in (number, bigint, numeric string, Core
// Integer or whole Rational), a BigInt out. The bare functions above are the
// internal BigInt-only ones.
export const publicInt = {
  isqrt: (n) => isqrt(atLeast(n, 0n, "n")),
  iroot: (n, k) => iroot(atLeast(n, 0n, "n"), atLeast(k, 1n, "k")),
  gcd: (a, b) => gcd(big(a, "a"), big(b, "b")),
  lcm: (a, b) => lcm(big(a, "a"), big(b, "b")),
  egcd: (a, b) => egcd(big(a, "a"), big(b, "b")),
};
