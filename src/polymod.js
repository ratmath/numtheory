// Polynomials over Z/n (n a prime, or a number being treated as one), as
// arrays of BigInts with the constant term first and no trailing zeros.
// Used by Schoof's algorithm and by the search for elliptic curve
// certificates. When n is not prime a division can meet a non-invertible
// leading coefficient; that throws a NotInvertible carrying the factor found.

export class NotInvertible extends Error {
  constructor(factor, n) {
    super(`${factor} divides the modulus ${n}`);
    this.factor = factor;
  }
}

/** The inverse of a mod n, or a NotInvertible carrying gcd(a, n). */
export function invMod(a, n) {
  a %= n;
  if (a < 0n) a += n;
  let [r0, r1, x0, x1] = [n, a, 0n, 1n];
  while (r1) {
    const q = r0 / r1;
    [r0, r1, x0, x1] = [r1, r0 - q * r1, x1, x0 - q * x1];
  }
  if (r0 !== 1n) throw new NotInvertible(r0, n);
  return x0 < 0n ? x0 + n : x0;
}

export const trim = (a) => {
  let k = a.length;
  while (k && a[k - 1] === 0n) k--;
  a.length = k;
  return a;
};
export const deg = (a) => a.length - 1;
export const pol = (cs, n) => trim(cs.map((c) => ((c % n) + n) % n));
export const eq = (a, b) => a.length === b.length && a.every((c, i) => c === b[i]);

export function add(a, b, n) {
  const out = new Array(Math.max(a.length, b.length));
  for (let i = 0; i < out.length; i++) {
    const s = (a[i] ?? 0n) + (b[i] ?? 0n);
    out[i] = s >= n ? s - n : s;
  }
  return trim(out);
}
export function sub(a, b, n) {
  const out = new Array(Math.max(a.length, b.length));
  for (let i = 0; i < out.length; i++) {
    const s = (a[i] ?? 0n) - (b[i] ?? 0n);
    out[i] = s < 0n ? s + n : s;
  }
  return trim(out);
}
export const scale = (a, c, n) => trim(a.map((x) => (x * c) % n));

export function mul(a, b, n) {
  if (!a.length || !b.length) return [];
  const out = new Array(a.length + b.length - 1).fill(0n);
  for (let i = 0; i < a.length; i++) {
    const ai = a[i];
    if (ai === 0n) continue;
    for (let j = 0; j < b.length; j++) out[i + j] += ai * b[j];
  }
  for (let i = 0; i < out.length; i++) out[i] %= n;
  return trim(out);
}

/** a = q·b + r with deg r < deg b. */
export function divmod(a, b, n) {
  if (!b.length) throw new RangeError("division by the zero polynomial");
  const r = a.slice();
  const db = deg(b);
  if (r.length <= db) return { q: [], r };
  const lead = b[db];
  const li = lead === 1n ? 1n : invMod(lead, n);
  const q = new Array(r.length - db).fill(0n);
  for (let i = r.length - 1; i >= db; i--) {
    const c = (r[i] * li) % n;
    q[i - db] = c;
    if (c === 0n) continue;
    for (let j = 0; j <= db; j++) {
      const k = i - db + j;
      r[k] = (r[k] - c * b[j]) % n;
    }
  }
  r.length = db;
  for (let i = 0; i < db; i++) if (r[i] < 0n) r[i] += n;
  return { q: trim(q), r: trim(r) };
}
export const rem = (a, b, n) => divmod(a, b, n).r;
export const mulmod = (a, b, h, n) => rem(mul(a, b, n), h, n);

/** The monic gcd. */
export function gcd(a, b, n) {
  while (b.length) [a, b] = [b, rem(a, b, n)];
  return a.length ? scale(a, invMod(a[a.length - 1], n), n) : a;
}

/** baseᵉ mod h. */
export function powmod(base, e, h, n) {
  let r = rem([1n], h, n);
  let b = rem(base, h, n);
  while (e > 0n) {
    if (e & 1n) r = mulmod(r, b, h, n);
    e >>= 1n;
    if (e > 0n) b = mulmod(b, b, h, n);
  }
  return r;
}

/** The inverse of a mod h; when it does not exist, throws { divisor } holding gcd(a, h) (monic, of positive degree). */
export class ZeroDivisor extends Error {
  constructor(divisor) {
    super("not invertible modulo the polynomial");
    this.divisor = divisor;
  }
}
export function inverse(a, h, n) {
  // extended Euclid on (h, a)
  let [r0, r1, s0, s1] = [h, rem(a, h, n), [], [1n]];
  while (r1.length) {
    const { q, r } = divmod(r0, r1, n);
    [r0, r1, s0, s1] = [r1, r, s1, sub(s0, mul(q, s1, n), n)];
  }
  if (deg(r0) !== 0) throw new ZeroDivisor(r0.length ? scale(r0, invMod(r0[r0.length - 1], n), n) : h);
  return rem(scale(s0, invMod(r0[0], n), n), h, n);
}

/** p(x) evaluated at x mod n. */
export function evaluate(a, x, n) {
  let r = 0n;
  for (let i = a.length - 1; i >= 0; i--) r = (r * x + a[i]) % n;
  return r;
}

/**
 * One root in Z/n of a polynomial that splits into distinct linear factors
 * (n prime), by equal-degree splitting; null when none is found.
 */
export function oneRoot(f, n) {
  let g = gcd(sub(powmod([0n, 1n], n, f, n), [0n, 1n], n), f, n);
  if (deg(g) < 1) return null;
  for (let r = 1n; deg(g) > 1 && r < 200n; r++) {
    const w = sub(powmod([r, 1n], (n - 1n) / 2n, g, n), [1n], n);
    const d = gcd(w, g, n);
    if (deg(d) >= 1 && deg(d) < deg(g)) g = deg(d) * 2 <= deg(g) ? d : divmod(g, d, n).q;
  }
  if (deg(g) !== 1) return null;
  return (((n - g[0]) % n) * invMod(g[1], n)) % n;
}
