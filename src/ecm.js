// Lenstra's elliptic curve method of factoring.
//
// Pollard's p − 1 method finds a prime factor p of n when the order p − 1 of
// the group (Z/p)* happens to be smooth. The elliptic curve method replaces
// that one group by the groups E(F_p) of many curves: their orders are spread
// over the Hasse interval around p + 1, so if one curve is unlucky the next
// may not be. The work depends on the size of the factor found, not on n.
//
// Curves are in Montgomery form By² = x³ + Ax² + x with Suyama's
// parametrisation by σ (which gives a group order divisible by 12), and
// arithmetic uses the x-coordinate only. Stage 1 multiplies a point by every
// prime power up to B1; stage 2 allows one more prime between B1 and B2.
// A factor shows up as a gcd with n.

import { big, gcd, mod, refuse } from "./int.js";

const sieveCache = { limit: 0, flags: null };
function primeFlags(limit) {
  if (sieveCache.limit >= limit) return sieveCache.flags;
  const f = new Uint8Array(limit + 1).fill(1);
  f[0] = f[1] = 0;
  for (let i = 2; i * i <= limit; i++) if (f[i]) for (let j = i * i; j <= limit; j += i) f[j] = 0;
  sieveCache.limit = limit;
  sieveCache.flags = f;
  return f;
}

/**
 * One curve of the elliptic curve method on n: Suyama's curve for σ ≥ 6, with
 * bounds B1 and B2. Returns { factor, stage, sigma, B1, B2 }: factor is a
 * proper factor of n or null, and stage says where it appeared (0: while
 * setting up the curve, 1 or 2).
 */
export function ecmCurve(n, sigma, { B1 = 2000, B2 = 100 * B1 } = {}) {
  n = big(n, "n");
  const s = big(sigma, "sigma");
  if (n < 2n || s < 6n) throw new RangeError("ecmCurve needs n ≥ 2 and σ ≥ 6");
  if (!Number.isSafeInteger(B1) || !Number.isSafeInteger(B2) || B1 < 2 || B2 < B1 || B2 > 200000000) throw new RangeError("the bounds must satisfy 2 ≤ B1 ≤ B2 ≤ 2·10⁸");
  const done = (factor, stage) => ({ factor, stage, sigma: s, B1, B2 });
  const proper = (g) => g > 1n && g < n;
  const u = mod(s * s - 5n, n);
  const v = mod(4n * s, n);
  // (A + 2)/4 = (v − u)³(3u + v) / (16u³v)
  const den = mod(16n * u * u * u * v, n);
  const g0 = gcd(den, n);
  if (g0 !== 1n) return done(proper(g0) ? g0 : null, 0);
  const a24 = mod((v - u) ** 3n * (3n * u + v) * inverse(den, n), n);
  const dbl = ([X, Z]) => {
    const a = (X + Z) * (X + Z) % n;
    const b = (X - Z) * (X - Z) % n;
    const c = a - b;
    return [(a * b) % n, mod(c * (b + a24 * c), n)];
  };
  // P + Q from P, Q and P − Q
  const add = ([X1, Z1], [X2, Z2], [X0, Z0]) => {
    const p = (X1 - Z1) * (X2 + Z2) % n;
    const q = (X1 + Z1) * (X2 - Z2) % n;
    return [mod(Z0 * (p + q) * (p + q), n), mod(X0 * (p - q) * (p - q), n)];
  };
  const ladder = (k, P) => {
    if (k === 1n) return P;
    let R0 = P;
    let R1 = dbl(P);
    const bits = k.toString(2);
    for (let i = 1; i < bits.length; i++) {
      if (bits[i] === "1") {
        R0 = add(R0, R1, P);
        R1 = dbl(R1);
      } else {
        R1 = add(R0, R1, P);
        R0 = dbl(R0);
      }
    }
    return R0;
  };
  const flags = primeFlags(B2 + 420);
  let Q = [(u * u * u) % n, (v * v * v) % n];
  for (let p = 2; p <= B1; p++) {
    if (!flags[p]) continue;
    let pk = p;
    while (pk * p <= B1) pk *= p;
    Q = ladder(BigInt(pk), Q);
  }
  let g = gcd(Q[1], n);
  if (g !== 1n) return done(proper(g) ? g : null, 1);
  if (B2 === B1) return done(null, null);
  // stage 2: primes q = i·D ± j in (B1, B2], comparing x(iD·Q) with x(j·Q)
  const D = 210;
  const S = new Array(D / 2 + 1);
  S[1] = Q;
  S[2] = dbl(Q);
  S[3] = add(S[2], S[1], S[1]);
  for (let j = 5; j <= D / 2; j += 2) S[j] = add(S[j - 2], S[2], S[j - 4]);
  const js = [];
  for (let j = 1; j <= D / 2; j += 2) if (j % 3 && j % 5 && j % 7) js.push(j);
  const T = ladder(BigInt(D), Q);
  let i = Math.max(1, Math.floor(B1 / D));
  let prev = i === 1 ? null : ladder(BigInt(i - 1), T);
  let cur = ladder(BigInt(i), T);
  let acc = 1n;
  for (; i * D - D / 2 <= B2; i++) {
    const [XR, ZR] = cur;
    const base = i * D;
    for (const j of js) {
      if ((base - j > B1 && flags[base - j]) || (base + j <= B2 && flags[base + j])) {
        acc = (acc * mod(XR * S[j][1] - S[j][0] * ZR, n)) % n;
      }
    }
    const next = prev === null ? dbl(T) : add(cur, T, prev);
    prev = cur;
    cur = next;
  }
  g = gcd(acc, n);
  return done(proper(g) ? g : null, proper(g) ? 2 : null);
}
function inverse(a, n) {
  let [r0, r1, x0, x1] = [n, a, 0n, 1n];
  while (r1) {
    const q = r0 / r1;
    [r0, r1, x0, x1] = [r1, r0 - q * r1, x1, x0 - q * x1];
  }
  return mod(x0, n);
}

/**
 * The elliptic curve method: a proper factor of n, or null when `curves`
 * curves with bounds B1 and B2 (default 100·B1) find none. The curves are
 * σ = sigma, sigma + 1, …, so a run can be repeated. As a guide, B1 = 2000
 * with 25 curves finds most 15-digit factors, B1 = 11000 with 90 curves
 * 20-digit ones, and B1 = 50000 with 300 curves 25-digit ones. Option
 * details: true returns the report of the successful curve instead.
 */
export function ecm(n, { B1 = 2000, B2 = 100 * B1, curves = 40, sigma = 6n, details = false } = {}) {
  n = big(n, "n");
  if (n < 4n) refuse("bad-argument", "nothing to factor below 4");
  for (const p of [2n, 3n]) if (n % p === 0n && n !== p) return details ? { factor: p, stage: 0, sigma: null, B1, B2 } : p;
  let s = big(sigma, "sigma");
  if (s < 6n) throw new RangeError("σ starts at 6");
  for (let c = 0; c < curves; c++, s++) {
    const r = ecmCurve(n, s, { B1, B2 });
    if (r.factor !== null) return details ? r : r.factor;
  }
  return null;
}
