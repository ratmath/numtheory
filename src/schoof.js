// Schoof's algorithm: the number of points of y² = x³ + Ax + B over F_p
// (p > 3) in time polynomial in log p.
//
// The Frobenius map φ(x, y) = (x^p, y^p) satisfies φ² − tφ + p = 0 on the
// curve, where #E(F_p) = p + 1 − t and |t| ≤ 2√p (Hasse). For a small prime l
// the equation is tested on the points of order l, which are the roots of the
// division polynomial ψ_l: t mod l is the one τ with
//     φ²(P) + [p mod l]P = [τ]φ(P)   for P of order l.
// Enough primes l determine t by the Chinese remainder theorem.
//
// Points of order l are handled all at once, as pairs (X(x), y·Y(x)) of
// polynomials modulo h = ψ_l. ψ_l need not be irreducible, so a division can
// meet a zero divisor; that hands over a proper factor of h, and the
// computation is run again modulo the factor. Any factor will do, because the
// equation holds on every point of order l and fixes τ on each one.

import * as P from "./polymod.js";

/** The polynomials f_0 … f_N with ψ_n = f_n (n odd), ψ_n = 2y·f_n (n even), mod p. */
export function divisionPolynomials(A, B, p, N) {
  const m = (c) => ((c % p) + p) % p;
  const mul = (a, b) => P.mul(a, b, p);
  const f = [[], [1n], [1n]];
  f[3] = P.pol([-A * A, 12n * B, 6n * A, 0n, 3n], p);
  f[4] = P.pol([2n * (-8n * B * B - A * A * A), -8n * A * B, -10n * A * A, 40n * B, 10n * A, 0n, 2n], p);
  const F = P.pol([4n * m(B), 4n * m(A), 0n, 4n], p);
  const F2 = mul(F, F);
  const cube = (a) => mul(a, mul(a, a));
  for (let n = 5; n <= N; n++) {
    const k = n >> 1;
    if (n % 2) {
      const a = mul(f[k + 2], cube(f[k]));
      const b = mul(f[k - 1], cube(f[k + 1]));
      f[n] = k % 2 === 0 ? P.sub(mul(F2, a), b, p) : P.sub(a, mul(F2, b), p);
    } else {
      f[n] = mul(f[k], P.sub(mul(f[k + 2], mul(f[k - 1], f[k - 1])), mul(f[k - 2], mul(f[k + 1], f[k + 1])), p));
    }
  }
  return f.slice(0, N + 1);
}

/** t mod l for an odd prime l ≠ p, working modulo the monic polynomial h (a factor of ψ_l). */
function traceModL(A, B, p, l, h) {
  const mulm = (a, b) => P.mulmod(a, b, h, p);
  const fx = P.rem(P.pol([B, A, 0n, 1n], p), h, p);
  const neg = (a) => P.sub([], a, p);
  const add = (S, T) => {
    if (S === null) return T;
    if (T === null) return S;
    let L;
    if (P.eq(S.X, T.X)) {
      if (!P.eq(S.Y, T.Y)) {
        if (P.eq(S.Y, neg(T.Y))) return null;
        throw new P.ZeroDivisor(P.gcd(P.sub(S.Y, T.Y, p), h, p));
      }
      if (!S.Y.length) return null;
      const den = P.scale(mulm(fx, S.Y), 2n, p);
      const num = P.add(P.scale(mulm(S.X, S.X), 3n, p), P.pol([A], p), p);
      L = mulm(num, P.inverse(den, h, p));
    } else {
      L = mulm(P.sub(T.Y, S.Y, p), P.inverse(P.sub(T.X, S.X, p), h, p));
    }
    const X = P.sub(P.sub(mulm(fx, mulm(L, L)), S.X, p), T.X, p);
    const Y = P.sub(mulm(L, P.sub(S.X, X, p)), S.Y, p);
    return { X, Y };
  };
  const same = (S, T) => (S === null || T === null ? S === T : P.eq(S.X, T.X) && P.eq(S.Y, T.Y));
  const xp = P.powmod([0n, 1n], p, h, p);
  const yp = P.powmod(fx, (p - 1n) / 2n, h, p);
  const xpp = P.powmod(xp, p, h, p);
  const ypp = mulm(P.powmod(yp, p, h, p), yp);
  // [p mod l](x, y)
  let Q = null;
  const base = { X: P.rem([0n, 1n], h, p), Y: P.rem([1n], h, p) };
  const k = Number(p % BigInt(l));
  for (const bit of k.toString(2)) {
    Q = add(Q, Q);
    if (bit === "1") Q = add(Q, base);
  }
  const S = add({ X: xpp, Y: ypp }, Q);
  const phi = { X: xp, Y: yp };
  let T = null;
  for (let tau = 0; tau < l; tau++) {
    if (same(S, T)) return tau;
    T = add(T, phi);
  }
  throw new Error("Schoof: no trace matched; this should not happen");
}

/** The trace of Frobenius t of y² = x³ + Ax + B over F_p (p > 3 prime), so that #E = p + 1 − t. Option `log` collects { l, t } per prime. */
export function schoofTrace(A, B, p, { log } = {}) {
  const f = P.pol([B, A, 0n, 1n], p);
  // l = 2: a point of order 2 exists exactly when f has a root in F_p
  const g = P.gcd(P.sub(P.powmod([0n, 1n], p, f, p), [0n, 1n], p), f, p);
  let t = P.deg(g) >= 1 ? 0n : 1n;
  let M = 2n;
  log?.push({ l: 2, t: Number(t) });
  const ls = [];
  for (let l = 3n, prod = 2n; prod * prod <= 16n * p; l += 2n) {
    if (l === p || [3n, 5n, 7n, 11n, 13n].some((q) => q < l && l % q === 0n) || isComposite(l)) continue;
    ls.push(Number(l));
    prod *= l;
  }
  const psi = ls.length ? divisionPolynomials(A, B, p, ls[ls.length - 1]) : [];
  for (const l of ls) {
    let h = psi[l];
    h = P.scale(h, P.invMod(h[h.length - 1], p), p);
    let tl;
    for (;;) {
      try {
        tl = traceModL(A, B, p, l, h);
        break;
      } catch (e) {
        if (!(e instanceof P.ZeroDivisor) || P.deg(e.divisor) < 1 || P.deg(e.divisor) >= P.deg(h)) throw e;
        h = e.divisor;
      }
    }
    log?.push({ l, t: tl });
    // CRT: t ≡ t (mod M), t ≡ tl (mod l)
    const L = BigInt(l);
    const inv = P.invMod(M % L, L);
    t = t + M * ((((BigInt(tl) - t) % L) + L) * inv % L);
    M *= L;
  }
  t %= M;
  if (t > M / 2n) t -= M;
  return t;
}
function isComposite(l) {
  for (let d = 3n; d * d <= l; d += 2n) if (l % d === 0n) return true;
  return false;
}
