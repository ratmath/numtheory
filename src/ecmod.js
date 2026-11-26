// Affine arithmetic on y² = x³ + ax + b modulo n, written to be strict: every
// division checks that the divisor is a unit mod n. When n is prime nothing
// can go wrong; when n is composite a failed division throws a NotInvertible
// carrying a proper factor of n. Points are [x, y]; null is the point at
// infinity. Used for point counting (n prime) and by the search for elliptic
// curve certificates (n a probable prime).

import { invMod, NotInvertible } from "./polymod.js";

const g = (a, b) => {
  if (a < 0n) a = -a;
  while (b) [a, b] = [b, a % b];
  return a;
};

export function ecAdd(P, Q, a, n) {
  if (P === null) return Q;
  if (Q === null) return P;
  const [x1, y1] = P;
  const [x2, y2] = Q;
  let lam;
  if (x1 === x2) {
    if ((y1 + y2) % n === 0n) return null;
    // same x, and y neither equal nor opposite: impossible modulo a prime
    if (y1 !== y2) throw new NotInvertible(g(y1 - y2, n), n);
    lam = ((3n * x1 * x1 + a) * invMod(2n * y1, n)) % n;
  } else {
    lam = ((y2 - y1) * invMod(x2 - x1, n)) % n;
  }
  let x3 = (lam * lam - x1 - x2) % n;
  if (x3 < 0n) x3 += n;
  let y3 = (lam * (x1 - x3) - y1) % n;
  if (y3 < 0n) y3 += n;
  return [x3, y3];
}
export const ecNeg = (P, n) => (P === null ? null : [P[0], P[1] === 0n ? 0n : n - P[1]]);

/** k·P for k ≥ 0. */
export function ecMul(k, P, a, n) {
  if (k < 0n) return ecMul(-k, ecNeg(P, n), a, n);
  let R = null;
  if (k === 0n || P === null) return R;
  const bits = k.toString(2);
  for (let i = 0; i < bits.length; i++) {
    R = ecAdd(R, R, a, n);
    if (bits[i] === "1") R = ecAdd(R, P, a, n);
  }
  return R;
}
