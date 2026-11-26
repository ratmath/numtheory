// An independent checker for primality certificates.
//
// This file imports nothing from the rest of the package: it has its own
// modular power, gcd and integer square root, so that a mistake in the code
// that finds a certificate cannot hide a mistake in the code that checks it.
//
// The theorems used:
//   small        n < 2^32 is checked by trial division up to √n.
//   pocklington  n − 1 = F·R with every prime q | F certified and a base a with
//                a^(n−1) ≡ 1 and gcd(a^((n−1)/q) − 1, n) = 1 (mod n). Then every
//                prime factor p of n is ≡ 1 (mod F), so p ≥ F + 1, and when
//                (F + 1)² > n, n is prime.
//   bls          the same conditions with F³ > n ≥ (F + 1)²; write
//                n − 1 = F(c2·F + c1), 0 ≤ c1 < F. Then n is prime exactly when
//                c1² − 4·c2 is not a square (Brillhart, Lehmer, Selfridge 1975).
//   elliptic     gcd(n, 6) = 1, a curve y² = x³ + ax + b with 4a³ + 27b² a unit
//                mod n, a point P on it, and m = k·q with q a certified prime
//                above (n^(1/4) + 1)², such that U = k·P is a finite point and
//                q·U is the point at infinity. If n had a prime factor
//                p ≤ √n, then U would have order q on the curve mod p, which
//                has at most (√p + 1)² ≤ (n^(1/4) + 1)² < q points (Hasse):
//                impossible, so n is prime (Goldwasser and Kilian 1986).
//                The curve arithmetic is done with affine coordinates and
//                every divisor is checked to be a unit mod n, so that each
//                step is the true group operation modulo every prime factor
//                of n at once.

const pw = (a, e, m) => {
  let r = 1n % m;
  a %= m;
  if (a < 0n) a += m;
  while (e > 0n) {
    if (e & 1n) r = (r * a) % m;
    a = (a * a) % m;
    e >>= 1n;
  }
  return r;
};
const g = (a, b) => {
  if (a < 0n) a = -a;
  if (b < 0n) b = -b;
  while (b) [a, b] = [b, a % b];
  return a;
};
const sq = (n) => {
  if (n < 2n) return n;
  let x = n;
  let y = (x + 1n) / 2n;
  while (y < x) {
    x = y;
    y = (x + n / x) / 2n;
  }
  return x;
};

// k·P on y² = x³ + ax + b mod n, strictly: null is the point at infinity, and
// a divisor that is not a unit mod n makes the certificate invalid.
function curveAdd(P, Q, a, n) {
  if (P === null) return Q;
  if (Q === null) return P;
  const [x1, y1] = P;
  const [x2, y2] = Q;
  let num;
  let den;
  if (x1 === x2) {
    if ((y1 + y2) % n === 0n) return null;
    if (y1 !== y2) fail("two points with the same x and unrelated y: n is composite");
    num = 3n * x1 * x1 + a;
    den = 2n * y1;
  } else {
    num = y2 - y1;
    den = x2 - x1;
  }
  den = ((den % n) + n) % n;
  // the inverse of den, by Euclid; it must exist
  let [r0, r1, s0, s1] = [n, den, 0n, 1n];
  while (r1) {
    const q = r0 / r1;
    [r0, r1, s0, s1] = [r1, r0 - q * r1, s1, s0 - q * s1];
  }
  if (r0 !== 1n) fail("a divisor in the curve arithmetic shares a factor with n");
  const lam = (((num * s0) % n) + n) % n;
  const x3 = (((lam * lam - x1 - x2) % n) + n) % n;
  const y3 = (((lam * (x1 - x3) - y1) % n) + n) % n;
  return [x3, y3];
}
function curveMul(k, P, a, n) {
  let R = null;
  for (const bit of k.toString(2)) {
    R = curveAdd(R, R, a, n);
    if (bit === "1") R = curveAdd(R, P, a, n);
  }
  return R;
}
function checkElliptic(c, depth) {
  const { n, a, b, x, y, m, q } = c;
  for (const v of [a, b, x, y, m, q]) if (typeof v !== "bigint") fail("an elliptic certificate needs BigInt a, b, x, y, m and q");
  if (g(n, 6n) !== 1n) fail("n is divisible by 2 or 3");
  for (const v of [a, b, x, y]) if (v < 0n || v >= n) fail("a, b, x and y must be residues mod n");
  if ((y * y - (x * x * x + a * x + b)) % n !== 0n) fail("the point is not on the curve");
  if (g((4n * a * a * a + 27n * b * b) % n, n) !== 1n) fail("the discriminant of the curve is not a unit mod n");
  if (q < 2n || m % q !== 0n) fail("q does not divide m");
  const s = sq(q);
  if (s < 2n || (s - 1n) ** 4n <= n) fail("q does not exceed (n^(1/4) + 1)²");
  if (!c.certificate || c.certificate.n !== q) fail(`the prime ${q} has no certificate of its own`);
  check(c.certificate, depth + 1);
  const U = curveMul(m / q, [x, y], a, n);
  if (U === null) fail("(m/q)·P is the point at infinity");
  if (curveMul(q, U, a, n) !== null) fail("m·P is not the point at infinity");
}

/**
 * Check a primality certificate from provePrime. Returns
 * { valid: true } or { valid: false, reason }. Nothing is taken on trust:
 * every prime in the certificate is itself checked.
 */
export function verifyCertificate(cert) {
  try {
    check(cert, 0);
    return { valid: true };
  } catch (e) {
    if (e instanceof Invalid) return { valid: false, reason: e.message };
    throw e;
  }
}
class Invalid extends Error {}
const fail = (m) => {
  throw new Invalid(m);
};

function check(c, depth) {
  if (!c || typeof c !== "object" || typeof c.n !== "bigint") fail("a certificate is an object with a BigInt n");
  if (depth > 2000) fail("the certificate is nested too deeply");
  const n = c.n;
  if (n < 2n) fail(`${n} is below 2`);
  if (c.kind === "small") {
    if (n >= 1n << 32n) fail(`${n} is too large for a trial-division certificate`);
    for (let p = 2n; p * p <= n; p += p === 2n ? 1n : 2n) if (n % p === 0n) fail(`${p} divides ${n}`);
    return;
  }
  if (c.kind === "elliptic") return checkElliptic(c, depth);
  if (c.kind !== "pocklington" && c.kind !== "bls") fail(`unknown kind ${c.kind}`);
  const m = n - 1n;
  if (!Array.isArray(c.factors) || !c.factors.length) fail("no factors of n − 1 are given");
  let F = 1n;
  for (const f of c.factors) {
    const { q, e, a } = f;
    if (typeof q !== "bigint" || typeof a !== "bigint" || !Number.isInteger(e) || e < 1) fail("each factor needs q, e and a");
    if (!f.certificate || f.certificate.n !== q) fail(`the factor ${q} has no certificate of its own`);
    check(f.certificate, depth + 1);
    F *= q ** BigInt(e);
    if (pw(a, m, n) !== 1n) fail(`${a}^(n−1) is not 1 mod n`);
    if (g(pw(a, m / q, n) - 1n, n) !== 1n) fail(`gcd(${a}^((n−1)/${q}) − 1, n) is not 1`);
  }
  if (F !== c.F) fail("F is not the product of the certified prime powers");
  if (m % F !== 0n || m / F !== c.R) fail("n − 1 is not F·R");
  if (c.kind === "pocklington") {
    if ((F + 1n) * (F + 1n) <= n) fail("(F + 1)² does not exceed n");
    return;
  }
  if (F * F * F <= n) fail("F³ does not exceed n");
  const { c1, c2 } = c;
  if (typeof c1 !== "bigint" || typeof c2 !== "bigint" || c1 < 0n || c1 >= F || c2 < 1n) fail("c1 and c2 are out of range");
  if (c.R !== c2 * F + c1) fail("R is not c2·F + c1");
  const d = c1 * c1 - 4n * c2;
  if (d >= 0n && sq(d) ** 2n === d) fail("c1² − 4·c2 is a square, so n is composite");
}

/** Check that a proof of compositeness holds: a proper factor, or a Miller–Rabin witness. */
export function verifyComposite(n, { factor, witness } = {}) {
  if (typeof n !== "bigint" || n < 4n) return { valid: false, reason: "n must be a BigInt ≥ 4" };
  if (factor !== undefined) {
    if (factor > 1n && factor < n && n % factor === 0n) return { valid: true };
    return { valid: false, reason: `${factor} is not a proper factor of ${n}` };
  }
  if (witness !== undefined) {
    if (n % 2n === 0n) return { valid: true };
    let d = n - 1n;
    let s = 0;
    while (d % 2n === 0n) {
      d /= 2n;
      s++;
    }
    let x = pw(witness, d, n);
    if (x === 1n || x === n - 1n) return { valid: false, reason: `${witness} is not a witness` };
    for (let i = 1; i < s; i++) {
      x = (x * x) % n;
      if (x === n - 1n) return { valid: false, reason: `${witness} is not a witness` };
    }
    return { valid: true };
  }
  return { valid: false, reason: "give a factor or a witness" };
}
