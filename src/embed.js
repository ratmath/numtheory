// Numerical side of the number-field code: the roots of the defining
// polynomial in fixed-point BigInt arithmetic, logarithms to match, and two
// lattice reductions (an exact integral LLL for small lattices, a
// floating-point one for the wide relation lattice whose output is checked
// exactly by the caller).

// ---- fixed point: a real x is the BigInt round(x·2^bits) ----

const mulF = (a, b, bits) => (a * b) >> bits;

/** Roots of a monic integer polynomial T (constant first): { real: [x], complex: [[re, im]] } in fixed point, im > 0. */
export function rootsFixed(T, r1, bits) {
  const n = T.length - 1;
  const c = T.map(Number);
  // Durand–Kerner in doubles
  const bound = 1 + Math.max(...c.slice(0, n).map(Math.abs));
  let z = Array.from({ length: n }, (_, k) => {
    const a = 0.4 + (2 * Math.PI * k) / n;
    const r = Math.min(bound, 1 + Math.abs(c[0]) ** (1 / n));
    return [r * Math.cos(a), r * Math.sin(a)];
  });
  const evalC = ([x, y]) => {
    let re = 1;
    let im = 0;
    for (let i = n - 1; i >= 0; i--) {
      const t = re * x - im * y + c[i];
      im = re * y + im * x;
      re = t;
    }
    return [re, im];
  };
  for (let it = 0; it < 2000; it++) {
    let move = 0;
    z = z.map((zi, i) => {
      let [nr, ni] = evalC(zi);
      let dr = 1;
      let di = 0;
      for (let j = 0; j < n; j++) {
        if (j === i) continue;
        const a = zi[0] - z[j][0];
        const b = zi[1] - z[j][1];
        const t = dr * a - di * b;
        di = dr * b + di * a;
        dr = t;
      }
      const den = dr * dr + di * di || 1e-300;
      const qr = (nr * dr + ni * di) / den;
      const qi = (ni * dr - nr * di) / den;
      move = Math.max(move, Math.abs(qr) + Math.abs(qi));
      return [zi[0] - qr, zi[1] - qi];
    });
    if (move < 1e-15 * bound) break;
  }
  z.sort((a, b) => Math.abs(a[1]) - Math.abs(b[1]));
  const one = 1n << BigInt(bits);
  const B = BigInt(bits);
  const toF = (x) => {
    // a double as fixed point
    const s = Math.round(x * 2 ** 50);
    return (BigInt(s) << B) >> 50n;
  };
  const Tb = T;
  // Newton in fixed-point complex arithmetic
  const newton = (re, im, realOnly) => {
    for (let it = 0; it < 12 + Math.ceil(Math.log2(bits)); it++) {
      // T and T′ by Horner
      let [pr, pi] = [one, 0n];
      let [qr, qi] = [0n, 0n];
      for (let i = n - 1; i >= 0; i--) {
        // q = q·z + p ; p = p·z + T[i]
        const nqr = mulF(qr, re, B) - mulF(qi, im, B) + pr;
        const nqi = mulF(qr, im, B) + mulF(qi, re, B) + pi;
        qr = nqr;
        qi = nqi;
        const npr = mulF(pr, re, B) - mulF(pi, im, B) + Tb[i] * one;
        const npi = mulF(pr, im, B) + mulF(pi, re, B);
        pr = npr;
        pi = npi;
      }
      const den = mulF(qr, qr, B) + mulF(qi, qi, B);
      if (den === 0n) break;
      const wr = ((mulF(pr, qr, B) + mulF(pi, qi, B)) << B) / den;
      const wi = ((mulF(pi, qr, B) - mulF(pr, qi, B)) << B) / den;
      re -= wr;
      im = realOnly ? 0n : im - wi;
      const tiny = (x) => (x < 0n ? -x : x) < 8n;
      if (tiny(wr) && (realOnly || tiny(wi))) break;
    }
    return [re, im];
  };
  const real = z.slice(0, r1).map(([x]) => newton(toF(x), 0n, true)[0]);
  real.sort((a, b) => (a < b ? -1 : a > b ? 1 : 0));
  const complex = z
    .slice(r1)
    .filter(([, y]) => y > 0)
    .map(([x, y]) => newton(toF(x), toF(y), false));
  if (complex.length !== (n - r1) / 2) throw new Error("the numerical roots did not separate into real and complex ones");
  return { real, complex };
}

const ln2cache = new Map();
function atanhF(t, bits) {
  // t + t³/3 + … for |t| ≤ 1/3
  const B = BigInt(bits);
  const t2 = mulF(t, t, B);
  let term = t;
  let sum = t;
  for (let k = 3n; ; k += 2n) {
    term = mulF(term, t2, B);
    if (term === 0n) break;
    sum += term / k;
  }
  return sum;
}
/** ln of a positive fixed-point number. */
export function lnFixed(x, bits) {
  if (x <= 0n) throw new RangeError("logarithm of a number that is not positive (more precision is needed)");
  const B = BigInt(bits);
  const one = 1n << B;
  if (!ln2cache.has(bits)) ln2cache.set(bits, 2n * atanhF(one / 3n, bits));
  const k = BigInt(x.toString(2).length) - B - 1n;
  const m = k >= 0n ? x >> k : x << -k; // in [1, 2)
  const t = ((m - one) << B) / (m + one);
  return 2n * atanhF(t, bits) + k * ln2cache.get(bits);
}
export const fixedToNumber = (x, bits) => {
  const neg = x < 0n;
  const a = neg ? -x : x;
  const shift = Math.max(0, a.toString(2).length - 60);
  const v = Number(a >> BigInt(shift)) * 2 ** (shift - bits);
  return neg ? -v : v;
};

// ---- exact integral LLL (Cohen, Algorithm 2.6.7), rows linearly independent ----

/** LLL-reduce linearly independent integer rows (BigInt), δ = 3/4; returns the reduced rows. */
export function lllInt(rows) {
  const b = rows.map((r) => r.slice());
  const n = b.length;
  if (n < 2) return b;
  const dot = (u, v) => {
    let s = 0n;
    for (let i = 0; i < u.length; i++) s += u[i] * v[i];
    return s;
  };
  const d = new Array(n + 1).fill(0n); // d[0] = 1, d[i] = det of the first i Gram rows
  d[0] = 1n;
  const lam = Array.from({ length: n }, () => new Array(n).fill(0n));
  const roundDiv = (a, m) => {
    // nearest integer to a/m, m > 0
    const t = 2n * a + m;
    const q = t / (2n * m);
    return t % (2n * m) !== 0n && t < 0n ? q - 1n : q;
  };
  d[1] = dot(b[0], b[0]);
  let k = 1;
  let kmax = 0;
  const red = (k, l) => {
    if (2n * (lam[k][l] < 0n ? -lam[k][l] : lam[k][l]) <= d[l + 1]) return;
    const q = roundDiv(lam[k][l], d[l + 1]);
    for (let i = 0; i < b[k].length; i++) b[k][i] -= q * b[l][i];
    lam[k][l] -= q * d[l + 1];
    for (let i = 0; i < l; i++) lam[k][i] -= q * lam[l][i];
  };
  while (k < n) {
    if (k > kmax) {
      kmax = k;
      for (let j = 0; j <= k; j++) {
        let u = dot(b[k], b[j]);
        for (let i = 0; i < j; i++) u = (d[i + 1] * u - lam[k][i] * lam[j][i]) / d[i];
        if (j < k) lam[k][j] = u;
        else {
          d[k + 1] = u;
          if (u === 0n) throw new RangeError("lllInt needs linearly independent rows");
        }
      }
    }
    red(k, k - 1);
    if (4n * d[k + 1] * d[k - 1] < 3n * d[k] * d[k] - 4n * lam[k][k - 1] * lam[k][k - 1]) {
      // swap k and k − 1
      const l = lam[k][k - 1];
      const B = (d[k - 1] * d[k + 1] + l * l) / d[k];
      [b[k], b[k - 1]] = [b[k - 1], b[k]];
      for (let j = 0; j < k - 1; j++) [lam[k][j], lam[k - 1][j]] = [lam[k - 1][j], lam[k][j]];
      for (let i = k + 1; i <= kmax; i++) {
        const t = lam[i][k];
        lam[i][k] = (d[k + 1] * lam[i][k - 1] - l * t) / d[k];
        lam[i][k - 1] = (B * t + l * lam[i][k]) / d[k + 1];
      }
      d[k] = B;
      k = Math.max(1, k - 1);
    } else {
      for (let l = k - 2; l >= 0; l--) red(k, l);
      k++;
    }
  }
  return b;
}

// ---- floating-point LLL for wide lattices of small integers ----

/**
 * LLL on rows of small JavaScript integers, Gram–Schmidt in doubles. Every
 * step is an exact unimodular row operation, so the result spans the same
 * lattice whatever the rounding does; only how reduced it is depends on the
 * floating point.
 */
export function lllFloat(rows, delta = 0.99) {
  const b = rows.map((r) => Float64Array.from(r));
  const n = b.length;
  if (n < 2) return b.map((r) => Array.from(r));
  const len = b[0].length;
  const mu = Array.from({ length: n }, () => new Float64Array(n));
  const Bn = new Float64Array(n);
  const dot = (u, v) => {
    let s = 0;
    for (let i = 0; i < len; i++) s += u[i] * v[i];
    return s;
  };
  const gso = (i) => {
    for (let j = 0; j < i; j++) {
      let s = dot(b[i], b[j]);
      for (let l = 0; l < j; l++) s -= mu[j][l] * mu[i][l] * Bn[l];
      mu[i][j] = Bn[j] > 0 ? s / Bn[j] : 0;
    }
    let s = dot(b[i], b[i]);
    for (let j = 0; j < i; j++) s -= mu[i][j] * mu[i][j] * Bn[j];
    Bn[i] = s;
  };
  gso(0);
  let k = 1;
  let guard = 0;
  while (k < n) {
    if (++guard > 4e6) break;
    gso(k);
    for (let pass = 0; pass < 4; pass++) {
      let changed = false;
      for (let j = k - 1; j >= 0; j--) {
        const q = Math.round(mu[k][j]);
        if (q === 0) continue;
        changed = true;
        const bk = b[k];
        const bj = b[j];
        for (let i = 0; i < len; i++) bk[i] -= q * bj[i];
        for (let l = 0; l < j; l++) mu[k][l] -= q * mu[j][l];
        mu[k][j] -= q;
      }
      if (!changed) break;
      gso(k);
    }
    if (Bn[k] >= (delta - mu[k][k - 1] * mu[k][k - 1]) * Bn[k - 1]) k++;
    else {
      [b[k], b[k - 1]] = [b[k - 1], b[k]];
      k = Math.max(1, k - 1);
      if (k === 1) gso(0);
    }
  }
  return b.map((r) => Array.from(r));
}
