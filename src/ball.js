// Rigorous real and complex ball arithmetic in fixed point, for the values of
// zeta and L functions. A real ball is a midpoint m and a radius r, both
// BigInts in units of 2^−p: it stands for every real within r of m. Each
// operation rounds its midpoint and adds a bound for that rounding to the
// radius, so the true value of an expression always lies in its ball.

const abs = (x) => (x < 0n ? -x : x);

export class Ctx {
  constructor(p) {
    this.p = BigInt(p);
    this.bits = p;
    this.one = 1n << this.p;
    this._pi = null;
    this._ln2 = null;
  }
  // ---- real balls: [m, r] ----
  int(n) {
    return [BigInt(n) << this.p, 0n];
  }
  /** The ball of a fraction a/b. */
  frac(a, b) {
    if (b < 0n) [a, b] = [-a, -b];
    const m = (a << this.p) / b;
    return [m, (a << this.p) % b === 0n ? 0n : 1n];
  }
  add([a, r], [b, s]) {
    return [a + b, r + s];
  }
  sub([a, r], [b, s]) {
    return [a - b, r + s];
  }
  neg([a, r]) {
    return [-a, r];
  }
  mul([a, r], [b, s]) {
    return [(a * b) >> this.p, ((abs(a) * s + abs(b) * r + r * s) >> this.p) + 2n];
  }
  /** Multiply by an exact integer. */
  scale([a, r], n) {
    return [a * n, r * abs(n)];
  }
  inv([a, r]) {
    const lo = abs(a) - r;
    if (lo <= 0n) throw new RangeError("division by a ball that contains zero");
    const m = (this.one * this.one) / a;
    // |1/x − 1/a| ≤ r / (|a|·lo)
    return [m, (r * this.one * this.one) / (abs(a) * lo) + 2n];
  }
  div(x, y) {
    return this.mul(x, this.inv(y));
  }
  /** A bound for |x|, as a fixed-point number. */
  mag([a, r]) {
    return abs(a) + r;
  }
  // ---- constants ----
  ln2() {
    return (this._ln2 ??= this.scale(this.#atanhInv(3n), 2n));
  }
  pi() {
    if (!this._pi) {
      // Machin: π = 16·atan(1/5) − 4·atan(1/239)
      this._pi = this.sub(this.scale(this.#atanInv(5n), 16n), this.scale(this.#atanInv(239n), 4n));
    }
    return this._pi;
  }
  #atanInv(k) {
    // atan(1/k) = Σ (−1)^j / ((2j+1) k^(2j+1))
    let term = this.one / k;
    let sum = term;
    let n = 1n;
    let count = 1n;
    const k2 = k * k;
    for (let sign = -1n; term > 0n; sign = -sign) {
      term /= k2;
      n += 2n;
      sum += (sign * term) / n;
      count++;
    }
    return [sum, 2n * count + 2n];
  }
  #atanhInv(k) {
    // atanh(1/k) = Σ 1 / ((2j+1) k^(2j+1))
    let term = this.one / k;
    let sum = term;
    let n = 1n;
    let count = 1n;
    const k2 = k * k;
    while (term > 0n) {
      term /= k2;
      n += 2n;
      sum += term / n;
      count++;
    }
    return [sum, 2n * count + 2n];
  }
  // ---- elementary functions ----
  exp([a, r]) {
    if (r > this.one) throw new RangeError("exp of a ball wider than 1: more precision is needed");
    const [l2] = this.ln2();
    const k = (2n * a + l2) / (2n * l2) - (a < 0n && (2n * a + l2) % (2n * l2) !== 0n ? 1n : 0n); // round(a / ln 2)
    // t = a − k·ln2, |t| ≤ 0.35 + tiny
    const t = this.sub([a, 0n], this.scale(this.ln2(), k));
    // Taylor at the midpoint of t
    let term = this.one;
    let sum = this.one;
    let count = 1n;
    for (let j = 1n; ; j++) {
      term = (term * t[0]) >> this.p;
      term /= j;
      if (term === 0n || term === -1n) break;
      sum += term;
      count++;
    }
    // rounding: 2 per term; the truncated tail is below one unit
    let out = [sum, 4n * count + 16n];
    // widen for the radius of t and of the input: e^(±ρ) − 1 ≤ ρ + ρ² for ρ ≤ 1
    const rho = t[1] + r;
    if (rho > 0n) out = this.mul(out, [this.one, rho + ((rho * rho) >> this.p) + 1n]);
    if (k >= 0n) return [out[0] << k, out[1] << k];
    return [out[0] >> -k, (out[1] >> -k) + 1n];
  }
  ln([a, r]) {
    const lo = a - r;
    if (lo <= 0n) throw new RangeError("logarithm of a ball that is not positive");
    const k = BigInt(a.toString(2).length) - this.p - 1n;
    const m = k >= 0n ? a >> k : a << -k; // in [1, 2), error below one unit when shifted right
    const t = ((m - this.one) << this.p) / (m + this.one); // in [0, 1/3]
    const t2 = (t * t) >> this.p;
    let term = t;
    let sum = t;
    let count = 1n;
    for (let j = 3n; ; j += 2n) {
      term = (term * t2) >> this.p;
      if (term === 0n) break;
      sum += term / j;
      count++;
    }
    // series rounding (3 per term, doubled), the shift of m, and the input radius r/lo
    const base = [2n * sum, 6n * count + 16n];
    const withK = this.add(base, this.scale(this.ln2(), k));
    return [withK[0], withK[1] + (r << this.p) / lo + 1n];
  }
  /** [cos x, sin x]. */
  cosSin([a, r]) {
    const half = this.pi(); // π
    const hp = [half[0] >> 1n, (half[1] >> 1n) + 1n]; // π/2
    const k = (2n * a + hp[0]) / (2n * hp[0]) - (a < 0n && (2n * a + hp[0]) % (2n * hp[0]) !== 0n ? 1n : 0n);
    const t = this.sub([a, 0n], this.scale(hp, k)); // |t| ≤ π/4 + tiny
    const x = t[0];
    const x2 = (x * x) >> this.p;
    let c = this.one;
    let s = x;
    let tc = this.one;
    let ts = x;
    let count = 2n;
    for (let j = 1n; ; j++) {
      tc = -((tc * x2) >> this.p) / ((2n * j - 1n) * (2n * j));
      ts = -((ts * x2) >> this.p) / ((2n * j) * (2n * j + 1n));
      if (tc === 0n && ts === 0n) break;
      c += tc;
      s += ts;
      count++;
    }
    const rad = 3n * count + 4n + t[1] + r; // sin and cos are 1-Lipschitz
    const q = Number(((k % 4n) + 4n) % 4n);
    const C = [c, rad];
    const S = [s, rad];
    return q === 0 ? [C, S] : q === 1 ? [this.neg(S), C] : q === 2 ? [this.neg(C), this.neg(S)] : [S, this.neg(C)];
  }
  // ---- complex balls: [re, im] ----
  cadd([a, b], [c, d]) {
    return [this.add(a, c), this.add(b, d)];
  }
  csub([a, b], [c, d]) {
    return [this.sub(a, c), this.sub(b, d)];
  }
  cmul([a, b], [c, d]) {
    return [this.sub(this.mul(a, c), this.mul(b, d)), this.add(this.mul(a, d), this.mul(b, c))];
  }
  cscale([a, b], n) {
    return [this.scale(a, n), this.scale(b, n)];
  }
  cinv([a, b]) {
    const n = this.inv(this.add(this.mul(a, a), this.mul(b, b)));
    return [this.mul(a, n), this.neg(this.mul(b, n))];
  }
  cdiv(x, y) {
    return this.cmul(x, this.cinv(y));
  }
  cexp([a, b]) {
    const e = this.exp(a);
    const [c, s] = this.cosSin(b);
    return [this.mul(e, c), this.mul(e, s)];
  }
  creal(x) {
    return [x, [0n, 0n]];
  }
  /** A bound for |z|, fixed point. */
  cmag([a, b]) {
    return this.mag(a) + this.mag(b);
  }
}
