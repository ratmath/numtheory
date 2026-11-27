// Dirichlet characters, numbered as in Conrey's scheme (the one the LMFDB
// uses): the character χ_q(m, ·) for each m coprime to q, with
// χ_q(m, n) = χ_q(n, m). A value is a root of unity e(k/N) = exp(2πi·k/N)
// and is kept exactly as the fraction k/N.

import { Rational } from "@ratmath/core";
import { refuse, big, gcd, powmod, mod, inverse } from "./int.js";
import { factor } from "./primes.js";
import { kronecker } from "./quadratic.js";

const MAX_MODULUS = 1000000n;

// ---- logarithm tables for prime powers ----

const tables = new Map();
function table(p, e) {
  const key = `${p}^${e}`;
  if (tables.has(key)) return tables.get(key);
  const pe = p ** e;
  let t;
  if (p === 2) {
    // n ≡ ε·5^b mod 2^e: sign[n] = 0 for ε = 1, 1 for ε = −1; log[n] = b
    const log = new Int32Array(pe);
    const sign = new Uint8Array(pe);
    const half = e >= 3 ? pe / 4 : 1;
    let x = 1;
    for (let b = 0; b < half; b++) {
      log[x] = b;
      sign[x] = 0;
      log[pe - x] = b;
      sign[pe - x] = 1;
      x = (x * 5) % pe;
    }
    if (e === 1) sign[1] = 0;
    t = { p, e, pe, log, sign, phi: pe / 2, half };
  } else {
    // the least g that is a primitive root mod p and mod p²
    const phi = pe - pe / p;
    const P = BigInt(p);
    const qs = factor(P - 1n).map((f) => f.prime);
    let g = 2n;
    for (; ; g++) {
      if (qs.some((q) => powmod(g, (P - 1n) / q, P) === 1n)) continue;
      if (powmod(g, P - 1n, P * P) === 1n) continue;
      break;
    }
    const log = new Int32Array(pe);
    let x = 1;
    const gn = Number(g);
    for (let b = 0; b < phi; b++) {
      log[x] = b;
      x = Number((BigInt(x) * BigInt(gn)) % BigInt(pe));
    }
    t = { p, e, pe, log, phi, g: gn };
  }
  tables.set(key, t);
  return t;
}

function setup(q) {
  q = big(q, "q");
  if (q < 1n) throw new RangeError("the modulus must be at least 1");
  if (q > MAX_MODULUS) refuse("too-large", `characters are tabulated for moduli up to ${MAX_MODULUS}`);
  const Q = Number(q);
  const parts = Q === 1 ? [] : factor(q).map(({ prime, exponent }) => table(Number(prime), exponent));
  return { Q, parts };
}

export class DirichletCharacter {
  constructor(q, m) {
    const { Q, parts } = setup(q);
    const M = Number(mod(big(m, "m"), BigInt(Q)));
    if (Q > 1 && gcd(BigInt(M), BigInt(Q)) !== 1n) refuse("bad-argument", `${m} is not coprime to ${Q}: characters mod ${Q} are numbered by the units mod ${Q}`);
    /** The modulus q. */
    this.modulus = Q;
    /** The Conrey number m: this is χ_q(m, ·). */
    this.number = Q === 1 ? 1 : M;
    this._parts = parts;
    // per prime power: the exponent of n ↦ value is (c·log n)/den (+ sign part for 2)
    this._c = parts.map((t) => {
      const r = M % t.pe;
      if (t.p === 2) return { a: t.log[r], s: t.sign[r] };
      return { a: t.log[r] };
    });
    // the order: lcm of the orders of the parts
    let N = 1n;
    this._c.forEach((c, i) => {
      const t = parts[i];
      let o;
      if (t.p === 2) {
        const oa = t.e >= 3 ? BigInt(t.half) / gcd(BigInt(t.half), BigInt(c.a)) : 1n;
        o = c.s && t.e >= 2 ? (oa > 2n ? oa : 2n) : oa;
      } else o = BigInt(t.phi) / gcd(BigInt(t.phi), BigInt(c.a));
      N = (N * o) / gcd(N, o);
    });
    /** The order of the character in the group of characters. */
    this.order = Number(N);
    Object.freeze(this);
  }
  /** k with χ(n) = e(k/order), or null when gcd(n, q) > 1. */
  _k(n) {
    const Q = this.modulus;
    const r = Number(mod(typeof n === "bigint" ? n : big(n, "n"), BigInt(Q)));
    if (Q === 1) return 0;
    let num = 0n; // exponent as num/den
    let den = 1n;
    for (let i = 0; i < this._parts.length; i++) {
      const t = this._parts[i];
      const x = r % t.pe;
      if (x % t.p === 0) return null;
      const c = this._c[i];
      let a;
      let b;
      if (t.p === 2) {
        // (1−ε_m)(1−ε_n)/8 + a·b/2^(e−2)
        const h = BigInt(t.half);
        a = (c.s && t.sign[x] ? h : 0n) + 2n * BigInt(c.a) * BigInt(t.log[x]);
        b = 2n * h;
        // (1−ε)(1−ε′)/8 = 1/2 when both are −1: as a fraction over 2h it is h
      } else {
        a = BigInt(c.a) * BigInt(t.log[x]);
        b = BigInt(t.phi);
      }
      num = num * b + a * den;
      den *= b;
    }
    const N = BigInt(this.order);
    const k = (num * N) / den;
    return Number(mod(k, N));
  }
  /** The exponent t with χ(n) = e(t) = exp(2πi·t), a fraction in [0, 1); null when χ(n) = 0. */
  exponent(n) {
    const k = this._k(n);
    return k === null ? null : new Rational(BigInt(k), BigInt(this.order));
  }
  /** χ(n) as { re, im } in doubles (exactly 0, ±1 or ±i when it is one of those). */
  value(n) {
    const k = this._k(n);
    if (k === null) return { re: 0, im: 0 };
    return unitValue(k, this.order);
  }
  /** χ(n) for a real character: 0n, 1n or −1n. */
  at(n) {
    if (this.order > 2) refuse("bad-argument", "at(n) is for real characters; use value(n) or exponent(n)");
    const k = this._k(n);
    return k === null ? 0n : k === 0 ? 1n : -1n;
  }
  /** The values χ(0), χ(1), …, χ(q − 1) as exponents (null for 0). */
  exponents() {
    return Array.from({ length: this.modulus }, (_, n) => this.exponent(n));
  }
  get isPrincipal() {
    return this.order === 1;
  }
  get isReal() {
    return this.order <= 2;
  }
  /** χ(−1) = 1. */
  get isEven() {
    return this._k(-1) === 0;
  }
  get parity() {
    return this.isEven ? "even" : "odd";
  }
  /** The conductor: the smallest modulus the character can be defined with. */
  get conductor() {
    let f = 1;
    this._parts.forEach((t, i) => {
      const c = this._c[i];
      if (t.p === 2) {
        if (t.e >= 3 && c.a % t.half !== 0) {
          let v = 0;
          let a = c.a;
          while (a % 2 === 0) {
            a /= 2;
            v++;
          }
          f *= 2 ** (t.e - v);
        } else if (c.s && t.e >= 2) f *= 4;
      } else if (c.a % t.phi !== 0) {
        let v = 0;
        let a = c.a;
        while (a % t.p === 0) {
          a /= t.p;
          v++;
        }
        f *= t.p ** (t.e - v);
      }
    });
    return f;
  }
  get isPrimitive() {
    return this.conductor === this.modulus;
  }
  /** The primitive character, modulo the conductor, that induces this one. */
  primitive() {
    const f = this.conductor;
    if (f === this.modulus) return this;
    if (f === 1) return new DirichletCharacter(1, 1);
    // match the values on the units mod f: search the characters of the same order
    const G = dirichletGroup(f);
    const probe = [];
    for (let n = 2; probe.length < 40 && n < 4 * f + 50; n++) if (gcd(BigInt(n), BigInt(this.modulus)) === 1n) probe.push(n);
    const want = probe.map((n) => this._k(n));
    for (const chi of G.characters) {
      if (chi.order !== this.order || chi.conductor !== f) continue;
      if (probe.every((n, i) => chi._k(n) === want[i]) && sameOnUnits(this, chi)) return chi;
    }
    throw new Error("no inducing character found");
  }
  /** The product χ·ψ, modulo the least common multiple of the moduli. */
  mul(psi) {
    const q = BigInt(this.modulus);
    const r = BigInt(psi.modulus);
    if (q === r) return new DirichletCharacter(q, (BigInt(this.number) * BigInt(psi.number)) % q);
    const L = (q * r) / gcd(q, r);
    return induced(this, L).mul(induced(psi, L));
  }
  pow(k) {
    return new DirichletCharacter(this.modulus, powmodSigned(BigInt(this.number), BigInt(k), BigInt(this.modulus)));
  }
  /** The complex conjugate character, also the inverse. */
  conjugate() {
    return this.pow(-1);
  }
  /** The Gauss sum τ(χ) = Σ χ(a)·e(a/q) in doubles; |τ|² = q for a primitive character. */
  gaussSum() {
    const q = this.modulus;
    let re = 0;
    let im = 0;
    for (let a = 1; a <= q; a++) {
      const k = this._k(a);
      if (k === null) continue;
      const ang = 2 * Math.PI * (k / this.order + a / q);
      re += Math.cos(ang);
      im += Math.sin(ang);
    }
    return { re, im };
  }
  /** The same character seen modulo a multiple of q (zero on the new common factors). */
  induce(modulus) {
    return induced(this, big(modulus, "modulus"));
  }
  equals(psi) {
    return psi instanceof DirichletCharacter && psi.modulus === this.modulus && psi.number === this.number;
  }
  toString() {
    return `χ_${this.modulus}(${this.number}, ·)`;
  }
}
function unitValue(k, N) {
  const g = Number(gcd(BigInt(k), BigInt(N)));
  const [a, b] = [k / g, N / g];
  if (b === 1) return { re: 1, im: 0 };
  if (b === 2) return { re: -1, im: 0 };
  if (b === 4) return { re: 0, im: a === 1 ? 1 : -1 };
  const ang = (2 * Math.PI * a) / b;
  return { re: Math.cos(ang), im: Math.sin(ang) };
}
function powmodSigned(a, k, q) {
  if (q === 1n) return 0n;
  return k >= 0n ? powmod(a, k, q) : powmod(inverse(a, q), -k, q);
}
function sameOnUnits(chi, psi) {
  // χ mod q and ψ mod f | q agree on every n coprime to q; q is small enough to check all when below 5000
  const q = chi.modulus;
  if (q > 5000) return true;
  for (let n = 1; n < q; n++) {
    const k = chi._k(n);
    if (k !== null && k !== psi._k(n)) return false;
  }
  return true;
}
function induced(chi, L) {
  const q = BigInt(chi.modulus);
  if (L % q !== 0n) refuse("bad-argument", `${L} is not a multiple of the modulus ${q}`);
  if (L === q) return chi;
  // the character mod L with the same values on units: found among those of the same order
  const G = dirichletGroup(L);
  const probe = [];
  for (let n = 2; probe.length < 40 && n < 4 * Number(L) + 50; n++) if (gcd(BigInt(n), L) === 1n) probe.push(n);
  const want = probe.map((n) => chi._k(n));
  for (const psi of G.characters) {
    if (psi.order !== chi.order) continue;
    if (probe.every((n, i) => psi._k(n) === want[i]) && sameOnUnits(psi, chi)) return psi;
  }
  throw new Error("no induced character found");
}

/** The character χ_q(m, ·) in Conrey's numbering; m must be coprime to q. */
export function dirichletCharacter(q, m = 1) {
  return new DirichletCharacter(q, m);
}

const groups = new Map();
/**
 * The group of Dirichlet characters modulo q: { modulus, order, characters,
 * principal, primitive, real, structure }. The group is isomorphic to the
 * units mod q; structure lists the orders of cyclic factors (one for each
 * odd prime power, up to two for the power of 2).
 */
export function dirichletGroup(q) {
  q = big(q, "q");
  const key = q.toString();
  if (groups.has(key)) return groups.get(key);
  const { Q, parts } = setup(q);
  if (Q > 20000) refuse("too-large", "listing every character stops at modulus 20000; build single ones with dirichletCharacter");
  const characters = [];
  for (let m = 1; m <= Math.max(1, Q - (Q > 1 ? 1 : 0)); m++) if (Q === 1 || gcd(BigInt(m), q) === 1n) characters.push(new DirichletCharacter(q, m));
  const structure = [];
  for (const t of parts) {
    if (t.p === 2) {
      if (t.e === 2) structure.push(2);
      if (t.e >= 3) structure.push(2, t.half);
    } else structure.push(t.phi);
  }
  const G = Object.freeze({
    modulus: Q,
    order: characters.length,
    characters,
    principal: characters[0],
    get primitive() {
      return characters.filter((c) => c.isPrimitive);
    },
    get real() {
      return characters.filter((c) => c.isReal);
    },
    structure,
    /** Σ_χ χ(a)·conj(χ(b)): φ(q) when a ≡ b are units mod q, else 0. */
    orthogonality(a, b) {
      let re = 0;
      let im = 0;
      for (const c of characters) {
        const ka = c._k(a);
        const kb = c._k(b);
        if (ka === null || kb === null) continue;
        const v = unitValue((((ka - kb) % c.order) + c.order) % c.order, c.order);
        re += v.re;
        im += v.im;
      }
      return { re: Math.round(re * 1e9) / 1e9, im: Math.round(im * 1e9) / 1e9 };
    },
    toString() {
      return `the ${characters.length} Dirichlet characters mod ${Q}`;
    },
  });
  groups.set(key, G);
  return G;
}

/**
 * The Kronecker symbol (D/·) as a Dirichlet character modulo |D|, for a
 * discriminant D ≡ 0 or 1 mod 4 that is not a square. For a fundamental
 * discriminant it is the primitive quadratic character of Q(√D).
 */
export function kroneckerCharacter(D) {
  D = big(D, "D");
  if (mod(D, 4n) > 1n) refuse("bad-argument", `${D} is not a discriminant: it must be 0 or 1 mod 4`);
  const q = D < 0n ? -D : D;
  setup(q);
  // a real character has a number m with m² ≡ 1; match the symbol on small n
  const Q = Number(q);
  const probe = [];
  for (let n = 2; probe.length < 60 && n < 4 * Q + 60; n++) if (gcd(BigInt(n), q) === 1n) probe.push(n);
  const want = probe.map((n) => BigInt(kronecker(D, BigInt(n))));
  for (let m = 1; m < Math.max(2, Q); m++) {
    if ((m * m) % Q !== 1 % Q) continue;
    const chi = new DirichletCharacter(q, m);
    if (probe.every((n, i) => chi.at(n) === want[i])) return chi;
  }
  throw new Error("the Kronecker symbol was not found among the real characters");
}
