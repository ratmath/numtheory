// Binary quadratic forms f(x, y) = a·x² + b·x·y + c·y² with integer
// coefficients, discriminant D = b² − 4ac.
//
// Two forms are (properly) equivalent when one becomes the other under an
// integer change of variables of determinant 1. Every class contains reduced
// forms: exactly one when D < 0 (Gauss and Lagrange), a cycle of them when
// D > 0. Equivalence is decided by reducing, and every equivalence comes with
// the matrix that proves it. The classes of primitive forms of one
// discriminant form a finite abelian group under Gauss's composition, the
// form class group; its order is the class number h(D).

import { big, egcd, fdiv, gcd, isqrt, isSquare, mod, refuse } from "./int.js";
import { sqrtsMod } from "./quadratic.js";
import { divisors } from "./arith.js";
import { factor } from "./primes.js";

const I = [[1n, 0n], [0n, 1n]];
const mul = ([[a, b], [c, d]], [[e, f], [g, h]]) => [[a * e + b * g, a * f + b * h], [c * e + d * g, c * f + d * h]];
const inv = ([[p, q], [r, s]]) => [[s, -q], [-r, p]];

export class Form {
  constructor(a, b, c) {
    this.a = big(a, "a");
    this.b = big(b, "b");
    this.c = big(c, "c");
    Object.freeze(this);
  }
  get discriminant() {
    return this.b * this.b - 4n * this.a * this.c;
  }
  /** f(x, y). */
  at(x, y) {
    x = big(x, "x");
    y = big(y, "y");
    return this.a * x * x + this.b * x * y + this.c * y * y;
  }
  /** f ∘ M: the form f(px + qy, rx + sy) for M = [[p, q], [r, s]]. */
  transform([[p, q], [r, s]]) {
    const { a, b, c } = this;
    return new Form(a * p * p + b * p * r + c * r * r, 2n * a * p * q + b * (p * s + q * r) + 2n * c * r * s, a * q * q + b * q * s + c * s * s);
  }
  isPrimitive() {
    return gcd(gcd(this.a, this.b), this.c) === 1n;
  }
  /** Positive definite: D < 0 and a > 0. */
  isPositiveDefinite() {
    return this.discriminant < 0n && this.a > 0n;
  }
  equals(g) {
    return g instanceof Form && g.a === this.a && g.b === this.b && g.c === this.c;
  }
  isReduced() {
    const D = this.discriminant;
    const { a, b, c } = this;
    if (D < 0n) return a > 0n && abs(b) <= a && a <= c && (b >= 0n || (abs(b) !== a && a !== c));
    if (isSquare(D)) return false;
    const s = isqrt(D);
    const A = abs(a);
    return b >= 1n && b <= s && 2n * A + b >= s + 1n && 2n * A - b <= s;
  }
  /** { form, matrix } with form reduced and this.transform(matrix) equal to it. */
  reduce() {
    return reduceForm(this);
  }
  /** The cycle of reduced forms of an indefinite form's class, with the matrix taking each to the next. */
  cycle() {
    if (this.discriminant <= 0n) refuse("bad-argument", "only indefinite forms (D > 0) have cycles");
    return cycleOf(reduceForm(this).form).map(({ form }) => form);
  }
  toString() {
    const t = (k, v, first) => {
      if (k === 0n) return "";
      const sg = k < 0n ? (first ? "−" : " − ") : first ? "" : " + ";
      const m = abs(k) === 1n && v ? "" : String(abs(k));
      return `${sg}${m}${v}`;
    };
    const parts = [t(this.a, "x²", true)];
    parts.push(t(this.b, "xy", !parts.join("")));
    parts.push(t(this.c, "y²", !parts.join("")));
    return parts.join("") || "0";
  }
  toArray() {
    return [this.a, this.b, this.c];
  }
}
const abs = (x) => (x < 0n ? -x : x);

/** A form from (a, b, c) or [a, b, c]. */
export function form(a, b, c) {
  if (Array.isArray(a)) return new Form(...a);
  if (a instanceof Form) return a;
  return new Form(a, b, c);
}
const asForm = (f) => (f instanceof Form ? f : form(f));

// ---- reduction ----

function reduceDefinite(f0) {
  let f = f0;
  let M = I;
  const apply = (N) => {
    f = f.transform(N);
    M = mul(M, N);
  };
  for (let guard = 0; guard < 100000; guard++) {
    if (abs(f.b) > f.a || f.b === -f.a) {
      // translate b into (−a, a]: b + 2ak
      const k = fdiv(f.a - f.b, 2n * f.a);
      apply([[1n, k], [0n, 1n]]);
      continue;
    }
    if (f.a > f.c || (f.a === f.c && f.b < 0n)) {
      apply([[0n, -1n], [1n, 0n]]);
      continue;
    }
    return { form: f, matrix: M };
  }
  throw new Error("reduction did not finish");
}

/** One step of the reduction operator ρ for D > 0: (a, b, c) ↦ (c, b', ·) with b' ≡ −b (mod 2c) normalised. */
function rho(f, s) {
  const c = f.c;
  const C = abs(c);
  // b' ≡ −b (mod 2|c|): in [s + 1 − 2|c|, s] when |c| < √D, else in (−|c|, |c|]
  let lo = C <= s ? s + 1n - 2n * C : -C + 1n;
  const bp = lo + mod(-f.b - lo, 2n * C);
  const t = (bp + f.b) / (2n * c);
  const N = [[0n, -1n], [1n, t]];
  return { form: f.transform(N), N };
}
function reduceIndefinite(f0) {
  const D = f0.discriminant;
  if (isSquare(D)) refuse("bad-argument", `the discriminant ${D} is a square: the form factors into linear forms`);
  const s = isqrt(D);
  let f = f0;
  let M = I;
  if (f.c === 0n || f.a === 0n) refuse("bad-argument", "a form with a = 0 or c = 0 has a square discriminant");
  for (let guard = 0; guard < 1000000; guard++) {
    if (f.isReduced()) return { form: f, matrix: M };
    const { form: g, N } = rho(f, s);
    f = g;
    M = mul(M, N);
  }
  throw new Error("reduction did not finish");
}
function reduceForm(f) {
  const D = f.discriminant;
  if (D === 0n) refuse("bad-argument", "the discriminant is 0");
  if (D < 0n) {
    if (f.a < 0n) refuse("bad-argument", "a negative definite form: reduce −f instead");
    return reduceDefinite(f);
  }
  return reduceIndefinite(f);
}
/** The ρ-cycle of a reduced indefinite form: [{ form, matrix }] where matrix takes the first form to this one. */
function cycleOf(f) {
  const s = isqrt(f.discriminant);
  const out = [{ form: f, matrix: I }];
  let g = f;
  let M = I;
  for (let guard = 0; guard < 1000000; guard++) {
    const { form: h, N } = rho(g, s);
    M = mul(M, N);
    if (h.equals(f)) return out;
    out.push({ form: h, matrix: M });
    g = h;
  }
  throw new Error("cycle too long");
}

/**
 * Proper equivalence: { equivalent: true, matrix } with f.transform(matrix)
 * equal to g (determinant 1), or { equivalent: false }.
 */
export function equivalent(f, g) {
  f = asForm(f);
  g = asForm(g);
  if (f.discriminant !== g.discriminant) return { equivalent: false, reason: "different discriminants" };
  const rf = reduceForm(f);
  const rg = reduceForm(g);
  if (f.discriminant < 0n) {
    if (!rf.form.equals(rg.form)) return { equivalent: false, reason: `they reduce to different forms, ${rf.form} and ${rg.form}` };
    return { equivalent: true, matrix: mul(rf.matrix, inv(rg.matrix)) };
  }
  const cyc = cycleOf(rf.form);
  const hit = cyc.find((x) => x.form.equals(rg.form));
  if (!hit) return { equivalent: false, reason: "the reduced form of the second is not in the cycle of the first" };
  return { equivalent: true, matrix: mul(mul(rf.matrix, hit.matrix), inv(rg.matrix)) };
}

// ---- composition and the class group ----

/** Gauss composition of two primitive forms of the same discriminant (Shanks's formulas, as in Cohen's Algorithm 5.4.7), reduced. */
export function compose(f1, f2) {
  f1 = asForm(f1);
  f2 = asForm(f2);
  const D = f1.discriminant;
  if (f2.discriminant !== D) refuse("bad-argument", "the forms have different discriminants");
  f1 = positiveA(f1);
  f2 = positiveA(f2);
  if (f1.a > f2.a) [f1, f2] = [f2, f1];
  const { a: a1, b: b1 } = f1;
  const { a: a2, b: b2, c: c2 } = f2;
  const s = (b1 + b2) / 2n;
  const n = b2 - s;
  let y1;
  let d;
  if (a2 % a1 === 0n) {
    y1 = 0n;
    d = a1;
  } else {
    const e = egcd(a2, a1);
    y1 = e.x;
    d = e.g;
  }
  let x2;
  let y2;
  let d1;
  if (s % d === 0n) {
    y2 = -1n;
    x2 = 0n;
    d1 = d;
  } else {
    const e = egcd(s, d);
    x2 = e.x;
    y2 = -e.y;
    d1 = e.g;
  }
  const v1 = a1 / d1;
  const v2 = a2 / d1;
  const r = mod(y1 * y2 * n - x2 * c2, v1);
  const b3 = b2 + 2n * v2 * r;
  const a3 = v1 * v2;
  const c3 = (b3 * b3 - D) / (4n * a3);
  return reduceForm(new Form(a3, b3, c3)).form;
}
/** An equivalent form with a > 0 (for D > 0, a neighbour in the cycle). */
function positiveA(f) {
  if (f.a > 0n) return f;
  if (f.discriminant < 0n) refuse("bad-argument", "negative definite forms are not in the class group");
  const r = reduceForm(f).form;
  const c = cycleOf(r).find((x) => x.form.a > 0n);
  if (!c) refuse("bad-argument", "no form with a > 0 in the cycle");
  return c.form;
}
/** The principal form of discriminant D: x² + (D mod 2)xy + ((D mod 2) − D)/4·y², reduced. */
export function principalForm(D) {
  const d = big(D, "D");
  if (mod(d, 4n) > 1n) refuse("bad-argument", `${d} is not a discriminant: it must be 0 or 1 mod 4`);
  const b = mod(d, 2n);
  return reduceForm(new Form(1n, b, (b * b - d) / 4n)).form;
}
/** The inverse class: (a, −b, c), reduced. */
export function inverseForm(f) {
  f = asForm(f);
  return reduceForm(new Form(f.a, -f.b, f.c)).form;
}

/** The reduced primitive forms of discriminant D (positive definite ones when D < 0), in a fixed order. */
export function reducedForms(D) {
  const d = big(D, "D");
  if (mod(d, 4n) > 1n) refuse("bad-argument", `${d} is not a discriminant: it must be 0 or 1 mod 4`);
  if (d === 0n || (d > 0n && isSquare(d))) refuse("bad-argument", "the discriminant must not be a square");
  const out = [];
  if (d < 0n) {
    const A = isqrt(-d / 3n);
    if (A > 10000000n) refuse("too-large", "the discriminant is too large to list its forms");
    for (let a = 1n; a <= A; a++) {
      for (let b = -a + 1n; b <= a; b++) {
        if (mod(b - d, 2n) !== 0n) continue;
        const num = b * b - d;
        if (num % (4n * a) !== 0n) continue;
        const c = num / (4n * a);
        const f = new Form(a, b, c);
        if (f.isReduced() && f.isPrimitive()) out.push(f);
      }
    }
    return out;
  }
  const s = isqrt(d);
  if (s > 10000000n) refuse("too-large", "the discriminant is too large to list its forms");
  for (let b = 1n; b <= s; b++) {
    if (mod(b - d, 2n) !== 0n) continue;
    const ac = (b * b - d) / 4n; // = a·c < 0
    for (const m of divisors(-ac)) {
      for (const a of [m, -m]) {
        const c = ac / a;
        const f = new Form(a, b, c);
        if (f.isReduced() && f.isPrimitive()) out.push(f);
      }
    }
  }
  return out.sort((x, y) => (x.b !== y.b ? (x.b < y.b ? -1 : 1) : x.a < y.a ? -1 : x.a > y.a ? 1 : 0));
}

/**
 * The form class group of discriminant D (narrow, when D > 0): classes of
 * primitive forms (positive definite ones when D < 0) under composition.
 * Returns { D, h, forms, identity, classOf(f), compose(f, g), inverse(f),
 * order(f), structure } where forms has one reduced form per class and
 * structure lists the invariant factors (the group is their direct product of
 * cyclic groups).
 */
export function classGroup(D) {
  const d = big(D, "D");
  let reps;
  let cycles = null;
  if (d < 0n) reps = reducedForms(d);
  else {
    const all = reducedForms(d);
    const seen = new Set();
    reps = [];
    cycles = [];
    for (const f of all) {
      const key = f.toArray().join(",");
      if (seen.has(key)) continue;
      const cyc = cycleOf(f).map((x) => x.form);
      for (const g of cyc) seen.add(g.toArray().join(","));
      const rep = cyc.find((g) => g.a > 0n) ?? f;
      reps.push(rep);
      cycles.push(cyc);
    }
  }
  const keyOf = (f) => {
    const r = reduceForm(asForm(f)).form;
    if (d < 0n) return reps.findIndex((g) => g.equals(r));
    return cycles.findIndex((c) => c.some((g) => g.equals(r)));
  };
  const classOf = (f) => {
    const k = keyOf(f);
    if (k < 0) refuse("bad-argument", "not a primitive form of this discriminant");
    return reps[k];
  };
  const identity = classOf(principalForm(d));
  const op = (f, g) => classOf(compose(f, g));
  const order = (f) => {
    let g = classOf(f);
    for (let k = 1n; ; k++) {
      if (g.equals(identity)) return k;
      g = op(g, f);
      if (k > BigInt(reps.length)) throw new Error("order exceeds the class number");
    }
  };
  // invariant factors from counts of elements killed by p^j
  const h = BigInt(reps.length);
  const structure = [];
  if (reps.length <= 5000) {
    const orders = reps.map(order);
    const parts = [];
    for (const { prime: p } of h > 1n ? factor(h, { prove: false }) : []) {
      const counts = [1n];
      for (let j = 1; ; j++) {
        const pj = p ** BigInt(j);
        counts.push(BigInt(orders.filter((o) => pj % o === 0n).length));
        if (counts[j] === counts[j - 1]) break;
      }
      // number of cyclic factors of order ≥ p^j = log_p(counts[j]/counts[j−1])
      const atLeast = [];
      for (let j = 1; j < counts.length; j++) {
        let q = counts[j] / counts[j - 1];
        let k = 0;
        while (q > 1n) {
          q /= p;
          k++;
        }
        atLeast.push(k);
      }
      const exps = [];
      for (let j = 0; j < atLeast.length; j++) {
        const exactly = atLeast[j] - (atLeast[j + 1] ?? 0);
        for (let i = 0; i < exactly; i++) exps.push(p ** BigInt(j + 1));
      }
      parts.push(exps.sort((x, y) => (x > y ? -1 : 1)));
    }
    const width = Math.max(0, ...parts.map((x) => x.length));
    for (let i = 0; i < width; i++) structure.push(parts.reduce((acc, ps) => acc * (ps[i] ?? 1n), 1n));
    structure.reverse();
  }
  return {
    D: d,
    h,
    forms: reps,
    identity,
    classOf,
    compose: op,
    inverse: (f) => classOf(inverseForm(f)),
    order,
    structure: reps.length <= 5000 ? (structure.length ? structure : [1n].slice(0, h === 1n ? 0 : 1)) : null,
  };
}
/** The class number h(D) (the narrow one when D > 0). */
export function classNumber(D) {
  return classGroup(D).h;
}

// ---- representations ----

/**
 * Every (x, y) with f(x, y) = n, for a positive definite f (finitely many):
 * from 4a·f = (2ax + by)² + |D|·y², so |y| ≤ √(4an/|D|).
 */
export function representations(f, n) {
  f = asForm(f);
  const N = big(n, "n");
  if (!f.isPositiveDefinite()) refuse("bad-argument", "representations lists them for positive definite forms; use represent for others");
  const D = -f.discriminant;
  const Y = isqrt((4n * f.a * N) / D);
  if (Y > 10000000n) refuse("too-large", "n is too large to search");
  const out = [];
  for (let y = -Y; y <= Y; y++) {
    // a x² + b y x + (c y² − n) = 0
    const disc = (f.b * y) ** 2n - 4n * f.a * (f.c * y * y - N);
    if (disc < 0n || !isSquare(disc)) continue;
    const r = isqrt(disc);
    for (const sg of r ? [1n, -1n] : [1n]) {
      const num = -f.b * y + sg * r;
      if (num % (2n * f.a) === 0n) out.push([num / (2n * f.a), y]);
    }
  }
  return out.sort((u, v) => (u[0] !== v[0] ? (u[0] < v[0] ? -1 : 1) : u[1] < v[1] ? -1 : 1));
}

/**
 * One proper representation (gcd(x, y) = 1) of n > 0 by f, or null. Found
 * by the classical criterion: f properly represents n exactly when f is
 * equivalent to a form (n, b, ·) with b² ≡ D (mod 4n), and then the first
 * column of the equivalence matrix is a representation. Works for any D.
 */
export function represent(f, n) {
  f = asForm(f);
  const N = big(n, "n");
  if (N <= 0n) refuse("bad-argument", "n must be positive (use −f for negative values)");
  const D = f.discriminant;
  for (const b of sqrtsMod(D, 4n * N)) {
    if (b >= 2n * N) continue; // b and b + 2n give equivalent forms
    const g = new Form(N, b, (b * b - D) / (4n * N));
    const e = equivalent(f, g);
    if (e.equivalent) {
      const [x, y] = [e.matrix[0][0], e.matrix[1][0]];
      if (f.at(x, y) !== N) throw new Error("internal check failed");
      return [x, y];
    }
  }
  return null;
}
