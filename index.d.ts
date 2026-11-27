// Type declarations for @ratmath/numtheory.

import type { Rational, RationalInterval } from "@ratmath/core";
import type { AlgebraicReal, NumberField, NumberFieldElement } from "@ratmath/algebraic";

/** A whole number: a JavaScript integer, a BigInt, decimal text, or a Rational with denominator 1. Non-whole JS numbers are refused. */
export type Whole = bigint | number | string | Rational;
/** An exact rational: a whole number, fraction text ("3/7", "0.25") or a Rational. */
export type RationalLike = bigint | number | string | Rational;

/** Thrown for refusals; `code` says why ("bad-argument", "too-large", "unproven", "unknown-valuation", "not-a-square", "diverges", "too-hard", "singular", "not-on-curve", "no-inverse", "budget", "pole", "precision", "unsupported"). */
export class NumberTheoryError extends Error {
  code: string;
}

// ---- integers

export function isqrt(n: bigint): bigint;
export function iroot(n: bigint, k: number | bigint): bigint;
export function gcd(a: bigint, b: bigint): bigint;
export function lcm(a: bigint, b: bigint): bigint;
export function egcd(a: bigint, b: bigint): { g: bigint; x: bigint; y: bigint };

// ---- primes, factoring, certificates

/** Below this bound isPrime is a proof (Baillie–PSW has no exceptions there). */
export const PROVEN_BELOW: bigint;
export function primes(limit: Whole): bigint[];
export function primePi(x: Whole): bigint;
export function nthPrime(k: Whole): bigint;
/** Proven below 2^64; above, by a certificate, or throws with code "unproven". */
export function isPrime(n: Whole): boolean;
/** Baillie–PSW: fast, not a proof above 2^64. */
export function isProbablePrime(n: bigint | number): boolean;
export function nextPrime(n: Whole): bigint;
export function prevPrime(n: Whole): bigint;
export function primesBetween(a: Whole, b: Whole): bigint[];
export function pollardRho(n: Whole, options?: { steps?: number; seed?: bigint }): bigint | null;
export function pollardPm1(n: Whole, bound?: number): bigint | null;
export function trialDivide(n: Whole, limit?: number): { found: Map<bigint, number>; rest: bigint };
export interface PrimeFactor {
  prime: bigint;
  exponent: number;
  /** True below 2^64 and for larger primes with a certificate. */
  proven: boolean;
}
/** Trial division, rho, p − 1, then the elliptic curve method; `effort` scales the number of curves. */
export function factor(n: Whole, options?: { prove?: boolean; effort?: number }): PrimeFactor[];
export function formatFactors(factors: Array<{ prime: bigint; exponent: number }>): string;

export interface SmallCertificate {
  n: bigint;
  kind: "small";
}
export interface CertifiedFactor {
  q: bigint;
  e: number;
  /** The base a with a^(n−1) ≡ 1 and gcd(a^((n−1)/q) − 1, n) = 1. */
  a: bigint;
  certificate: Certificate;
}
export interface PocklingtonCertificate {
  n: bigint;
  kind: "pocklington";
  /** n − 1 = F·R with F fully factored. */
  F: bigint;
  R: bigint;
  factors: CertifiedFactor[];
}
export interface BLSCertificate {
  n: bigint;
  kind: "bls";
  F: bigint;
  R: bigint;
  c1: bigint;
  c2: bigint;
  factors: CertifiedFactor[];
}
/** A curve y² = x³ + ax + b mod n and a point (x, y) with m·P = 0 and (m/q)·P ≠ 0, for a certified prime q > (n^(1/4) + 1)². */
export interface EllipticCertificate {
  n: bigint;
  kind: "elliptic";
  /** The complex multiplication discriminant the curve was built from (information only; not checked). */
  discriminant?: number;
  a: bigint;
  b: bigint;
  x: bigint;
  y: bigint;
  m: bigint;
  q: bigint;
  certificate: Certificate;
}
export type Certificate = SmallCertificate | PocklingtonCertificate | BLSCertificate | EllipticCertificate;
export type PrimalityProof =
  | { n: bigint; prime: true; proven: true; certificate: Certificate }
  | { n: bigint; prime: false; proven: true; factor?: bigint; witness?: bigint; reason?: string }
  | { n: bigint; prime: boolean; proven: false; reason: string };
export function provePrime(
  n: Whole,
  options?: { budget?: number; elliptic?: boolean; maxDiscriminant?: number; maxClassNumber?: number },
): PrimalityProof;
export function millerRabinWitness(n: bigint, a: bigint): boolean;
/** The independent checker: shares no code with the search. */
export function verifyCertificate(cert: Certificate): { valid: true } | { valid: false; reason: string };
export function verifyComposite(n: bigint, proof: { factor?: bigint; witness?: bigint }): { valid: true } | { valid: false; reason: string };

// ---- elliptic curves

/** A point of a curve: { x, y }, or the zero of the group { infinity: true }. T is Rational over Q and bigint over F_p. */
export type CurvePoint<T> = { readonly x: T; readonly y: T; readonly infinity?: undefined } | { readonly infinity: true };
export type PointLike<T> = CurvePoint<T> | [RationalLike, RationalLike] | { x: RationalLike; y: RationalLike } | null;
export interface Reduction {
  prime: bigint;
  type: "good" | "multiplicative" | "additive";
  /** For multiplicative reduction: are the tangents at the node defined over F_p? */
  split?: boolean;
  discriminantValuation: number;
  /** 0, 1, or 2 for additive reduction at p ≥ 5; null for additive reduction at 2 or 3. */
  conductorExponent: number | null;
}
export interface ModelChange<T> {
  curve: EllipticCurve<Rational>;
  to(P: PointLike<T>): CurvePoint<Rational>;
  from(P: PointLike<Rational>): CurvePoint<T>;
}
/** y² + a1·xy + a3·y = x³ + a2·x² + a4·x + a6 over Q (T = Rational) or F_p (T = bigint). */
export class EllipticCurve<T = Rational | bigint> {
  readonly a1: T;
  readonly a2: T;
  readonly a3: T;
  readonly a4: T;
  readonly a6: T;
  readonly b2: T;
  readonly b4: T;
  readonly b6: T;
  readonly b8: T;
  readonly c4: T;
  readonly c6: T;
  readonly discriminant: T;
  readonly j: T;
  /** null over Q, the prime p over F_p. */
  readonly modulus: bigint | null;
  readonly field: string;
  readonly coefficients: [T, T, T, T, T];
  readonly isShort: boolean;
  readonly infinity: { readonly infinity: true };
  toString(): string;
  point(x: RationalLike, y: RationalLike): CurvePoint<T>;
  has(P: PointLike<T>): boolean;
  equal(P: PointLike<T>, Q: PointLike<T>): boolean;
  neg(P: PointLike<T>): CurvePoint<T>;
  add(P: PointLike<T>, Q: PointLike<T>): CurvePoint<T>;
  sub(P: PointLike<T>, Q: PointLike<T>): CurvePoint<T>;
  double(P: PointLike<T>): CurvePoint<T>;
  mul(k: Whole, P: PointLike<T>): CurvePoint<T>;
  /** A bigint; over Q, Infinity for a point of infinite order. */
  order(P: PointLike<T>): bigint | number;

  // over Q
  readonly isIntegral: boolean;
  shortModel(): ModelChange<T> & { A: bigint; B: bigint };
  minimalModel(): ModelChange<T> & { u: Rational; r: Rational; s: Rational; t: Rational };
  badPrimes(): bigint[];
  reduction(p: Whole): Reduction;
  /** Refuses with "too-hard" when the reduction at 2 or 3 is additive. */
  conductor(): bigint;
  modulo(p: Whole): EllipticCurve<bigint>;
  torsion(): { order: number; structure: number[]; points: CurvePoint<T>[]; generators: CurvePoint<T>[] };
  naiveHeight(P: PointLike<T>): bigint;
  /** An enclosure of ĥ(P) = lim log H(x(2ⁿP))/4ⁿ; exact (0) for torsion points. */
  height(P: PointLike<T>, options?: { doublings?: number }): { lower: Rational; upper: Rational; exact: boolean };
  rationalPoints(bound?: Whole): CurvePoint<T>[];

  // over F_p
  count(options?: { method?: "auto" | "naive" | "bsgs" | "schoof"; log?: Array<{ l: number; t: number }> }): bigint;
  trace(): bigint;
  liftX(x: RationalLike): CurvePoint<T>[];
  points(): CurvePoint<T>[];
  randomPoint(seed?: Whole): CurvePoint<T>;
  structure(): { order: bigint; structure: bigint[] };
  divisionPolynomial(n: number): bigint[];
}
/** [a, b] for y² = x³ + ax + b, or [a1, a2, a3, a4, a6]; over Q, or over F_p with option mod. */
export function ellipticCurve(coefficients: RationalLike[], options?: { mod?: undefined | null }): EllipticCurve<Rational>;
export function ellipticCurve(coefficients: RationalLike[], options: { mod: Whole }): EllipticCurve<bigint>;

export interface EcmReport {
  factor: bigint | null;
  /** 0: found while setting up the curve; 1 or 2: the stage; null: not found. */
  stage: 0 | 1 | 2 | null;
  sigma: bigint | null;
  B1: number;
  B2: number;
}
/** One Suyama curve (σ ≥ 6) of the elliptic curve method. */
export function ecmCurve(n: Whole, sigma: Whole, options?: { B1?: number; B2?: number }): EcmReport;
/** The elliptic curve method: a proper factor, or null. */
export function ecm(n: Whole, options?: { B1?: number; B2?: number; curves?: number; sigma?: Whole; details?: false }): bigint | null;
export function ecm(n: Whole, options: { B1?: number; B2?: number; curves?: number; sigma?: Whole; details: true }): EcmReport | null;
/** The Hilbert class polynomial H_D (constant term first), or null when its rounding check fails. */
export function classPolynomial(D: number | bigint): bigint[] | null;

// ---- arithmetic functions

export function divisors(n: Whole): bigint[];
export function sigma(n: Whole, k?: number): bigint;
export function tau(n: Whole): bigint;
export function mobius(n: Whole): bigint;
export function totient(n: Whole): bigint;
export function jordan(n: Whole, k: number): bigint;
export function carmichael(n: Whole): bigint;
export function omega(n: Whole): bigint;
export function bigOmega(n: Whole): bigint;
export function liouville(n: Whole): bigint;
export function radical(n: Whole): bigint;
export function isSquarefree(n: Whole): boolean;
export function aliquotSum(n: Whole): bigint;
export function classify(n: Whole): "perfect" | "abundant" | "deficient";
export function isPerfect(n: Whole): boolean;
export function aliquotSequence(n: Whole, options?: { steps?: number; max?: bigint }): { terms: bigint[]; ending: string };
export function amicablePairs(limit: Whole): Array<[bigint, bigint]>;
export function perfectNumbers(count: Whole): bigint[];
export function lucasLehmer(p: Whole): boolean;
export type ArithmeticFunction = (n: bigint) => bigint;
export function dirichlet(f: ArithmeticFunction, g: ArithmeticFunction): ArithmeticFunction;
export function mobiusInversion(F: ArithmeticFunction): ArithmeticFunction;
export const unit: ArithmeticFunction;
export const one: ArithmeticFunction;
export const id: ArithmeticFunction;
export function isMultiplicative(f: ArithmeticFunction, limit?: number): boolean;

// ---- congruences

export type IntPoly = Whole[];
export function solveLinear(
  a: Whole,
  b: Whole,
  m: Whole,
): { solutions: bigint[]; modulus: bigint; g: bigint; reduced?: { x: bigint; modulus: bigint }; reason?: string };
export function crt(residues: Whole[], moduli: Whole[]): { x: bigint; modulus: bigint } | null;
export function solveLinearSystem(equations: Array<[Whole, Whole, Whole]>): { x: bigint; modulus: bigint } | null;
export function polyMod(poly: IntPoly, x: Whole, m: Whole): bigint;
export function hensel(poly: IntPoly, r: Whole, p: Whole, k: Whole): bigint;
export function rootsMod(poly: IntPoly, m: Whole): bigint[];
export function hasPrimitiveRoot(n: Whole): boolean;
export function primitiveRootMod(n: Whole): bigint | null;
export function order(a: Whole, n: Whole): bigint;
export function discreteLog(g: Whole, h: Whole, n: Whole): bigint | null;
export function valuation(n: Whole, p: Whole): number;

// ---- reciprocity

export function jacobi(a: Whole, n: Whole): -1 | 0 | 1;
export function legendre(a: Whole, p: Whole): -1 | 0 | 1;
export function kronecker(a: Whole, n: Whole): -1 | 0 | 1;
export function jacobiSteps(a: Whole, n: Whole): { value: -1 | 0 | 1; steps: string[] };
export function eulerCriterion(a: Whole, p: Whole): { power: bigint; value: -1 | 0 | 1 };
export function gaussLemma(a: Whole, p: Whole): { value: -1 | 1; mu: number; large: bigint[] };
export function quadraticResidues(n: Whole): bigint[];
export function sqrtsMod(a: Whole, n: Whole): bigint[];
export function sqrtMod(a: bigint, p: bigint): bigint | null;

// ---- sums of squares

export function primeAsTwoSquares(p: Whole): [bigint, bigint];
export function twoSquares(n: Whole): [bigint, bigint] | null;
export function twoSquaresAll(n: Whole): Array<[bigint, bigint]>;
export function twoSquaresReason(n: Whole): { possible: boolean; reason: string };
export function r2(n: Whole): bigint;
export function threeSquares(n: Whole): [bigint, bigint, bigint] | null;
export function fourSquares(n: Whole): [bigint, bigint, bigint, bigint];
export function r4(n: Whole): bigint;
export function rk(n: Whole, k: Whole): bigint;
export function pythagoreanTriples(limit: Whole): Array<[bigint, bigint, bigint]>;
export function triplesWithHypotenuse(c: Whole): Array<[bigint, bigint, bigint]>;
export function triplesWithLeg(a: Whole): Array<[bigint, bigint, bigint]>;

// ---- Diophantine equations

export function linearDiophantine(
  coefficients: Whole[],
  c: Whole,
): { solvable: true; particular: bigint[]; basis: bigint[][]; g: bigint } | { solvable: false; g: bigint; reason: string };
export function nonNegativeSolutions(a: Whole, b: Whole, c: Whole): { solutions: Array<[bigint, bigint]> };
export function frobenius(coins: Whole[]): bigint;
export function coinWays(amount: Whole, coins: Whole[]): bigint;
export type PellResult =
  | { kind: "infinite"; D: bigint; N: bigint; unit: [bigint, bigint]; fundamental: Array<[bigint, bigint]>; solutions(count?: number): Array<[bigint, bigint]> }
  | { kind: "finite"; solutions: Array<[bigint, bigint]> };
/** x² − D·y² = N. */
export function pellSolutions(D: Whole, N: Whole): PellResult;
export function normEquation(D: Whole, n: Whole): PellResult;
export function egyptianGreedy(q: RationalLike): bigint[];
export function egyptianShortest(q: RationalLike, options?: { maxDenominator?: bigint; maxTerms?: number }): bigint[] | null;

// ---- binary quadratic forms

/** A 2×2 integer matrix [[p, q], [r, s]]. */
export type Matrix2 = [[bigint, bigint], [bigint, bigint]];
/** The form a·x² + b·xy + c·y². */
export class Form {
  constructor(a: Whole, b: Whole, c: Whole);
  readonly a: bigint;
  readonly b: bigint;
  readonly c: bigint;
  readonly discriminant: bigint;
  at(x: Whole, y: Whole): bigint;
  transform(matrix: Matrix2): Form;
  isPrimitive(): boolean;
  isPositiveDefinite(): boolean;
  equals(g: Form): boolean;
  isReduced(): boolean;
  /** form is reduced and this.transform(matrix) equals it. */
  reduce(): { form: Form; matrix: Matrix2 };
  cycle(): Form[];
  toString(): string;
  toArray(): [bigint, bigint, bigint];
}
export type FormLike = Form | [Whole, Whole, Whole];
export function form(a: Whole | [Whole, Whole, Whole], b?: Whole, c?: Whole): Form;
export function equivalent(f: FormLike, g: FormLike): { equivalent: true; matrix: Matrix2 } | { equivalent: false; reason: string };
export function compose(f: FormLike, g: FormLike): Form;
export function principalForm(D: Whole): Form;
export function inverseForm(f: FormLike): Form;
export function reducedForms(D: Whole): Form[];
export interface ClassGroup {
  D: bigint;
  h: bigint;
  forms: Form[];
  identity: Form;
  classOf(f: FormLike): Form;
  compose(f: FormLike, g: FormLike): Form;
  inverse(f: FormLike): Form;
  order(f: FormLike): bigint;
  /** Invariant factors: the group is the product of cyclic groups of these orders. */
  structure: bigint[];
}
export function classGroup(D: Whole): ClassGroup;
export function classNumber(D: Whole): bigint;
export function representations(f: FormLike, n: Whole): Array<[bigint, bigint]>;
export function represent(f: FormLike, n: Whole): [bigint, bigint] | null;

// ---- p-adic numbers

/** A p-adic number computed lazily: approx(N) is a rational r with v_p(x − r) ≥ N. */
export class Padic {
  readonly p: bigint;
  /** A lower bound on the valuation. */
  readonly vmin: number;
  /** The exact rational value, when known. */
  readonly exact: Rational | null;
  approx(N: number): Rational;
  /** Throws "unknown-valuation" when x agrees with 0 to the search limit. */
  valuation(limit?: number): number;
  norm(): Rational;
  isZero(): boolean;
  mod(N: number | bigint): bigint;
  digits(count?: number): { valuation: number; digits: bigint[] };
  toString(count?: number): string;
  equalsTo(y: Padic | RationalLike, N: number): boolean;
  add(y: Padic | RationalLike): Padic;
  neg(): Padic;
  sub(y: Padic | RationalLike): Padic;
  mul(y: Padic | RationalLike): Padic;
  inv(): Padic;
  div(y: Padic | RationalLike): Padic;
  pow(k: Whole): Padic;
  isSquare(): boolean;
  shift(k: number | bigint): Padic;
  sqrt(): Padic;
  exp(): Padic;
  log(): Padic;
}
export function padic(x: Padic | RationalLike, p: Whole): Padic;
export function henselRoot(coefficients: Array<Padic | RationalLike>, r0: Whole, p: Whole): Padic;
export function teichmuller(a: Whole, p: Whole): Padic;
export function expansion(x: RationalLike, p: Whole): { valuation: number; preperiod: bigint[]; period: bigint[] };
export function padicValuation(x: RationalLike, p: Whole): number;
export function padicNorm(x: RationalLike, p: Whole): Rational;

// ---- continued fractions and Pell (moved from algebraic)

export function continuedFraction(x: AlgebraicReal | RationalLike, terms?: number): bigint[];
export function convergents(terms: bigint[]): Rational[];
export function fromContinuedFraction(terms: bigint[], period?: bigint[] | null): AlgebraicReal;
export function periodicContinuedFraction(x: AlgebraicReal): { preperiod: bigint[]; period: bigint[] };
export function pell(D: bigint | number): { x: bigint; y: bigint; negative: { x: bigint; y: bigint } | null };
export function fundamentalUnit(d: bigint | number): { unit: NumberFieldElement; norm: number };
export function fromQuadratic(p: bigint, q: bigint, d: bigint, r: bigint): AlgebraicReal;

// ---- the modular kernel, re-exported from algebraic

export { factorInteger, modPow, modInverse, multiplicativeOrder, primitiveRoot, lll, recognize, discreteLogMod } from "@ratmath/algebraic";

// ---- number fields: rings of integers, primes and ideals

export { NumberField, NumberFieldElement } from "@ratmath/algebraic";
/** A field: a NumberField, or the coefficients (constant first) of an irreducible polynomial. */
export type FieldLike = NumberField | Order | Array<bigint | number | string | Rational | [bigint, bigint]>;
/** An element of the field: a NumberFieldElement, a whole number or a Rational. */
export type ElementLike = NumberFieldElement | bigint | number | Rational;

export interface Round2Step {
  p: bigint;
  /** The exponent of p in the discriminant of the defining polynomial. */
  exponent: number;
  /** The order grew by index p^indexExponent at p. */
  indexExponent: number;
  rounds: number;
  note: string;
}
export class Order {
  readonly field: NumberField;
  readonly degree: number;
  /** The field discriminant. */
  readonly discriminant: bigint;
  /** [O_K : Z[θ]] for the integral generator θ = scale·α used inside. */
  readonly index: bigint;
  readonly scale: bigint;
  /** The discriminant of the minimal polynomial of θ. */
  readonly polyDiscriminant: bigint;
  /** [r₁, r₂]. */
  readonly signature: [number, number];
  /** What the Round 2 algorithm did at each prime whose square divides the polynomial discriminant. */
  readonly steps: Round2Step[];
  /** The integral basis, ω₁ = 1 first. */
  readonly basis: NumberFieldElement[];
  readonly traceMatrix: bigint[][];
  readonly ramifiedPrimes: bigint[];
  readonly unitIdeal: Ideal;
  coordinates(x: ElementLike): Rational[];
  element(coordinates: Array<bigint | number | Rational>): NumberFieldElement;
  contains(x: ElementLike): boolean;
  primesAbove(p: Whole): PrimeIdeal[];
  describePrime(p: Whole): { p: bigint; kind: "inert" | "splits" | "splits completely" | "ramified" | "totally ramified"; primes: PrimeIdeal[]; text: string };
  ideal(...generators: Array<ElementLike | Ideal | PrimeIdeal>): Ideal;
  toString(): string;
}
export class PrimeIdeal {
  readonly order: Order;
  /** The rational prime below. */
  readonly p: bigint;
  /** The ramification index. */
  readonly e: number;
  /** The residue degree. */
  readonly f: number;
  readonly norm: bigint;
  /** p and a second generator (p alone when the prime is inert). */
  readonly generators: NumberFieldElement[];
  /** The exponent of this prime in an element, an ideal or a prime ideal. */
  valuation(x: ElementLike | Ideal | PrimeIdeal): number;
  toIdeal(): Ideal;
  equals(other: PrimeIdeal): boolean;
  toString(): string;
}
export class Ideal {
  readonly order: Order;
  readonly norm: Rational;
  readonly isIntegral: boolean;
  readonly isOne: boolean;
  readonly isPrime: boolean;
  /** A Z-basis. */
  readonly basis: NumberFieldElement[];
  /** The smallest positive rational in the ideal. */
  readonly minimum: Rational;
  mul(other: Ideal | PrimeIdeal | ElementLike): Ideal;
  add(other: Ideal | PrimeIdeal | ElementLike): Ideal;
  div(other: Ideal | PrimeIdeal | ElementLike): Ideal;
  pow(k: number | bigint): Ideal;
  inverse(): Ideal;
  contains(x: ElementLike): boolean;
  equals(other: Ideal | PrimeIdeal | ElementLike): boolean;
  factor(): Array<{ prime: PrimeIdeal; exponent: number }>;
  toString(): string;
}
/** The ring of integers by the Round 2 algorithm; cached on the field. */
export function ringOfIntegers(K: FieldLike): Order;
export function integralBasis(K: FieldLike): NumberFieldElement[];
export function fieldDiscriminant(K: FieldLike): bigint;
export function primeDecomposition(K: FieldLike, p: Whole): PrimeIdeal[];

// ---- class groups and unit groups

export interface FieldOptions {
  /** false: skip the certification (the class number is then only known to divide the one shown). Default true. */
  certify?: boolean;
  /** The largest Minkowski bound accepted (default 2000). */
  maxBound?: number;
  /** The largest prime the saturation step goes to (default 20000). */
  maxIndexPrime?: number;
  /** How many elements may be tested for relations (default 400000). */
  budget?: number;
  /** How many times to gather more relations when certification fails (default 6). */
  rounds?: number;
}
export interface FieldCertificate {
  method: string;
  minkowskiBound?: bigint;
  relations?: number;
  /** R ≥ this (Friedman). */
  regulatorLowerBound?: number;
  /** The index of the units found in all units is at most this. */
  unitIndexBound?: number;
  /** The primes dividing the class number found, each given the full saturation test. */
  classNumberPrimes?: number[];
  /** How many primes passed a saturation test. */
  primesChecked?: number;
  /** How many relation elements were found as ℓ-th roots named by a failed saturation test. */
  rootsTaken?: number;
  narrowClassNumber?: bigint;
  formClassNumber?: bigint;
  note?: string;
}
export class ClassGroup {
  readonly field: NumberField;
  readonly order: Order;
  /** The class number. */
  readonly h: bigint;
  /** Invariant factors d₁ | d₂ | …; empty for the trivial group. */
  readonly structure: bigint[];
  /** One ideal of small norm per cyclic factor. */
  readonly generators: Ideal[];
  readonly proven: boolean;
  readonly certificate: FieldCertificate;
  readonly minkowskiBound: bigint;
  readonly factorBase: PrimeIdeal[];
  /** Exponents on the generators, each reduced modulo its invariant factor. */
  classOf(I: Ideal | PrimeIdeal | ElementLike): bigint[];
  isPrincipal(I: Ideal | PrimeIdeal | ElementLike): { principal: boolean; generator: NumberFieldElement | null };
  orderOf(I: Ideal | PrimeIdeal | ElementLike): bigint;
  toString(): string;
}
export class UnitGroup {
  readonly field: NumberField;
  readonly order: Order;
  /** r₁ + r₂ − 1. */
  readonly rank: number;
  readonly torsionOrder: number;
  readonly rootOfUnity: NumberFieldElement;
  readonly fundamentalUnits: NumberFieldElement[];
  /** To double precision. */
  readonly regulator: number;
  readonly proven: boolean;
  readonly certificate: FieldCertificate;
  isUnit(x: ElementLike): boolean;
  /** u = ζ^torsion · ∏ εᵢ^exponents[i]. */
  exponents(u: ElementLike): { torsion: number; exponents: bigint[] };
  toString(): string;
}
export function idealClassGroup(K: FieldLike, options?: FieldOptions): ClassGroup;
export function unitGroup(K: FieldLike, options?: FieldOptions): UnitGroup;
export function fieldClassNumber(K: FieldLike, options?: FieldOptions): bigint;
export function regulator(K: FieldLike, options?: FieldOptions): number;
/** ⌊the Minkowski bound⌋. */
export function minkowskiBound(K: FieldLike): bigint;
/** The ideal a·Z + ((−b + √D)/2)·Z of the quadratic field of discriminant D for a form (a, b, c). */
export function formToIdeal(f: Form | [Whole, Whole, Whole], K?: FieldLike): Ideal;

// ---- counting primes

/** π(x) by Meissel's formula with its pieces: π = φ + a − 1 − P2. For x ≤ 10¹². */
export function meissel(x: Whole): { pi: bigint; y: bigint; a: bigint; phi: bigint; P2: bigint };
/** Legendre's φ(x, a): the n ≤ x with no prime factor among the first a primes. */
export function legendrePhi(x: Whole, a: Whole): bigint;

// ---- Dirichlet characters

export class DirichletCharacter {
  constructor(q: Whole, m: Whole);
  readonly modulus: number;
  /** The Conrey number m of χ_q(m, ·). */
  readonly number: number;
  readonly order: number;
  readonly conductor: number;
  readonly isPrimitive: boolean;
  readonly isPrincipal: boolean;
  readonly isReal: boolean;
  readonly isEven: boolean;
  readonly parity: "even" | "odd";
  /** t with χ(n) = exp(2πi·t), in [0, 1); null when χ(n) = 0. */
  exponent(n: Whole): Rational | null;
  /** χ(n) in doubles. */
  value(n: Whole): { re: number; im: number };
  /** χ(n) for a real character: 0n, 1n or −1n. */
  at(n: Whole): bigint;
  exponents(): Array<Rational | null>;
  primitive(): DirichletCharacter;
  mul(other: DirichletCharacter): DirichletCharacter;
  pow(k: number | bigint): DirichletCharacter;
  conjugate(): DirichletCharacter;
  induce(modulus: Whole): DirichletCharacter;
  gaussSum(): { re: number; im: number };
  equals(other: DirichletCharacter): boolean;
  toString(): string;
}
export function dirichletCharacter(q: Whole, m?: Whole): DirichletCharacter;
export interface DirichletGroup {
  modulus: number;
  order: number;
  characters: DirichletCharacter[];
  principal: DirichletCharacter;
  readonly primitive: DirichletCharacter[];
  readonly real: DirichletCharacter[];
  /** Orders of cyclic factors. */
  structure: number[];
  /** Σ_χ χ(a)·conj(χ(b)). */
  orthogonality(a: Whole, b: Whole): { re: number; im: number };
  toString(): string;
}
export function dirichletGroup(q: Whole): DirichletGroup;
/** The Kronecker symbol (D/·) as a character mod |D|. */
export function kroneckerCharacter(D: Whole): DirichletCharacter;

// ---- zeta and L values

export function bernoulli(n: Whole): Rational;
export function bernoulliPolynomial(n: Whole, x: RationalLike): Rational;
export function generalizedBernoulli(k: Whole, chi: DirichletCharacter): Rational | NumberFieldElement;
/** A Rational for n ≤ 0; { coefficient, piPower } standing for coefficient·π^n for even n > 0. */
export function zetaExact(n: Whole): Rational | { coefficient: Rational; piPower: number; toString(): string };
/** L(n, χ) at an integer n ≤ 0: a Rational for a real character, else an element of a cyclotomic field. */
export function lValueExact(chi: DirichletCharacter, n: Whole): Rational | NumberFieldElement;
/** A point of the plane: [σ, t] or { re, im }, given exactly. */
export type PlanePoint = [RationalLike, RationalLike] | { re: RationalLike; im?: RationalLike };
export interface PlaneEnclosure {
  re: RationalInterval;
  im: RationalInterval;
}
export interface EnclosureOptions {
  /** The enclosure is narrower than 2^−bits (default 64). */
  bits?: number;
}
export function zetaValue(s: RationalLike, options?: EnclosureOptions): RationalInterval;
export function zetaValue(s: PlanePoint, options?: EnclosureOptions): PlaneEnclosure;
export function hurwitzZeta(s: RationalLike, a: RationalLike, options?: EnclosureOptions): RationalInterval;
export function hurwitzZeta(s: PlanePoint, a: RationalLike, options?: EnclosureOptions): PlaneEnclosure;
/** A RationalInterval when s is real and χ is real; { re, im } otherwise. */
export function lValue(chi: DirichletCharacter, s: RationalLike | PlanePoint, options?: EnclosureOptions): RationalInterval | PlaneEnclosure;

// ---- the plan

export interface PlannedPhase {
  phase: number;
  contents: string;
  /** Rough size in source lines, as text. */
  estimate: string;
}
/** Empty: every planned phase is built. */
export const plan: readonly PlannedPhase[];
