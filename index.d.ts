// Type declarations for @ratmath/numtheory.

import type { Rational } from "@ratmath/core";
import type { AlgebraicReal, NumberFieldElement } from "@ratmath/algebraic";

/** A whole number: a JavaScript integer, a BigInt, decimal text, or a Rational with denominator 1. Non-whole JS numbers are refused. */
export type Whole = bigint | number | string | Rational;
/** An exact rational: a whole number, fraction text ("3/7", "0.25") or a Rational. */
export type RationalLike = bigint | number | string | Rational;

/** Thrown for refusals; `code` says why ("bad-argument", "too-large", "unproven", "unknown-valuation", "not-a-square", "diverges", "too-hard"). */
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
export function factor(n: Whole, options?: { prove?: boolean }): PrimeFactor[];
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
export type Certificate = SmallCertificate | PocklingtonCertificate | BLSCertificate;
export type PrimalityProof =
  | { n: bigint; prime: true; proven: true; certificate: Certificate }
  | { n: bigint; prime: false; proven: true; factor?: bigint; witness?: bigint; reason?: string }
  | { n: bigint; prime: boolean; proven: false; reason: string };
export function provePrime(n: Whole, options?: { budget?: number }): PrimalityProof;
export function millerRabinWitness(n: bigint, a: bigint): boolean;
/** The independent checker: shares no code with the search. */
export function verifyCertificate(cert: Certificate): { valid: true } | { valid: false; reason: string };
export function verifyComposite(n: bigint, proof: { factor?: bigint; witness?: bigint }): { valid: true } | { valid: false; reason: string };

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

// ---- still to come

export interface PlannedPhase {
  phase: number;
  contents: string;
  /** Rough size in source lines, as text. */
  estimate: string;
}
/** The phases not built yet (4 to 7). */
export const plan: readonly PlannedPhase[];
