// @ratmath/numtheory: number theory for RatMath, with proofs attached.
//
// Phase 1: primes with certificates, arithmetic functions, congruences,
// reciprocity, sums of squares, continued fractions and Pell (moved here from
// algebraic), and the modular toolkit re-exported from algebraic.
// Phase 2: Diophantine equations and binary quadratic forms.
// Phase 3: p-adic numbers.
// Phase 4: elliptic curves over Q and F_p, the elliptic curve method of
// factoring, and elliptic curve primality certificates.
// Phase 5: rings of integers of number fields (integral bases, discriminants)
// and the factorisation of primes and ideals into prime ideals.
// Phase 6: ideal class groups and unit groups, certified.
// Phase 7: π(x) by Meissel's formula, Dirichlet characters, and exact and
// enclosed values of zeta and L functions.

export { NumberTheoryError } from "./src/int.js";
import { publicInt } from "./src/int.js";
export const { isqrt, iroot, gcd, lcm, egcd } = publicInt;
export {
  primes, primePi, nthPrime, isPrime, isProbablePrime, nextPrime, prevPrime, primesBetween,
  pollardRho, pollardPm1, trialDivide, factor, formatFactors, provePrime, millerRabinWitness, PROVEN_BELOW,
} from "./src/primes.js";
export { verifyCertificate, verifyComposite } from "./src/verify.js";
export { EllipticCurve, ellipticCurve } from "./src/elliptic.js";
export { ecm, ecmCurve } from "./src/ecm.js";
export { classPolynomial } from "./src/ecpp.js";
export {
  divisors, sigma, tau, mobius, totient, jordan, carmichael, omega, bigOmega, liouville, radical, isSquarefree,
  aliquotSum, classify, isPerfect, aliquotSequence, amicablePairs, perfectNumbers, lucasLehmer,
  dirichlet, mobiusInversion, unit, one, id, isMultiplicative,
} from "./src/arith.js";
export {
  solveLinear, crt, solveLinearSystem, polyMod, hensel, rootsMod, hasPrimitiveRoot, primitiveRootMod, order,
  discreteLog, discreteLogMod, valuation,
} from "./src/congruence.js";
export { jacobi, legendre, kronecker, jacobiSteps, eulerCriterion, gaussLemma, quadraticResidues, sqrtsMod, sqrtMod } from "./src/quadratic.js";
export {
  primeAsTwoSquares, twoSquares, twoSquaresAll, twoSquaresReason, r2, threeSquares, fourSquares, r4, rk,
  pythagoreanTriples, triplesWithHypotenuse, triplesWithLeg,
} from "./src/squares.js";
export {
  linearDiophantine, nonNegativeSolutions, frobenius, coinWays, pellSolutions, normEquation, egyptianGreedy, egyptianShortest,
} from "./src/diophantine.js";
export {
  Form, form, equivalent, compose, principalForm, inverseForm, reducedForms, classGroup, classNumber, representations, represent,
} from "./src/forms.js";
export { Padic, padic, henselRoot, teichmuller, expansion, padicValuation, padicNorm } from "./src/padic.js";
export { continuedFraction, convergents, fromContinuedFraction, periodicContinuedFraction, pell, fundamentalUnit, fromQuadratic } from "./src/cf.js";
// the modular kernel stays in algebraic, which needs it for finite fields (roadmap 4.3)
export { factorInteger, modPow, modInverse, multiplicativeOrder, primitiveRoot, lll, recognize } from "@ratmath/algebraic";
export { Order, PrimeIdeal, Ideal, ringOfIntegers, integralBasis, fieldDiscriminant, primeDecomposition } from "./src/order.js";
export { ClassGroup, UnitGroup, idealClassGroup, unitGroup, fieldClassNumber, regulator, minkowskiBound, formToIdeal } from "./src/classgroup.js";
export { meissel, legendrePhi } from "./src/primecount.js";
export { DirichletCharacter, dirichletCharacter, dirichletGroup, kroneckerCharacter } from "./src/characters.js";
export { bernoulli, bernoulliPolynomial, generalizedBernoulli, zetaExact, lValueExact, zetaValue, hurwitzZeta, lValue } from "./src/lfunctions.js";
// number fields come from algebraic; re-exported so that this package is enough for the pages on them
export { NumberField, NumberFieldElement } from "@ratmath/algebraic";
export { plan } from "./src/plan.js";
