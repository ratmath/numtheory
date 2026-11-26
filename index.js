// @ratmath/numtheory: number theory for RatMath, with proofs attached.
//
// Phase 1: primes with certificates, arithmetic functions, congruences,
// reciprocity, sums of squares, continued fractions and Pell (moved here from
// algebraic), and the modular toolkit re-exported from algebraic.
// Phase 2: Diophantine equations and binary quadratic forms.
// Phase 3: p-adic numbers.
// `plan` lists what is still to come.

export { NumberTheoryError, isqrt, iroot, gcd, lcm, egcd } from "./src/int.js";
export {
  primes, primePi, nthPrime, isPrime, isProbablePrime, nextPrime, prevPrime, primesBetween,
  pollardRho, pollardPm1, trialDivide, factor, formatFactors, provePrime, millerRabinWitness, PROVEN_BELOW,
} from "./src/primes.js";
export { verifyCertificate, verifyComposite } from "./src/verify.js";
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
export { plan } from "./src/plan.js";
