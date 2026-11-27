# Changelog

## 0.3.0 (unreleased)

- `egcd`, `isqrt` and `iroot` took BigInts only (`egcd(252, 198)` threw
  "Invalid mix of BigInt and other type"), and `gcd`/`lcm` gave back numbers
  for numbers. All five now read any whole number (number, BigInt, numeric
  string, Core `Integer` or whole `Rational`) and return BigInts, like the rest
  of the package. Found by the Course Test (ratmath/shelf#7).

Phases 5, 6 and 7 of the roadmap: number fields and L-functions. The plan is
now fully built.

- `ringOfIntegers(K)` (Round 2): `basis`, `discriminant`, `index`,
  `coordinates`, `contains`, `traceMatrix`; also `integralBasis` and
  `fieldDiscriminant`. A field is a `NumberField` from algebraic (re-exported)
  or the coefficients of a monic irreducible polynomial, constant term first.
- Prime ideals above any prime, including those dividing the index:
  `primesAbove`, `describePrime`, `primeDecomposition`, `ramifiedPrimes`, with
  `e`, `f`, two generators and `valuation`.
- Ideals: `O.ideal(...)`, `mul`, `add`, `inverse`, `div`, `pow`, `norm`,
  `contains`, `factor`, and `sup`.
- `idealClassGroup(K)`: class number, structure, generators, `classOf`,
  `orderOf`, `isPrincipal` with a generator. Proven without the Riemann
  hypothesis: the Minkowski bound for generators, and saturation tests against
  Friedman's regulator bound for the index. A test that fails names a hidden
  ℓ-th power, whose root is taken in the field and added to the relations.
  Quadratic fields are checked against binary quadratic forms. `proven` and
  `certificate` say what was shown.
- `unitGroup(K)`: rank, roots of unity, fundamental units, `regulator`,
  `isUnit`, `exponents`. Also `fieldClassNumber`, `regulator`,
  `minkowskiBound`, `formToIdeal`.
- `primePi` now counts by Meissel's method above 10⁷, up to 10¹²; `meissel`
  and `legendrePhi` are exported.
- Dirichlet characters with Conrey numbering: `dirichletCharacter`,
  `dirichletGroup`, `kroneckerCharacter`; order, conductor, parity,
  primitivity, the inducing character, Gauss sums, products.
- Exact values: `bernoulli`, `bernoulliPolynomial`, `generalizedBernoulli`,
  `zetaExact`, `lValueExact`.
- Enclosures: `zetaValue`, `hurwitzZeta` and `lValue` at real or complex
  arguments, by Euler–Maclaurin summation in ball arithmetic with the
  remainder bounded.
- Seven new docs pages, and two new panels on the explore page.

## 0.2.0 (unreleased)

Phase 4 of the roadmap: elliptic curves.

- `ellipticCurve` over Q and over F_p, in short or long Weierstrass form: the
  group law, orders of points, and `toString`.
- Over Q: `torsion` (Nagell–Lutz, checked against Mazur's list), `shortModel`,
  `minimalModel` with its change of variables, `badPrimes`, `reduction`
  (good, split or non-split multiplicative, additive), `conductor` (refused
  when the reduction at 2 or 3 is additive), `rationalPoints`, and `height`:
  an enclosure of the canonical height with bounds proved for the curve.
- Over F_p: `count` by a character sum, by baby steps and giant steps with
  Mestre's twist, and by Schoof's algorithm; `trace`, `points`, `liftX`,
  `structure`, `divisionPolynomial`.
- `ecm` and `ecmCurve`: Lenstra's elliptic curve method (Montgomery curves,
  Suyama's parametrisation, two stages). `factor` uses it after rho and
  p − 1, and now factors F₇ = 2^128 + 1.
- Elliptic curve primality certificates (Atkin–Morain): `provePrime` builds a
  curve by complex multiplication when n − 1 cannot be factored far enough,
  and `verifyCertificate` checks the new `"elliptic"` kind independently.
  `classPolynomial` gives Hilbert class polynomials.
- `isqrt` no longer fails above 2^1024.
- `@ratmath/core` is now a peer dependency (`>=0.10.0 <0.12.0`): the package
  uses the program's own copy of Core instead of bringing its own, so Core 0.11
  works and there is never a second copy whose `Rational` is not the same class.

## 0.1.0 (unreleased)

Phases 1 to 3 of the roadmap.

- Primes: sieves, `primePi`, `nthPrime`, `isPrime` (proven below 2^64, by
  certificate above, otherwise refused), `nextPrime`, `prevPrime`,
  `primesBetween`.
- Certificates: `provePrime` gives Pocklington or Brillhart–Lehmer–Selfridge
  certificates (recursive), factors or Miller–Rabin witnesses;
  `verifyCertificate` and `verifyComposite` check them independently.
- Factoring: trial division, Pollard rho (Brent), Pollard p − 1, `factor` with
  proven prime factors.
- Arithmetic functions: divisors, σ_k, τ, φ, Jordan, Carmichael λ, μ, ω, Ω,
  Liouville, radical; perfect, amicable and aliquot sequences; Lucas–Lehmer;
  Dirichlet convolution and Möbius inversion.
- Congruences: linear congruences, CRT for any moduli, systems, polynomial
  roots mod m with Hensel lifting, primitive roots, orders, discrete logarithms.
- Reciprocity: Legendre, Jacobi and Kronecker symbols, the Jacobi computation
  step by step, Euler's criterion, Gauss's lemma, square roots mod n.
- Sums of squares: two squares through Gaussian integers, three and four
  squares, r₂, r₄, r_k, Pythagorean triples.
- Diophantine: linear equations with a lattice basis of solutions, Frobenius
  numbers, coin counting, generalised Pell equations by Nagell's bounds,
  Egyptian fractions.
- Continued fractions and Pell's equation, moved from `@ratmath/algebraic`.
- Binary quadratic forms: reduction with its matrix, equivalence, cycles,
  Gauss composition, class groups with their structure, representations.
- p-adic numbers: lazy exact arithmetic, valuations and norms, periodic
  expansions, square roots, Hensel roots, exp, log, Teichmüller
  representatives.
- Docs: eleven topic pages, explore and reference.
