# Changelog

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
