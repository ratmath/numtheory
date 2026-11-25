// What this package is planned to contain, phase by phase, from the RatMath
// library roadmap. Nothing here is built yet; a phase is removed from this list
// when it ships, and the list goes away when it is empty.

export const plan = Object.freeze([
  {
    "phase": 1,
    "contents": "Primes with certificates, arithmetic functions, congruences, reciprocity, sums of squares, re-exports",
    "estimate": "1,400"
  },
  {
    "phase": 2,
    "contents": "Diophantine equations and binary quadratic forms",
    "estimate": "900"
  },
  {
    "phase": 3,
    "contents": "p-adic numbers",
    "estimate": "500"
  },
  {
    "phase": 4,
    "contents": "Elliptic curves; ECM and curve-based certificates",
    "estimate": "1,400"
  },
  {
    "phase": 5,
    "contents": "Integral bases, discriminants, prime ideal factorisation",
    "estimate": "900"
  },
  {
    "phase": 6,
    "contents": "Class groups and unit groups",
    "estimate": "1,500"
  },
  {
    "phase": 7,
    "contents": "π(x), characters, ζ and L values",
    "estimate": "800"
  }
].map(Object.freeze));
