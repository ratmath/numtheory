# @ratmath/numtheory

Primes, congruences and equations in integers for [RatMath](https://ratmath.com/), part of the
[RatMath Shelf](https://ratmath.com/shelf/).

**Status: set up, nothing built yet.** This repository holds the package
skeleton and the docs site. The plan is in the
[library roadmap](https://github.com/ratmath/ratmath/blob/main/docs/design/library-roadmap.md#4-a-number-theory-package-ratmathnumtheory).

From a first course to algebraic number theory: primes that come with a certificate, congruences and reciprocity, Diophantine equations, p-adic numbers, elliptic curves, and the arithmetic of number fields.

## Planned

| Phase | Contents |
|---|---|
| 1 | Primes with certificates, arithmetic functions, congruences, reciprocity, sums of squares, re-exports |
| 2 | Diophantine equations and binary quadratic forms |
| 3 | p-adic numbers |
| 4 | Elliptic curves; ECM and curve-based certificates |
| 5 | Integral bases, discriminants, prime ideal factorisation |
| 6 | Class groups and unit groups |
| 7 | π(x), characters, ζ and L values |

It will depend on: Core and algebraic (the modular kernel, continued fractions, number fields), analytic for the analytic part.

## Layout

- `index.js`, `src/`: the library. For now it exports only `plan`, the list
  above, so that the tests and the docs build have something to load.
- `tests/`: `bun run test`.
- `docs/src/*.html`: the docs pages (overview, explore, reference, then one
  page per topic as topics are built). `bun run docs:build` writes `docs/`,
  typesets the mathematics and runs the code box at the end of every page.
  `docs/*.html` and `docs/assets/numtheory.js` are build output.
- `docs/assets/runbox.js` and `scripts/docs-kit.js` are shared with the other
  RatMath docs sites and should be kept in step with them.

## Development

```sh
bun install
bun run test
bun run docs:build
```

Each topic that is built gets its own docs page with usage, the mathematics
explained, and a runnable code box at the end, and is then listed on the
[Shelf page](https://ratmath.com/shelf/).

MIT licensed.
