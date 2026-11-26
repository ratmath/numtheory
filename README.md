# @ratmath/numtheory

Primes, congruences and equations in integers for [RatMath](https://ratmath.com/), part of the
[RatMath Shelf](https://ratmath.com/shelf/).

**Status: phases 1 to 4 built.** Primes with certificates, factoring,
arithmetic functions, congruences, reciprocity, sums of squares, linear and Pell
equations, continued fractions, binary quadratic forms, p-adic numbers, and elliptic
curves with the factoring method and primality certificates built on them. The
rest of the plan is in the
[library roadmap](https://github.com/ratmath/ratmath/blob/main/docs/design/library-roadmap.md#4-a-number-theory-package-ratmathnumtheory).
Docs: [ratmath.com/numtheory](https://ratmath.com/numtheory/).

```js
import { provePrime, verifyCertificate, factor, formatFactors, pellSolutions, classGroup, padic } from "@ratmath/numtheory";

const r = provePrime(2n ** 127n - 1n);     // a Pocklington certificate
verifyCertificate(r.certificate);          // { valid: true }, from an independent checker
formatFactors(factor(2n ** 64n + 1n));     // "274177 · 67280421310721"
pellSolutions(2, 7).solutions(4);          // [[3n, 1n], [5n, 3n], [13n, 9n], [27n, 19n]]
classGroup(-3299).structure;               // [3n, 9n]
padic(-1, 5n).sqrt().toString(12);         // "…223032431212"
formatFactors(factor(2n ** 128n + 1n));    // "59649589127497217 · 5704689200685129054721", by elliptic curves
```

"Prime" means proven: below 2^64 by Baillie–PSW, which has no exceptions
there, and above by a Pocklington, Brillhart–Lehmer–Selfridge or elliptic
curve certificate that `src/verify.js`, which imports nothing from the rest of the package,
checks. When no certificate can be found `isPrime` refuses with the reason.

## Phases

| Phase | Contents | Status |
|---|---|---|
| 1 | Primes with certificates, arithmetic functions, congruences, reciprocity, sums of squares, re-exports | built |
| 2 | Diophantine equations and binary quadratic forms | built |
| 3 | p-adic numbers | built |
| 4 | Elliptic curves; ECM and curve-based certificates | built |
| 5 | Integral bases, discriminants, prime ideal factorisation | planned |
| 6 | Class groups and unit groups | planned |
| 7 | π(x), characters, ζ and L values | planned |

It depends on Core and algebraic (vendored as a git submodule under
`vendor/algebraic`; clone with `--recursive`). Continued fractions, Pell's
equation and the Chinese remainder theorem moved here from algebraic, which
keeps the modular kernel it needs for finite fields; those functions are
re-exported.

## Layout

- `index.js`, `src/`: the library; `index.d.ts` the types.
- `tests/`: `bun run test`.
- `docs/src/*.html`: the docs pages (overview, explore, reference, then one
  page per topic as topics are built). `bun run docs:build` writes `docs/`,
  typesets the mathematics and runs the code box at the end of every page.
  `docs/*.html` and `docs/assets/numtheory.js` are build output.
- `docs/assets/runbox.js` and `scripts/docs-kit.js` are shared with the other
  RatMath docs sites and should be kept in step with them.

## Development

```sh
git submodule update --init
bun install
bun run test
bun run docs:build
```

Each topic that is built gets its own docs page with usage, the mathematics
explained, and a runnable code box at the end, and is then listed on the
[Shelf page](https://ratmath.com/shelf/).

MIT licensed.
