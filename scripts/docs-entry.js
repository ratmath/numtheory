// Docs bundle entry: the library plus the Core pieces the pages use directly.
export * from "../index.js";
export { Rational, Integer, RationalInterval } from "@ratmath/core";
// for the continued-fractions page: real algebraic numbers to expand
export { fromRational, sqrt, nthRoot, format } from "@ratmath/algebraic";
