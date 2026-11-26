// What this package is planned to contain, phase by phase, from the RatMath
// library roadmap. Phases 1 to 4 are built; a phase is removed from this list
// when it ships, and the list goes away when it is empty.

export const plan = Object.freeze([
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
