// Shared page furniture for the numtheory docs: the contents column, the theme
// switch and the code boxes. Add each topic page to PAGES as it is written.

import * as A from "./numtheory.js?v=c60683044c";
import { mountRunboxes } from "./runbox.js?v=c60683044c";

export { A };

export const PAGES = [
  { href: "index.html", title: "Overview" },
  { href: "reference.html", title: "Reference" },
  { href: "explore.html", title: "Explore" },
  { group: "Primes" },
  { href: "primes.html", title: "Primes and certificates", blurb: "Sieves, π(x), Baillie–PSW, and certificates checked by an independent verifier" },
  { href: "factoring.html", title: "Factoring", blurb: "Trial division, Pollard's rho and p − 1, with proven prime factors" },
  { href: "arithmetic-functions.html", title: "Arithmetic functions", blurb: "σ, τ, φ, μ, λ; perfect and amicable numbers; Dirichlet convolution and Möbius inversion" },
  { group: "Congruences" },
  { href: "congruences.html", title: "Congruences", blurb: "Linear congruences, CRT for any moduli, Hensel lifting, primitive roots and indices" },
  { href: "reciprocity.html", title: "Quadratic reciprocity", blurb: "Legendre, Jacobi and Kronecker symbols, with the reciprocity steps written out" },
  { href: "squares.html", title: "Sums of squares", blurb: "Two, three and four squares, Jacobi's counts, Pythagorean triples" },
  { group: "Diophantine equations" },
  { href: "diophantine.html", title: "Linear equations", blurb: "Every integer solution as a lattice, Frobenius numbers, Egyptian fractions" },
  { href: "continued-fractions.html", title: "Continued fractions", blurb: "Exact expansions, best approximations, periodic quadratic irrationals, units" },
  { href: "pell.html", title: "Pell equations", blurb: "x² − Dy² = N: every class of solutions, by Nagell's bounds" },
  { href: "quadratic-forms.html", title: "Quadratic forms", blurb: "Reduction with proof, cycles, Gauss composition, class groups and their structure" },
  { group: "p-adic numbers" },
  { href: "padic.html", title: "p-adic numbers", blurb: "Lazy exact digits, valuations, square and Hensel roots, exp and log" },
];

export const h = (tag, attrs = {}, ...kids) => {
  const el = document.createElement(tag);
  for (const [k, v] of Object.entries(attrs)) {
    if (v === null || v === undefined || v === false) continue;
    if (k === "class") el.className = v;
    else if (k.startsWith("on")) el.addEventListener(k.slice(2), v);
    else el.setAttribute(k, v === true ? "" : v);
  }
  for (const kid of kids.flat(Infinity)) if (kid !== null && kid !== undefined && kid !== false) el.append(kid instanceof Node ? kid : document.createTextNode(String(kid)));
  return el;
};

// ---- chrome -------------------------------------------------------------------

function theme() {
  let saved = null;
  try {
    saved = localStorage.getItem("ratmath-theme");
  } catch {
    /* storage may be blocked */
  }
  const apply = (t) => {
    if (t === "auto") document.documentElement.removeAttribute("data-theme");
    else document.documentElement.setAttribute("data-theme", t);
    try {
      localStorage.setItem("ratmath-theme", t);
    } catch {
      /* ignore */
    }
    for (const b of document.querySelectorAll(".theme button")) b.classList.toggle("on", b.dataset.t === t);
  };
  if (saved) apply(saved);
  return h("div", { class: "theme", role: "group", "aria-label": "Theme" },
    ["auto", "light", "dark"].map((t) => h("button", { "data-t": t, class: (saved ?? "auto") === t ? "on" : "", onclick: () => apply(t) }, t)));
}

function chrome() {
  const here = location.pathname.split("/").pop() || "index.html";
  const aside = document.querySelector("aside");
  if (!aside) return;
  const nav = h("nav", { class: "toc", "aria-label": "Contents" });
  const sel = h("select", { id: "toc-select", "aria-label": "Page", onchange: (e) => (location.href = e.target.value) });
  let group = null;
  for (const p of PAGES) {
    if (p.group) {
      nav.append(h("div", { class: "group" }, p.group));
      group = h("optgroup", { label: p.group });
      sel.append(group);
      continue;
    }
    nav.append(h("a", { href: p.href, class: p.href === here ? "on" : null, "aria-current": p.href === here ? "page" : null }, p.title));
    (group ?? sel).append(h("option", { value: p.href, selected: p.href === here }, p.title));
  }
  aside.append(
    h("header", { class: "brand" }, h("a", { href: "index.html" }, h("div", { class: "mark" }, "RatMath ", h("span", {}, "numtheory"))), h("small", {}, "Primes, congruences and equations in integers")),
    sel, nav, theme(),
  );
}

/** Topic cards for the overview page (pages that have a blurb). */
export function topicCards(el) {
  el.append(...PAGES.filter((p) => p.blurb).map((p) => h("a", { href: p.href }, h("strong", {}, p.title), h("span", {}, p.blurb))));
}

chrome();
mountRunboxes(document, { scope: A });
