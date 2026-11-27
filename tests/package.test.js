import { expect, test } from "bun:test";
import * as lib from "../index.js";

test("the package loads and has nothing left in its plan", () => {
  expect(Array.isArray(lib.plan)).toBe(true);
  expect(lib.plan.length).toBe(0);
  expect(Object.isFrozen(lib.plan)).toBe(true);
});
