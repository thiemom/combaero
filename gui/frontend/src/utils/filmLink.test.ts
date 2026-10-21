import { describe, expect, it } from "vitest";
import { defaultFilmSide, filmDuctId } from "./filmLink";

const flow = (source: string, target: string) => ({ source, target, data: {} });
const thermal = (source: string, target: string) => ({
	source,
	target,
	data: { type: "thermal" },
});

describe("defaultFilmSide", () => {
	it("picks the duct fed by the plenum the holes discharge into", () => {
		const edges = [
			flow("holes", "plenum"),
			flow("plenum", "hot"),
			flow("cin", "cold"),
		];
		expect(
			defaultFilmSide({ source: "hot", target: "cold" }, "holes", edges),
		).toBe("a");
		expect(
			defaultFilmSide({ source: "cold", target: "hot" }, "holes", edges),
		).toBe("b");
	});

	it("picks a duct directly downstream of the holes", () => {
		const edges = [flow("holes", "hot")];
		expect(
			defaultFilmSide({ source: "cold", target: "hot" }, "holes", edges),
		).toBe("b");
	});

	it("ignores thermal edges and falls back to A", () => {
		const edges = [thermal("holes", "hot")];
		expect(
			defaultFilmSide({ source: "cold", target: "hot" }, "holes", edges),
		).toBe("a");
	});
});

describe("filmDuctId", () => {
	it("maps the stored side to the wall edge's ends", () => {
		const w = { source: "x", target: "y" };
		expect(filmDuctId(w, "a")).toBe("x");
		expect(filmDuctId(w, undefined)).toBe("x");
		expect(filmDuctId(w, "b")).toBe("y");
	});
});
