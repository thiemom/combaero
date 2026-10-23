import { describe, expect, it } from "vitest";
import { draggedOffset, labelPosition, pathThroughBox } from "./labelOffset";

describe("labelPosition", () => {
	it("is the auto midpoint when no offset is stored", () => {
		expect(labelPosition(100, 50, undefined)).toEqual([100, 50]);
	});

	it("adds the stored offset, so the box follows the wall", () => {
		const offset = { dx: 30, dy: -20 };
		expect(labelPosition(100, 50, offset)).toEqual([130, 30]);
		expect(labelPosition(400, 10, offset)).toEqual([430, -10]);
	});
});

describe("draggedOffset", () => {
	it("keeps a click a click under the threshold", () => {
		expect(draggedOffset({ dx: 0, dy: 0 }, 2, 1, 1)).toBeNull();
	});

	it("converts screen pixels to flow units through the zoom", () => {
		expect(draggedOffset({ dx: 0, dy: 0 }, 40, -20, 2)).toEqual({
			dx: 20,
			dy: -10,
		});
	});

	it("adds to the offset the drag started from", () => {
		expect(draggedOffset({ dx: 5, dy: 5 }, 10, 0, 0.5)).toEqual({
			dx: 25,
			dy: 5,
		});
	});
});

describe("pathThroughBox", () => {
	const points = (d: string) =>
		[...d.matchAll(/-?\d+(\.\d+)?,-?\d+(\.\d+)?/g)].map((m) =>
			m[0].split(",").map(Number),
		);

	it("starts at the source, passes through the box and ends at the target", () => {
		const p = points(pathThroughBox(0, 0, "right", 200, 0, "left", 100, 80));
		expect(p[0]).toEqual([0, 0]);
		expect(p[3]).toEqual([100, 80]);
		expect(p[6]).toEqual([200, 0]);
	});

	it("leaves each handle the way it faces", () => {
		const p = points(pathThroughBox(0, 0, "right", 200, 0, "left", 100, 80));
		expect(p[1][0]).toBeGreaterThan(0); // out of the source to the right
		expect(p[1][1]).toBe(0);
		expect(p[5][0]).toBeLessThan(200); // into the target from the left
		expect(p[5][1]).toBe(0);
	});

	it("is smooth at the box: the two controls are mirrored about it", () => {
		const p = points(pathThroughBox(0, 0, "bottom", 50, 300, "top", -40, 150));
		expect(p[2][0] + p[4][0]).toBeCloseTo(2 * p[3][0]);
		expect(p[2][1] + p[4][1]).toBeCloseTo(2 * p[3][1]);
	});
});
