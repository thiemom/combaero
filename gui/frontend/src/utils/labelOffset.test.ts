import { describe, expect, it } from "vitest";
import { draggedOffset, labelPosition } from "./labelOffset";

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
