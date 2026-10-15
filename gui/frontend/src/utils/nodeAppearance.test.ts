import { describe, expect, it } from "vitest";
import { nodeTypes } from "../components/flowTypes";
import { NODE_APPEARANCE } from "./nodeAppearance";

describe("node appearance (#321 item 5)", () => {
	it("covers every node type the canvas registers, and nothing else", () => {
		expect(Object.keys(NODE_APPEARANCE).sort()).toEqual(
			Object.keys(nodeTypes).sort(),
		);
	});

	it("gives every type its own glyph", () => {
		// Mass/pressure boundary, effusion plate/liner and ejector/momentum
		// chamber each shared one before the table existed.
		const icons = Object.values(NODE_APPEARANCE).map((a) => a.icon);
		expect(new Set(icons).size).toBe(icons.length);
	});
});
