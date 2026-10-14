import { describe, expect, it } from "vitest";
import { nodeTypeLabel } from "./nodeUtils";

describe("nodeTypeLabel", () => {
	it("names the tee as the palette does, not by its internal type (#321)", () => {
		expect(nodeTypeLabel("mpce_tee")).toBe("Tee");
	});
	it("turns every underscore into a space", () => {
		expect(nodeTypeLabel("lossless_connection")).toBe("lossless connection");
		expect(nodeTypeLabel("a_b_c")).toBe("a b c");
	});
	it("is empty for a missing type", () => {
		expect(nodeTypeLabel(undefined)).toBe("");
		expect(nodeTypeLabel(null)).toBe("");
	});
});
