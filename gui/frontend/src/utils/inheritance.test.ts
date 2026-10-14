import { describe, expect, it } from "vitest";
import type { PreviewItem } from "../api";
import {
	inheritedField,
	inheritedNote,
	type PreviewItems,
	unresolvedReason,
} from "./inheritance";

const items: PreviewItems = {
	ch_u: {
		status: "resolved",
		reason: "",
		fields: { diameter: { value: 0.1, source: "user" } },
	},
	ch_inh: {
		status: "resolved",
		reason: "",
		fields: { diameter: { value: 0.1, source: "channel 'ch_u'" } },
	},
	ch_bad: {
		status: "unresolved",
		reason: "ChannelElement 'ch_bad': cannot determine diameter",
		fields: { diameter: { value: null, source: "" } },
	},
} satisfies Record<string, PreviewItem>;

describe("inheritedField", () => {
	it("returns the resolved value and source", () => {
		expect(inheritedField(items, "ch_inh", "diameter")).toEqual({
			value: 0.1,
			source: "channel 'ch_u'",
		});
	});

	it("is undefined for a missing node or field", () => {
		expect(inheritedField(items, "nope", "diameter")).toBeUndefined();
		expect(inheritedField(items, "ch_u", "F0")).toBeUndefined();
		expect(inheritedField(null, "ch_u", "diameter")).toBeUndefined();
	});
});

describe("unresolvedReason", () => {
	it("gives the reason only for an unresolved node", () => {
		expect(unresolvedReason(items, "ch_bad")).toContain("cannot determine");
		expect(unresolvedReason(items, "ch_u")).toBeUndefined();
		expect(unresolvedReason(items, "nope")).toBeUndefined();
	});
});

describe("inheritedNote", () => {
	it("phrases each source shape so it reads correctly", () => {
		expect(inheritedNote("channel 'ch_u'")).toBe(
			"Inherited from channel 'ch_u'",
		);
		expect(inheritedNote("node 'p1'")).toBe("Inherited from node 'p1'");
		// Already a phrase: capitalized, not wrapped ("Inherited from inherited ...").
		expect(inheritedNote("inherited from ch_u")).toBe("Inherited from ch_u");
		expect(inheritedNote("derived from F_C / psi")).toBe(
			"Derived from F_C / psi",
		);
	});

	it("is empty for user-set, defaulted, or missing sources", () => {
		expect(inheritedNote("user")).toBe("");
		expect(inheritedNote("default")).toBe("");
		expect(inheritedNote("")).toBe("");
		expect(inheritedNote(undefined)).toBe("");
	});
});
