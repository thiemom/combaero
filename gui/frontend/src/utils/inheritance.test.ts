import { describe, expect, it } from "vitest";
import type { PreviewItem } from "../api";
import {
	inheritedField,
	inheritedSourceLabel,
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

describe("inheritedSourceLabel", () => {
	it("labels inherited sources and blanks user/default/empty", () => {
		expect(inheritedSourceLabel("channel 'ch_u'")).toBe("channel 'ch_u'");
		expect(inheritedSourceLabel("inherited from ch_u")).toBe(
			"inherited from ch_u",
		);
		expect(inheritedSourceLabel("user")).toBe("");
		expect(inheritedSourceLabel("default")).toBe("");
		expect(inheritedSourceLabel("")).toBe("");
		expect(inheritedSourceLabel(undefined)).toBe("");
	});
});
