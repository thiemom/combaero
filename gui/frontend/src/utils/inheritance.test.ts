import { describe, expect, it } from "vitest";
import type { PreviewItem } from "../api";
import {
	inheritedField,
	inheritedNote,
	orificePipeNote,
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

// #1039: an in-pipe orifice whose bore is wider than its DOWNSTREAM pipe (the
// to-side channel that sets beta in reverse, #788) is refused by the resolve
// pass, though the upstream pipe is wider than the bore. The note follows the
// preview's verdict, not the upstream diameter alone.
describe("orificePipeNote", () => {
	const downstreamRefusal =
		"OrificeElement 'orf': bore 0.12 m is larger than the pipe it sits in " +
		"(channel 'c_out', D = 0.1 m). There is no such geometry. A bore equal " +
		"to the pipe is no orifice (zero loss); for a lumped Cd*A restriction " +
		"set in_pipe=False (no pipe diameter, beta = 0).";
	const preview = (
		status: "resolved" | "unresolved",
		reason = "",
	): PreviewItems => ({
		orf: {
			status,
			reason,
			fields: { pipe_diameter: { value: 0.2, source: "channel 'c_in'" } },
		},
	});
	const base = {
		nodeId: "orf",
		inPipe: true,
		pipeD: 0.2,
		passThrough: false,
		boreDiameter: 0.12,
	};

	it("shows the refusal against the downstream pipe", () => {
		expect(
			orificePipeNote({
				...base,
				preview: preview("unresolved", downstreamRefusal),
			}),
		).toEqual({ kind: "refused", reason: downstreamRefusal });
	});

	it("shows the upstream refusal once, as the solver's reason", () => {
		const upstream = downstreamRefusal.replace("c_out", "c_in");
		expect(
			orificePipeNote({
				...base,
				pipeD: 0.1,
				preview: preview("unresolved", upstream),
			}),
		).toEqual({ kind: "refused", reason: upstream });
	});

	it("leaves a chamber refusal to the orifice-family note", () => {
		expect(
			orificePipeNote({
				...base,
				preview: preview(
					"unresolved",
					"OrificeElement 'orf': bore 0.01 m^2 is larger than the momentum chamber it feeds",
				),
			}),
		).toBeNull();
	});

	it("flags a bore wider than the upstream pipe before the preview has", () => {
		expect(orificePipeNote({ ...base, pipeD: 0.1, preview: null })).toEqual({
			kind: "larger",
		});
	});

	it("reports the pass-through and says nothing for a fitting bore", () => {
		expect(
			orificePipeNote({
				...base,
				passThrough: true,
				preview: preview("resolved"),
			}),
		).toEqual({ kind: "equal" });
		expect(orificePipeNote({ ...base, preview: preview("resolved") })).toBe(
			null,
		);
		expect(
			orificePipeNote({ ...base, inPipe: false, preview: preview("resolved") }),
		).toBe(null);
	});
});
