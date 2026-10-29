import { describe, expect, it } from "vitest";
import { nodeTypes } from "../components/flowTypes";
import useStore from "../store/useStore";
import {
	CARD_ROW_DEFAULTS,
	type CardKind,
	cardRowKeys,
	isCardKind,
	isCardRow,
	toggleCardRow,
} from "./cardRows";
import { resolveKey } from "./diagnostics";
import { expectedResultKeys } from "./quantities";

const kinds = Object.keys(CARD_ROW_DEFAULTS) as CardKind[];

describe("CARD_ROW_DEFAULTS against what the backend emits (#775)", () => {
	it("every curated key is emitted for its kind", () => {
		// resultKeys.generated.json is the backend's per-kind key list: a
		// renamed or dropped key fails here instead of blanking the card.
		const missing = kinds.flatMap((kind) => {
			const emitted = new Set(expectedResultKeys(kind));
			return CARD_ROW_DEFAULTS[kind]
				.filter((k) => !emitted.has(k))
				.map((k) => `${kind}: ${k}`);
		});
		expect(missing).toEqual([]);
	});

	it("every curated kind is a palette kind", () => {
		for (const kind of kinds) expect(nodeTypes, kind).toHaveProperty(kind);
	});

	it("each kind has a short list without repeats", () => {
		for (const kind of kinds) {
			const rows = CARD_ROW_DEFAULTS[kind];
			expect(rows.length, kind).toBeGreaterThan(0);
			expect(rows.length, kind).toBeLessThanOrEqual(3);
			expect(new Set(rows).size, kind).toBe(rows.length);
		}
	});

	it("a node state key resolves on a node result", () => {
		const result = { state: { P: 2e5, Tt: 500 } };
		expect(resolveKey(result, CARD_ROW_DEFAULTS.plenum[0])).toBe(2e5);
		expect(resolveKey(result, CARD_ROW_DEFAULTS.plenum[1])).toBe(500);
	});
});

describe("card rows the user changes", () => {
	it("the default stands until the kind is overridden", () => {
		expect(cardRowKeys("channel", {})).toEqual(CARD_ROW_DEFAULTS.channel);
		expect(cardRowKeys("channel", { orifice: ["Re"] })).toEqual(
			CARD_ROW_DEFAULTS.channel,
		);
		expect(cardRowKeys("channel", { channel: ["Re"] })).toEqual(["Re"]);
	});

	it("a pin goes first and an unpin removes it", () => {
		const pinned = toggleCardRow(["m_dot", "dP"], "Re");
		expect(pinned).toEqual(["Re", "m_dot", "dP"]);
		expect(toggleCardRow(pinned, "m_dot")).toEqual(["Re", "dP"]);
	});

	it("a bare telemetry key matches its state. card key", () => {
		expect(isCardRow(["state.P"], "P")).toBe(true);
		expect(toggleCardRow(["state.P", "state.Tt"], "P")).toEqual(["state.Tt"]);
	});

	it("isCardKind is false for kinds without card rows", () => {
		expect(isCardKind("channel")).toBe(true);
		expect(isCardKind("probe")).toBe(false);
		expect(isCardKind("toString")).toBe(false);
		expect(isCardKind(undefined)).toBe(false);
	});

	it("the store keeps an override per kind and null resets it", () => {
		const { setCardRows } = useStore.getState();
		setCardRows("channel", ["Re", "m_dot"]);
		setCardRows("orifice", ["Cd"]);
		expect(useStore.getState().cardRows).toEqual({
			channel: ["Re", "m_dot"],
			orifice: ["Cd"],
		});
		setCardRows("channel", null);
		expect(useStore.getState().cardRows).toEqual({ orifice: ["Cd"] });
		setCardRows("orifice", null);
		expect(useStore.getState().cardRows).toEqual({});
	});
});
