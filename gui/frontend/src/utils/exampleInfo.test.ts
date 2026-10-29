/** @vitest-environment jsdom */
import { describe, expect, it } from "vitest";
import useStore from "../store/useStore";
import { labelText, parseExampleInfo } from "./exampleInfo";

// gui/examples/, the files the backend serves (#793).
const files = import.meta.glob<any>("../../../examples/*.json", {
	eager: true,
	import: "default",
});
const shipped: [string, any][] = Object.entries(files);

describe("shipped worked examples", () => {
	it("are there", () => {
		expect(shipped.length).toBeGreaterThanOrEqual(2);
	});

	it.each(shipped)("%s carries a full description and a label", (_, data) => {
		const info = parseExampleInfo(data.example);
		expect(info).not.toBeNull();
		expect(info?.source).not.toBe("");
		expect(info?.reproduced).not.toBe("");
		expect(info?.whereToLook).not.toBe("");
		expect(info?.labelNote).not.toBe("");
		expect(info?.script).toMatch(/^python\/examples\/papers\/.+\.py$/);
		expect(data.nodes.length).toBeGreaterThan(0);
		expect(data.edges.length).toBeGreaterThan(0);
	});
});

describe("parseExampleInfo", () => {
	it("is null for a plain network file", () => {
		expect(parseExampleInfo(undefined)).toBeNull();
		expect(parseExampleInfo({})).toBeNull();
	});

	it("never shows a label the file did not make", () => {
		expect(parseExampleInfo({ title: "x", label: "validated" })).toBeNull();
		expect(parseExampleInfo({ title: "x" })).toBeNull();
	});

	it("names the two labels as the validation policy does", () => {
		expect(labelText("fidelity")).toMatch(/own data/);
		expect(labelText("cross-source")).toMatch(/Cross-source/);
	});
});

describe("the example panel's state", () => {
	const example = {
		nodes: [{ id: "n1", position: { x: 0, y: 0 }, data: {} }],
		edges: [],
		example: { title: "T", label: "fidelity", source: "S" },
	};

	it("is set by loading an example and cleared by loading a plain file", () => {
		useStore.getState().loadNetwork(example);
		expect(useStore.getState().exampleInfo?.title).toBe("T");
		useStore.getState().loadNetwork({ nodes: [], edges: [] });
		expect(useStore.getState().exampleInfo).toBeNull();
	});

	it("is cleared by dismissing it and by Clear", () => {
		useStore.getState().loadNetwork(example);
		useStore.getState().dismissExampleInfo();
		expect(useStore.getState().exampleInfo).toBeNull();
		useStore.getState().loadNetwork(example);
		useStore.getState().clearCanvas();
		expect(useStore.getState().exampleInfo).toBeNull();
	});
});
