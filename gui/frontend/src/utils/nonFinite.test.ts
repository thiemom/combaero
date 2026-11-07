/** @vitest-environment jsdom */
import type { Node } from "reactflow";
import { afterEach, beforeEach, describe, expect, it, vi } from "vitest";
import useStore from "../store/useStore";
import { describeNonFinite, findNonFinite } from "./nonFinite";
import { validateNetwork } from "./validation";

// #1038: a hand-edited file with `"pressure_ratio": 1e999` parses as
// Infinity (JSON has no Infinity token, but an overflowing literal is valid
// JSON). Carried into the canvas it was saved as null, which the reloaded
// file's client validation read as the default 4 and the backend refused.

const FILE_TEXT = `{
 "version": "1.1",
 "nodes": [
  {"id": "comp", "type": "compressor", "position": {"x": 0, "y": 0},
   "data": {"label": "comp", "pressure_ratio": 1e999}}
 ],
 "edges": []
}`;

const compressor = (data: Record<string, unknown>): Node => ({
	id: "comp",
	type: "compressor",
	position: { x: 0, y: 0 },
	data: { label: "comp", ...data },
});

describe("findNonFinite", () => {
	it("finds an overflowing literal parsed from a file", () => {
		const parsed = JSON.parse(FILE_TEXT);
		expect(parsed.nodes[0].data.pressure_ratio).toBe(Infinity);
		const found = findNonFinite(parsed);
		expect(found).toEqual([
			{ owner: "node 'comp'", field: "pressure_ratio", value: Infinity },
		]);
		expect(describeNonFinite(found)).toBe(
			"node 'comp' pressure_ratio = Infinity",
		);
	});

	it("walks nested data, edges and solver settings, skipping results", () => {
		const found = findNonFinite({
			nodes: [
				{
					id: "n",
					data: {
						composition: { custom_fractions: { N2: Number.NaN } },
						areas: [0.1, -Infinity],
						result: { P: Infinity },
					},
				},
			],
			edges: [{ id: "e", data: { R_wall: Infinity } }],
			solverSettings: { timeout: Infinity },
		});
		expect(found.map((f) => `${f.owner} ${f.field}`)).toEqual([
			"node 'n' composition.custom_fractions.N2",
			"node 'n' areas[1]",
			"edge 'e' R_wall",
			"solver settings timeout",
		]);
	});

	it("passes a finite network and null values", () => {
		expect(
			findNonFinite({
				nodes: [{ id: "n", data: { D: null, L: 1, label: "n" } }],
				edges: [],
				solverSettings: { timeout: 30, ambient_T: null },
			}),
		).toEqual([]);
	});
});

describe("loadNetwork refuses a non-finite number (#1038)", () => {
	beforeEach(() => {
		useStore.setState({ nodes: [], edges: [], loadNotice: null });
	});

	it("keeps the canvas and names the node and field", () => {
		const before = [compressor({ pressure_ratio: 3 })];
		useStore.setState({ nodes: before });
		useStore.getState().loadNetwork(JSON.parse(FILE_TEXT));
		const s = useStore.getState();
		expect(s.nodes).toBe(before);
		expect(s.loadNotice).toContain("Not loaded");
		expect(s.loadNotice).toContain("node 'comp' pressure_ratio = Infinity");
	});

	it("loads the same file with a finite value", () => {
		useStore
			.getState()
			.loadNetwork(JSON.parse(FILE_TEXT.replace("1e999", "5")));
		const s = useStore.getState();
		expect(s.nodes).toHaveLength(1);
		expect(s.nodes[0].data.pressure_ratio).toBe(5);
		expect(s.loadNotice).toBeNull();
	});
});

describe("saveNetwork never writes null for a number (#1038)", () => {
	let blobs: Blob[];
	beforeEach(() => {
		blobs = [];
		// jsdom has neither; the legacy download path uses both.
		URL.createObjectURL = vi.fn((b: Blob) => {
			blobs.push(b);
			return "blob:x";
		});
		URL.revokeObjectURL = vi.fn();
		useStore.setState({ edges: [], loadNotice: null });
	});
	afterEach(() => {
		vi.restoreAllMocks();
	});

	it("refuses to save an infinite pressure ratio and says so", async () => {
		useStore.setState({ nodes: [compressor({ pressure_ratio: Infinity })] });
		await useStore.getState().saveNetwork("net");
		expect(blobs).toHaveLength(0);
		const notice = useStore.getState().loadNotice;
		expect(notice).toContain("Not saved");
		expect(notice).toContain("node 'comp' pressure_ratio = Infinity");
	});

	it("saves a finite value as that number", async () => {
		useStore.setState({ nodes: [compressor({ pressure_ratio: 6 })] });
		await useStore.getState().saveNetwork("net");
		expect(blobs).toHaveLength(1);
		const saved = JSON.parse(await blobs[0].text());
		expect(saved.nodes[0].data.pressure_ratio).toBe(6);
		expect(useStore.getState().loadNotice).toBeNull();
	});
});

describe("validation reads null as the backend does (#1038)", () => {
	const prErrors = (data: Record<string, unknown>) =>
		validateNetwork([compressor(data)]).filter((e) =>
			e.includes("Pressure ratio"),
		);

	it("refuses an explicit null pressure ratio", () => {
		expect(prErrors({ pressure_ratio: null })).toHaveLength(1);
	});

	it("reads an absent pressure ratio as the default", () => {
		expect(prErrors({})).toHaveLength(0);
	});

	it("refuses null for the other plain-float fields it defaults", () => {
		const tee: Node = {
			id: "tee",
			type: "mpce_tee",
			position: { x: 0, y: 0 },
			data: { theta_deg: null, psi: null },
		};
		const errs = validateNetwork([tee]).filter((e) => e.includes("is null"));
		expect(errs).toHaveLength(2);
	});

	it("leaves inheritable nulls alone", () => {
		const ch: Node = {
			id: "ch",
			type: "channel",
			position: { x: 0, y: 0 },
			data: { D: null, L: 1 },
		};
		expect(validateNetwork([ch])).toEqual([]);
	});
});
