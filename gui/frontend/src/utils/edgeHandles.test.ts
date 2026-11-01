/** @vitest-environment jsdom */
import { describe, expect, it } from "vitest";
import useStore from "../store/useStore";
import { repairEdgeHandles } from "./edgeHandles";
import table from "./nodeHandles.generated.json";

const HANDLES = table as Record<
	string,
	{ inlets: (string | null)[]; outlets: (string | null)[]; dynamic: boolean }
>;

// gui/examples/, the files the backend serves (#793).
const files = import.meta.glob<any>("../../../examples/*.json", {
	eager: true,
	import: "default",
});
const shipped: [string, any][] = Object.entries(files);

const flow = (
	id: string,
	source: string,
	target: string,
	sh?: string | null,
	th?: string | null,
) => ({
	id,
	source,
	target,
	sourceHandle: sh,
	targetHandle: th,
	type: "default",
});

describe("repairEdgeHandles (#901)", () => {
	it.each(shipped)("%s needs no repair", (_, data) => {
		const r = repairEdgeHandles(data.nodes, data.edges);
		expect(r.repaired).toEqual([]);
		expect(r.unresolved).toEqual([]);
	});

	it.each(shipped)(
		"%s: every positional flow edge is reconnected after its handles are broken",
		(_, data) => {
			const kind = new Map<string, string>(
				data.nodes.map((n: any) => [n.id, n.type]),
			);
			const broken = data.edges.map((e: any) =>
				e.type === "thermal"
					? e
					: { ...e, sourceHandle: "bogus", targetHandle: "bogus" },
			);
			const r = repairEdgeHandles(data.nodes, broken);
			expect(r.edges.length).toBe(data.edges.length);
			for (const e of r.edges) {
				if (e.type === "thermal") continue;
				for (const [node, key, side] of [
					[e.source, "sourceHandle", "outlets"],
					[e.target, "targetHandle", "inlets"],
				] as const) {
					const spec = HANDLES[kind.get(node) ?? ""];
					if (!spec || spec.dynamic) continue;
					const handles = spec[side];
					if (handles.every((h) => h === null || h.startsWith("flow-"))) {
						expect(handles).toContain((e as any)[key] ?? null);
					}
				}
			}
		},
	);

	it("moves a plenum end to its flow-right / flow-left handle", () => {
		const nodes = [
			{ id: "p", type: "plenum" },
			{ id: "c", type: "channel" },
			{ id: "b", type: "pressure_boundary" },
		];
		const r = repairEdgeHandles(nodes, [
			flow("in", "b", "p", "flow-source", "flow-target"),
			flow("out", "p", "c", "flow-source", "flow-target"),
		]);
		expect(r.edges[0].targetHandle).toBe("flow-left");
		expect(r.edges[1].sourceHandle).toBe("flow-right");
		expect(r.repaired).toHaveLength(2);
		expect(r.unresolved).toEqual([]);
	});

	it("gives a lossless connection the id-less handle it has", () => {
		const nodes = [
			{ id: "l", type: "lossless_connection" },
			{ id: "c", type: "channel" },
		];
		const r = repairEdgeHandles(nodes, [
			flow("e", "l", "c", "flow-source", "flow-target"),
		]);
		expect(r.edges[0].sourceHandle).toBeNull();
		expect(r.edges[0].targetHandle).toBe("flow-target");
	});

	it("does not guess a named port, and reports it", () => {
		const nodes = [
			{ id: "j", type: "ejector" },
			{ id: "b", type: "pressure_boundary" },
		];
		const r = repairEdgeHandles(nodes, [
			flow("e", "b", "j", "flow-source", "flow-target"),
		]);
		expect(r.edges[0].targetHandle).toBe("flow-target");
		expect(r.repaired).toEqual([]);
		expect(r.unresolved).toEqual(["e: j has no handle flow-target"]);
	});

	it("leaves an edge without handle ids alone (the canvas uses the first)", () => {
		const nodes = [
			{ id: "p", type: "plenum" },
			{ id: "c", type: "channel" },
		];
		const edges = [flow("e", "p", "c"), flow("f", "p", "c", null, null)];
		const r = repairEdgeHandles(nodes, edges);
		expect(r.edges).toEqual(edges);
		expect(r.repaired).toEqual([]);
	});

	it("leaves correct edges, thermal edges and computed-id nodes alone", () => {
		const nodes = [
			{ id: "a", type: "channel" },
			{ id: "b", type: "channel" },
			{ id: "h", type: "heat_source" },
		];
		const edges = [
			flow("ok", "a", "b", "flow-source", "flow-target"),
			{ ...flow("w", "a", "b", "bogus", "bogus"), type: "thermal" },
			flow("hs", "h", "a", "anything", "flow-target"),
		];
		const r = repairEdgeHandles(nodes, edges);
		expect(r.edges).toEqual(edges);
	});

	it("is reported in the load notice", () => {
		useStore.getState().loadNetwork({
			nodes: [
				{ id: "p", type: "plenum", position: { x: 0, y: 0 }, data: {} },
				{ id: "c", type: "channel", position: { x: 0, y: 0 }, data: {} },
			],
			edges: [flow("out", "p", "c", "flow-source", "flow-target")],
		});
		const s = useStore.getState();
		expect(s.edges[0].sourceHandle).toBe("flow-right");
		expect(s.loadNotice).toMatch(/Reconnected 1 flow edge end/);
	});
});
