/** @vitest-environment jsdom */
import type { Edge, Node } from "reactflow";
import { beforeEach, describe, expect, it } from "vitest";
import useStore, { PROBE_HANDLE } from "./useStore";

// What the Probe node and the Inspector's Live Telemetry read is each node's
// data.result. These pin when that result is current, kept or dropped.

const node = (id: string, type: string, data: any = {}): Node => ({
	id,
	type,
	position: { x: 0, y: 0 },
	data,
});

const solved = (nodeIds: string[], extra: any = {}) => ({
	success: true,
	message: "",
	node_results: {},
	element_results: Object.fromEntries(
		nodeIds.map((id) => [id, { m_dot: 0.1, P_in: 2e5, success: true }]),
	),
	edge_results: {},
	...extra,
});

describe("results after a solve", () => {
	beforeEach(() => {
		useStore.setState({
			nodes: [
				node("a", "channel", { result: { m_dot: 9, success: true } }),
				node("b", "channel", { result: { m_dot: 9, success: true } }),
				node("p", "probe", { target_id: "a", slot1_key: "m_dot" }),
			],
			edges: [
				{
					id: "w",
					source: "a",
					target: "b",
					type: "thermal",
					data: { type: "thermal", result: { Q: 5 } },
				} as Edge,
			],
			solveResults: null,
			displaySettings: ["m_dot"],
			knownFields: ["m_dot"],
		});
	});

	it("a failed request drops every stored result instead of showing it as current", () => {
		useStore.getState().setSolveResults({
			success: false,
			message: "Backend unreachable",
		});
		const { nodes, edges } = useStore.getState();
		expect(nodes.find((n) => n.id === "a")?.data.result).toBeUndefined();
		expect(nodes.find((n) => n.id === "b")?.data.result).toBeUndefined();
		expect(edges[0].data.result).toBeUndefined();
		// The probe's own settings are untouched.
		expect(nodes.find((n) => n.id === "p")?.data.target_id).toBe("a");
	});

	it("a node the response has no result for loses its old one", () => {
		useStore.getState().setSolveResults(solved(["a"]));
		const { nodes } = useStore.getState();
		expect(nodes.find((n) => n.id === "a")?.data.result.m_dot).toBe(0.1);
		expect(nodes.find((n) => n.id === "b")?.data.result).toBeUndefined();
	});

	it("clearing the status banner (null) keeps the results", () => {
		useStore.getState().setSolveResults(null);
		expect(useStore.getState().nodes[0].data.result.m_dot).toBe(9);
		expect(useStore.getState().edges[0].data.result.Q).toBe(5);
	});

	it("a wall result takes the solve's success flag", () => {
		useStore
			.getState()
			.setSolveResults(
				solved(["a", "b"], { success: false, edge_results: { w: { Q: 1 } } }),
			);
		expect(useStore.getState().edges[0].data.result).toEqual({
			Q: 1,
			success: false,
		});
	});

	it("a field the user unchecked stays unchecked at the next solve", () => {
		const { setSolveResults, setDisplaySettings } = useStore.getState();
		setSolveResults(solved(["a", "b"]));
		expect(useStore.getState().displaySettings).toContain("P_in");
		setDisplaySettings(
			useStore.getState().displaySettings.filter((k) => k !== "P_in"),
		);
		setSolveResults(solved(["a", "b"]));
		expect(useStore.getState().displaySettings).not.toContain("P_in");
		// A field never seen before is still promoted.
		setSolveResults(
			solved(["a"], {
				element_results: { a: { m_dot: 1, dPt: 3, success: true } },
			}),
		);
		expect(useStore.getState().displaySettings).toContain("dPt");
		expect(useStore.getState().displaySettings).not.toContain("P_in");
	});
});

describe("pointing a probe at a node", () => {
	beforeEach(() => {
		useStore.setState({
			nodes: [
				node("a", "channel", { result: { m_dot: 1, success: true } }),
				node("p", "probe", {}),
			],
			edges: [],
			solveResults: { success: true },
		});
	});

	it("a connection drawn to the probe sets its target, not a flow edge", () => {
		useStore.getState().onConnect({
			source: "a",
			sourceHandle: null,
			target: "p",
			targetHandle: PROBE_HANDLE,
		});
		const { nodes, edges } = useStore.getState();
		expect(edges).toEqual([]);
		expect(nodes.find((n) => n.id === "p")?.data.target_id).toBe("a");
		// Display-only: the solution it shows is kept.
		expect(nodes.find((n) => n.id === "a")?.data.result.m_dot).toBe(1);
	});

	it("re-targeting through the Inspector keeps every result", () => {
		// Inspector's target select calls updateNodeData(..., true).
		useStore.getState().updateNodeData("p", { target_id: "a" }, true);
		expect(useStore.getState().nodes[0].data.result.m_dot).toBe(1);
	});
});
