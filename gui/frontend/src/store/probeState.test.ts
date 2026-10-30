/** @vitest-environment jsdom */
import type { Edge, Node } from "reactflow";
import { beforeEach, describe, expect, it } from "vitest";
import { PROBE_LINK_HANDLE, probeLinkId } from "../utils/probeLink";
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

// #856: a probe is display-only (graph_builder drops it before the solve), so
// nothing done to one may clear the flow results or the solve's status.
describe("probe operations keep the solve (#856)", () => {
	const solve = { success: true, ledger: { ok: true } };
	beforeEach(() => {
		useStore.setState({
			nodes: [
				node("a", "channel", { result: { m_dot: 1, success: true } }),
				node("b", "plenum", { result: { state: { P: 2e5 }, success: true } }),
				node("p", "probe", { target_id: "a", slot1_key: "m_dot" }),
			],
			edges: [
				{
					id: "f",
					source: "a",
					target: "b",
					data: { result: { m_dot: 1 } },
				} as Edge,
				{
					id: "w",
					source: "a",
					target: "b",
					type: "thermal",
					data: { type: "thermal", result: { Q: 5 } },
				} as Edge,
			],
			solveResults: solve,
			selectedProbeLinks: [],
		});
	});

	const expectSolveKept = () => {
		const { nodes, edges, solveResults } = useStore.getState();
		expect(nodes.find((n) => n.id === "a")?.data.result.m_dot).toBe(1);
		expect(nodes.find((n) => n.id === "b")?.data.result.state.P).toBe(2e5);
		expect(edges.find((e) => e.id === "f")?.data.result.m_dot).toBe(1);
		expect(edges.find((e) => e.id === "w")?.data.result.Q).toBe(5);
		expect(solveResults).toBe(solve);
	};

	it("adding a probe keeps every result, and it reads its target at once", () => {
		useStore.getState().onNodesChange([
			{
				type: "add",
				item: node("p2", "probe", { target_id: "b", slot1_key: "state.P" }),
			},
		]);
		expectSolveKept();
		const { nodes } = useStore.getState();
		const p2 = nodes.find((n) => n.id === "p2");
		const target = nodes.find((n) => n.id === p2?.data.target_id);
		expect(target?.data.result.state.P).toBe(2e5);
	});

	it("removing a probe keeps every result", () => {
		useStore.getState().onNodesChange([{ type: "remove", id: "p" }]);
		expect(useStore.getState().nodes.map((n) => n.id)).toEqual(["a", "b"]);
		expectSolveKept();
	});

	it("selecting, moving or measuring a probe keeps every result", () => {
		useStore.getState().onNodesChange([
			{ type: "select", id: "p", selected: true },
			{ type: "select", id: "a", selected: false },
			{ type: "position", id: "p", position: { x: 9, y: 9 } },
			{ type: "dimensions", id: "p", dimensions: { width: 140, height: 60 } },
		]);
		expectSolveKept();
	});

	it("retargeting, to a node or a wall, keeps every result", () => {
		const { setProbeTarget } = useStore.getState();
		setProbeTarget("p", "b");
		expect(useStore.getState().nodes[2].data.target_id).toBe("b");
		setProbeTarget("p", "w");
		expect(useStore.getState().nodes[2].data.target_id).toBe("w");
		expectSolveKept();
	});

	it("refuses a flow edge, another probe or itself as the target", () => {
		const { setProbeTarget } = useStore.getState();
		setProbeTarget("p", "f");
		setProbeTarget("p", "p");
		expect(useStore.getState().nodes[2].data.target_id).toBe("a");
	});

	it("a drag from the probe's link handle that ends on a handle sets its target", () => {
		useStore.getState().onConnect({
			source: "p",
			sourceHandle: PROBE_LINK_HANDLE,
			target: "b",
			targetHandle: null,
		});
		expect(useStore.getState().edges).toHaveLength(2);
		expect(useStore.getState().nodes[2].data.target_id).toBe("b");
		expectSolveKept();
	});

	it("selecting the link and deleting it clears target_id, nothing else", () => {
		const link = probeLinkId("p");
		const { onEdgesChange } = useStore.getState();
		onEdgesChange([{ type: "select", id: link, selected: true }]);
		expect(useStore.getState().selectedProbeLinks).toEqual([link]);
		onEdgesChange([{ type: "remove", id: link }]);
		const { nodes, edges, selectedProbeLinks } = useStore.getState();
		expect(nodes[2].data.target_id).toBeUndefined();
		expect(nodes[2].data.slot1_key).toBe("m_dot");
		expect(edges.map((e) => e.id)).toEqual(["f", "w"]);
		expect(selectedProbeLinks).toEqual([]);
		expectSolveKept();
	});

	it("deleting a probe with its link keeps every result", () => {
		// React Flow removes a deleted node's edges first, then the node.
		useStore
			.getState()
			.onEdgesChange([{ type: "remove", id: probeLinkId("p") }]);
		useStore.getState().onNodesChange([{ type: "remove", id: "p" }]);
		expect(useStore.getState().nodes.map((n) => n.id)).toEqual(["a", "b"]);
		expectSolveKept();
	});

	it("adding or removing a flow node still clears them", () => {
		useStore.getState().onNodesChange([{ type: "remove", id: "b" }]);
		const { nodes, edges, solveResults } = useStore.getState();
		expect(nodes[0].data.result).toBeUndefined();
		expect(edges[0].data.result).toBeUndefined();
		expect(solveResults).toBeNull();
	});
});
