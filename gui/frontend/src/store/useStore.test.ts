/** @vitest-environment jsdom */
import type { Edge, EdgeChange, Node, NodeChange } from "reactflow";
import { beforeEach, describe, expect, it } from "vitest";
import useStore, { stripResults } from "./useStore";

// Regression coverage for the topology result invalidation shipped with the
// GUI integrity fixes: stored per-node/edge results warm-start the next
// solve, so structural edits must drop them while position/selection
// changes keep them.

const nodeWithResult = (id: string): Node => ({
	id,
	position: { x: 0, y: 0 },
	data: { label: id, diameter: 0.05, result: { P: 101325, success: true } },
});

const edgeWithResult = (id: string, source: string, target: string): Edge => ({
	id,
	source,
	target,
	data: { result: { m_dot: 0.2 } },
});

describe("stripResults", () => {
	it("removes data.result and preserves sibling keys", () => {
		const stripped = stripResults([nodeWithResult("n1")]);
		expect(stripped[0].data.result).toBeUndefined();
		expect(stripped[0].data.label).toBe("n1");
		expect(stripped[0].data.diameter).toBe(0.05);
	});

	it("returns items without a result unchanged", () => {
		const plain: Node = {
			id: "n1",
			position: { x: 0, y: 0 },
			data: { label: "n1" },
		};
		const noData: Node = {
			id: "n2",
			position: { x: 0, y: 0 },
			data: undefined,
		};
		const stripped = stripResults([plain, noData]);
		expect(stripped[0]).toBe(plain);
		expect(stripped[1]).toBe(noData);
	});
});

describe("useStore topology result invalidation", () => {
	beforeEach(() => {
		useStore.setState({
			nodes: [nodeWithResult("n1"), nodeWithResult("n2")],
			edges: [edgeWithResult("e1", "n1", "n2")],
			solveResults: { success: true },
		});
	});

	it("keeps results on position changes", () => {
		const changes: NodeChange[] = [
			{
				id: "n1",
				type: "position",
				position: { x: 10, y: 20 },
				dragging: false,
			},
		];
		useStore.getState().onNodesChange(changes);
		const { nodes, edges } = useStore.getState();
		expect(nodes[0].position).toEqual({ x: 10, y: 20 });
		expect(nodes[0].data.result).toBeDefined();
		expect(edges[0].data.result).toBeDefined();
	});

	it("keeps results on selection changes", () => {
		const changes: NodeChange[] = [
			{ id: "n1", type: "select", selected: true },
		];
		useStore.getState().onNodesChange(changes);
		expect(useStore.getState().nodes[0].data.result).toBeDefined();
		// The solve's status too (#856): a click used to drop it, and with it
		// the Inspector's second-law and Mach flags for the clicked node.
		expect(useStore.getState().solveResults).toEqual({ success: true });
	});

	it("strips all node and edge results when a node is removed", () => {
		const changes: NodeChange[] = [{ id: "n2", type: "remove" }];
		useStore.getState().onNodesChange(changes);
		const { nodes, edges, solveResults } = useStore.getState();
		expect(nodes.map((n) => n.id)).toEqual(["n1"]);
		expect(nodes[0].data.result).toBeUndefined();
		expect(edges[0].data.result).toBeUndefined();
		expect(solveResults).toBeNull();
	});

	it("strips all node and edge results when a node is added", () => {
		const changes: NodeChange[] = [
			{ type: "add", item: { id: "n3", position: { x: 5, y: 5 }, data: {} } },
		];
		useStore.getState().onNodesChange(changes);
		const { nodes, edges } = useStore.getState();
		// applyNodeChanges prepends added items.
		expect(nodes.map((n) => n.id)).toEqual(["n3", "n1", "n2"]);
		expect(nodes[1].data.result).toBeUndefined();
		expect(nodes[2].data.result).toBeUndefined();
		expect(edges[0].data.result).toBeUndefined();
	});

	it("strips all node and edge results when an edge is removed", () => {
		const changes: EdgeChange[] = [{ id: "e1", type: "remove" }];
		useStore.getState().onEdgesChange(changes);
		const { nodes, edges } = useStore.getState();
		expect(edges).toEqual([]);
		expect(nodes[0].data.result).toBeUndefined();
		expect(nodes[1].data.result).toBeUndefined();
	});

	it("strips all node and edge results on a new connection", () => {
		useStore.getState().onConnect({
			source: "n2",
			target: "n1",
			sourceHandle: null,
			targetHandle: null,
		});
		const { nodes, edges, solveResults } = useStore.getState();
		expect(edges).toHaveLength(2);
		expect(nodes[0].data.result).toBeUndefined();
		expect(nodes[1].data.result).toBeUndefined();
		expect(edges[0].data.result).toBeUndefined();
		expect(solveResults).toBeNull();
	});
});

describe("loadNetwork legacy migration", () => {
	it("maps the retired incompressible_warmstart strategy to outletref_warmstart", () => {
		useStore.getState().loadNetwork({
			nodes: [nodeWithResult("n1")],
			edges: [],
			solverSettings: {
				...useStore.getState().solverSettings,
				init_strategy: "incompressible_warmstart",
			},
		});
		expect(useStore.getState().solverSettings.init_strategy).toBe(
			"outletref_warmstart",
		);
	});

	it("moves the old omega_rpm to shaft_speed_rpm (#773)", () => {
		const { shaft_speed_rpm: _drop, ...oldSettings } =
			useStore.getState().solverSettings;
		useStore.getState().loadNetwork({
			nodes: [
				{
					id: "v1",
					type: "vortex",
					position: { x: 0, y: 0 },
					data: { label: "v1", omega_rpm: 3000 },
				},
			],
			edges: [],
			solverSettings: { ...oldSettings, omega_rpm: 1500 },
		});
		const { nodes, solverSettings } = useStore.getState();
		expect(nodes[0].data.shaft_speed_rpm).toBe(3000);
		expect(nodes[0].data.omega_rpm).toBeUndefined();
		expect(solverSettings.shaft_speed_rpm).toBe(1500);
		expect("omega_rpm" in solverSettings).toBe(false);
	});

	it("moves a Mynard tee's P_jct guess to Pt_jct, not the ideal tee's (#773)", () => {
		const tee = (id: string, junction_model?: string): Node => ({
			id,
			type: "mpce_tee",
			position: { x: 0, y: 0 },
			data: { label: id, junction_model, initial_guess: { P_jct: 1.2e5 } },
		});
		useStore.getState().loadNetwork({
			nodes: [tee("m"), tee("k", "constant_k"), tee("i", "ideal_momentum")],
			edges: [],
			solverSettings: useStore.getState().solverSettings,
		});
		const guess = (id: string) =>
			useStore.getState().nodes.find((n) => n.id === id)?.data.initial_guess;
		expect(guess("m")).toEqual({ Pt_jct: 1.2e5 });
		expect(guess("k")).toEqual({ Pt_jct: 1.2e5 });
		expect(guess("i")).toEqual({ P_jct: 1.2e5 });
	});

	it("gives a file saved before #864 the stall handover, and keeps an opt-out", () => {
		const { stall_handover: _drop, ...oldSettings } =
			useStore.getState().solverSettings;
		useStore.getState().loadNetwork({
			nodes: [nodeWithResult("n1")],
			edges: [],
			solverSettings: oldSettings,
		});
		expect(useStore.getState().solverSettings.stall_handover).toBe(true);
		useStore.getState().loadNetwork({
			nodes: [nodeWithResult("n1")],
			edges: [],
			solverSettings: { ...oldSettings, stall_handover: false },
		});
		expect(useStore.getState().solverSettings.stall_handover).toBe(false);
	});

	it("leaves current strategies untouched", () => {
		useStore.getState().loadNetwork({
			nodes: [nodeWithResult("n1")],
			edges: [],
			solverSettings: {
				...useStore.getState().solverSettings,
				init_strategy: "analytical_pt_prop",
			},
		});
		expect(useStore.getState().solverSettings.init_strategy).toBe(
			"analytical_pt_prop",
		);
	});
});

describe("loadNetwork drops walls on surfaceless nodes (#720)", () => {
	const settings = () => useStore.getState().solverSettings;
	const node = (id: string, type: string): Node => ({
		id,
		type,
		position: { x: 0, y: 0 },
		data: { label: id },
	});

	it("removes the wall and says so", () => {
		useStore.getState().loadNetwork({
			nodes: [node("ch", "channel"), node("or", "orifice")],
			edges: [
				{ id: "w1", source: "or", target: "ch", data: { type: "thermal" } },
				{ id: "f1", source: "ch", target: "or" },
			],
			solverSettings: settings(),
		});
		const { edges, loadNotice } = useStore.getState();
		expect(edges.map((e) => e.id)).toEqual(["f1"]);
		expect(loadNotice).toContain("w1");
		useStore.getState().dismissLoadNotice();
		expect(useStore.getState().loadNotice).toBeNull();
	});

	it("leaves a clean file without a notice", () => {
		useStore.getState().loadNetwork({
			nodes: [node("ch", "channel"), node("ch2", "channel")],
			edges: [
				{ id: "w1", source: "ch", target: "ch2", data: { type: "thermal" } },
			],
			solverSettings: settings(),
		});
		expect(useStore.getState().edges.map((e) => e.id)).toEqual(["w1"]);
		expect(useStore.getState().loadNotice).toBeNull();
	});
});

describe("the restore note (#596) goes with the restored canvas", () => {
	beforeEach(() => {
		useStore.getState().setRestoredAt("23/10/2026, 21:28");
	});

	it("is cleared by loading a file", () => {
		useStore
			.getState()
			.loadNetwork({ nodes: [nodeWithResult("n1")], edges: [] });
		expect(useStore.getState().restoredAt).toBeNull();
	});

	it("is cleared by Clear, with the canvas and its results", () => {
		useStore
			.getState()
			.loadNetwork({ nodes: [nodeWithResult("n1")], edges: [] });
		useStore.getState().setRestoredAt("23/10/2026, 21:28");
		useStore.getState().clearCanvas();
		const s = useStore.getState();
		expect(s.restoredAt).toBeNull();
		expect(s.nodes).toEqual([]);
		expect(s.edges).toEqual([]);
		expect(s.solveResults).toBeNull();
	});
});
