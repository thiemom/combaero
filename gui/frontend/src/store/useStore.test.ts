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
