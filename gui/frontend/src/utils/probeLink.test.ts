/** @vitest-environment jsdom */
import { type Edge, internalsSymbol, type Node, Position } from "reactflow";
import { describe, expect, it } from "vitest";
import {
	borderPoint,
	dropTargetId,
	edgeLabelPoint,
	isProbeTarget,
	isStructuralNodeChange,
	PROBE_HANDLE,
	PROBE_LINK_HANDLE,
	PROBE_LINK_TYPE,
	probeLinkEdges,
	probeLinkId,
	probeOfLink,
} from "./probeLink";

const node = (id: string, type: string, data: any = {}): Node => ({
	id,
	type,
	position: { x: 0, y: 0 },
	data,
});

const wall: Edge = {
	id: "w",
	source: "a",
	target: "b",
	sourceHandle: "a-thermal",
	targetHandle: "b-thermal",
	type: "thermal",
	data: { type: "thermal", label_offset: { dx: 5, dy: -5 } },
};
const flow: Edge = { id: "f", source: "a", target: "b", data: {} };

describe("probeLinkEdges: the link is derived from target_id", () => {
	it("draws one display-only link per probe with a target on the canvas", () => {
		const nodes = [
			node("a", "channel"),
			node("b", "plenum"),
			node("p", "probe", { target_id: "a" }),
			node("q", "probe", { target_id: "w" }),
			node("r", "probe", {}),
		];
		const links = probeLinkEdges(nodes, [wall, flow]);
		expect(links.map((l) => l.id)).toEqual([
			probeLinkId("p"),
			probeLinkId("q"),
		]);
		const [toNode, toWall] = links;
		expect(toNode).toMatchObject({
			type: PROBE_LINK_TYPE,
			source: "p",
			sourceHandle: PROBE_LINK_HANDLE,
			target: "p",
			targetHandle: PROBE_HANDLE,
			data: { kind: "node", targetId: "a" },
		});
		expect(toWall.data).toEqual({
			kind: "edge",
			targetId: "w",
			source: "a",
			target: "b",
			sourceHandle: "a-thermal",
			targetHandle: "b-thermal",
			labelOffset: { dx: 5, dy: -5 },
		});
	});

	it("moves with target_id and disappears when it is cleared or dangling", () => {
		const nodes = [
			node("a", "channel"),
			node("b", "plenum"),
			node("p", "probe", { target_id: "b" }),
		];
		expect(probeLinkEdges(nodes, [])[0].data?.targetId).toBe("b");
		nodes[2] = node("p", "probe", {});
		expect(probeLinkEdges(nodes, [])).toEqual([]);
		// A deleted target: no link, the card says "target deleted".
		nodes[2] = node("p", "probe", { target_id: "gone" });
		expect(probeLinkEdges(nodes, [])).toEqual([]);
	});

	it("marks the selected links", () => {
		const nodes = [
			node("a", "channel"),
			node("p", "probe", { target_id: "a" }),
		];
		expect(probeLinkEdges(nodes, [])[0].selected).toBe(false);
		expect(probeLinkEdges(nodes, [], [probeLinkId("p")])[0].selected).toBe(
			true,
		);
	});

	it("link ids round-trip and no other edge parses as one", () => {
		expect(probeOfLink(probeLinkId("node_1"))).toBe("node_1");
		expect(probeOfLink("edge_1")).toBeNull();
	});
});

describe("isProbeTarget", () => {
	const nodes = [node("a", "channel"), node("p", "probe"), node("q", "probe")];
	it("takes any node but a probe, and a thermal wall", () => {
		expect(isProbeTarget("p", "a", nodes, [wall, flow])).toBe(true);
		expect(isProbeTarget("p", "w", nodes, [wall, flow])).toBe(true);
	});
	it("refuses itself, another probe, a flow edge and a missing id", () => {
		expect(isProbeTarget("p", "p", nodes, [wall, flow])).toBe(false);
		expect(isProbeTarget("p", "q", nodes, [wall, flow])).toBe(false);
		expect(isProbeTarget("p", "f", nodes, [wall, flow])).toBe(false);
		expect(isProbeTarget("p", "nope", nodes, [wall, flow])).toBe(false);
	});
});

describe("isStructuralNodeChange", () => {
	const before = [node("a", "channel"), node("p", "probe")];
	it("a probe added or removed is not structural", () => {
		expect(
			isStructuralNodeChange(
				[{ type: "add", item: node("p2", "probe") }],
				before,
			),
		).toBe(false);
		expect(isStructuralNodeChange([{ type: "remove", id: "p" }], before)).toBe(
			false,
		);
	});
	it("any other kind added or removed is", () => {
		expect(
			isStructuralNodeChange(
				[{ type: "add", item: node("h", "heat_source") }],
				before,
			),
		).toBe(true);
		expect(
			isStructuralNodeChange(
				[
					{ type: "remove", id: "p" },
					{ type: "remove", id: "a" },
				],
				before,
			),
		).toBe(true);
	});
});

describe("dropTargetId", () => {
	it("finds a node, a wall's label box or an edge path under the pointer", () => {
		document.body.innerHTML = `
			<div class="react-flow__node" data-id="n1"><span id="in-node">x</span></div>
			<div data-film-wall="w1"><span id="in-box">Q</span></div>
			<svg><g class="react-flow__edge" data-testid="rf__edge-w2"><path id="in-edge"/></g></svg>
			<div id="pane"></div>`;
		expect(dropTargetId(document.getElementById("in-node"))).toBe("n1");
		expect(dropTargetId(document.getElementById("in-box"))).toBe("w1");
		expect(dropTargetId(document.getElementById("in-edge"))).toBe("w2");
		expect(dropTargetId(document.getElementById("pane"))).toBeNull();
		expect(dropTargetId(null)).toBeNull();
	});
});

describe("link geometry", () => {
	const box = { x: 0, y: 0, width: 100, height: 40 };
	it("leaves a node's box toward the other end", () => {
		expect(borderPoint(box, { x: 250, y: 20 })).toEqual({ x: 100, y: 20 });
		expect(borderPoint(box, { x: 50, y: -100 })).toEqual({ x: 50, y: 0 });
		// An end inside the box: its centre.
		expect(borderPoint(box, { x: 60, y: 25 })).toEqual({ x: 50, y: 20 });
	});

	it("ends on a wall at its label point, offset included", () => {
		const measured = (id: string, x: number, pos: Position, kind: string) =>
			({
				...node(id, "plenum"),
				// Positions are card centres (#865): this box spans x .. x + 100.
				position: { x: x + 50, y: 20 },
				positionAbsolute: { x: x + 50, y: 20 },
				width: 100,
				height: 40,
				[internalsSymbol]: {
					handleBounds: {
						[kind]: [
							{
								id: `${id}-h`,
								position: pos,
								x: 0,
								y: 15,
								width: 10,
								height: 10,
							},
						],
					},
				},
			}) as unknown as Node;
		// Source handle's anchor at x = 0 + 10 (Right edge of handle), target's
		// at x = 300 (Left): the straight bezier's midpoint is (155, 20).
		const a = measured("a", 0, Position.Right, "source");
		const b = measured("b", 300, Position.Left, "target");
		expect(edgeLabelPoint(a, b, "a-h", "b-h")).toEqual({ x: 155, y: 20 });
		expect(edgeLabelPoint(a, b, "a-h", "b-h", { dx: 5, dy: -5 })).toEqual({
			x: 160,
			y: 15,
		});
	});

	it("falls back to the node centres before React Flow measured handles", () => {
		// Card centres (#865): boxes at (0, 0) and (200, 100), 100 x 40.
		const a = {
			...node("a", "plenum"),
			position: { x: 50, y: 20 },
			width: 100,
			height: 40,
		};
		const b = {
			...node("b", "plenum"),
			position: { x: 250, y: 120 },
			width: 100,
			height: 40,
		};
		expect(edgeLabelPoint(a, b, null, null)).toEqual({ x: 150, y: 70 });
	});
});

describe("a fuel boundary's link to the combustor it controls (#892)", () => {
	it("is derived from combustor_id and marked as a control link", () => {
		const nodes = [
			node("cb", "combustor"),
			node("f", "fuel_boundary", { combustor_id: "cb" }),
			node("g", "fuel_boundary", {}),
			node("h", "fuel_boundary", { combustor_id: "gone" }),
		];
		const links = probeLinkEdges(nodes, []);
		expect(links.map((l) => l.id)).toEqual([probeLinkId("f")]);
		expect(links[0].data).toEqual({
			kind: "node",
			targetId: "cb",
			control: true,
		});
		expect(probeOfLink(links[0].id)).toBe("f");
	});

	it("may point at a combustor only", () => {
		const nodes = [
			node("cb", "combustor"),
			node("pl", "plenum"),
			node("f", "fuel_boundary"),
		];
		expect(isProbeTarget("f", "cb", nodes, [wall])).toBe(true);
		expect(isProbeTarget("f", "pl", nodes, [wall])).toBe(false);
		expect(isProbeTarget("f", "w", nodes, [wall])).toBe(false);
		expect(isProbeTarget("f", "f", nodes, [wall])).toBe(false);
	});

	it("is structural: adding a fuel boundary changes the solved system", () => {
		expect(
			isStructuralNodeChange(
				[{ type: "add", item: node("f", "fuel_boundary") }],
				[],
			),
		).toBe(true);
	});
});
