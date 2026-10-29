// A probe's link to its target (#856): a thin dotted, muted line from the
// probe to the node or thermal wall it reads. The probe's data.target_id is
// the one source of truth; the link is derived from it on every render and
// never stored, so saved files need no migration and the backend never sees
// it (graph_builder drops probes before the solve anyway).

import {
	type Edge,
	getBezierPath,
	internalsSymbol,
	type Node,
	type NodeChange,
	Position,
} from "reactflow";
import { labelPosition } from "./labelOffset";

/** The probe's hidden target handle. A connection drawn from a node to it
 *  points the probe at that node (the original way to attach a probe). */
export const PROBE_HANDLE = "probe-in";
/** The probe's link handle: drag from it onto a node or a thermal wall. */
export const PROBE_LINK_HANDLE = "probe-out";
/** The derived link's edge type (components/ProbeLinkEdge.tsx). */
export const PROBE_LINK_TYPE = "probe_link";
/** Muted stone, apart from flow (dark grey) and heat or thermal (orange/red). */
export const PROBE_LINK_COLOR = "#a8a29e";

const LINK_PREFIX = "probe-link:";

/** Canvas kinds that never enter the solve: adding, removing or moving one
 *  changes nothing the results solved. Only the probe: a heat source is an
 *  EnergyBoundary, and every other kind is a node or element (#856). */
export const DISPLAY_ONLY_KINDS: ReadonlySet<string> = new Set(["probe"]);

export const probeLinkId = (probeId: string): string =>
	`${LINK_PREFIX}${probeId}`;

/** The probe a link edge id belongs to, or null for any other edge. */
export function probeOfLink(edgeId: string): string | null {
	return edgeId.startsWith(LINK_PREFIX)
		? edgeId.slice(LINK_PREFIX.length)
		: null;
}

const isThermal = (e: Edge): boolean => e.data?.type === "thermal";

/** Whether a probe may read this node or edge: any node but a probe, or a
 *  thermal wall (the only edges a solve gives results for). */
export function isProbeTarget(
	probeId: string,
	targetId: string,
	nodes: Node[],
	edges: Edge[],
): boolean {
	if (targetId === probeId) return false;
	const node = nodes.find((n) => n.id === targetId);
	if (node) return !DISPLAY_ONLY_KINDS.has(node.type ?? "");
	const edge = edges.find((e) => e.id === targetId);
	return !!edge && isThermal(edge);
}

/** What a link edge carries: the target node, or the wall's ends and box
 *  offset so the link can find the wall's label point. */
export type ProbeLinkData =
	| { kind: "node"; targetId: string }
	| {
			kind: "edge";
			targetId: string;
			source: string;
			target: string;
			sourceHandle: string | null;
			targetHandle: string | null;
			labelOffset?: { dx: number; dy: number };
	  };

/**
 * One display-only edge per probe whose target is on the canvas. It loops on
 * the probe (link handle -> hidden probe-in handle) so React Flow needs no
 * handle on the target; ProbeLinkEdge draws the line to the target itself.
 * A probe whose target was deleted gets no link (its card says so).
 */
export function probeLinkEdges(
	nodes: Node[],
	edges: Edge[],
	selectedLinks: readonly string[] = [],
): Edge<ProbeLinkData>[] {
	const nodeIds = new Set(nodes.map((n) => n.id));
	const edgeById = new Map(edges.map((e) => [e.id, e]));
	const links: Edge<ProbeLinkData>[] = [];
	for (const p of nodes) {
		if (p.type !== "probe") continue;
		const targetId: string | undefined = p.data?.target_id;
		if (!targetId || targetId === p.id) continue;
		let data: ProbeLinkData | null = null;
		if (nodeIds.has(targetId)) {
			data = { kind: "node", targetId };
		} else {
			const wall = edgeById.get(targetId);
			if (wall) {
				data = {
					kind: "edge",
					targetId,
					source: wall.source,
					target: wall.target,
					sourceHandle: wall.sourceHandle ?? null,
					targetHandle: wall.targetHandle ?? null,
					labelOffset: isThermal(wall) ? wall.data?.label_offset : undefined,
				};
			}
		}
		if (!data) continue;
		const id = probeLinkId(p.id);
		links.push({
			id,
			type: PROBE_LINK_TYPE,
			source: p.id,
			sourceHandle: PROBE_LINK_HANDLE,
			target: p.id,
			targetHandle: PROBE_HANDLE,
			selected: selectedLinks.includes(id),
			data,
		});
	}
	return links;
}

/** A batch of node changes that alters the solved system: an add or remove
 *  of any kind but a display-only one. `before` is the node list the changes
 *  apply to (a removal names only the id). */
export function isStructuralNodeChange(
	changes: NodeChange[],
	before: Node[],
): boolean {
	const kindOf = new Map(before.map((n) => [n.id, n.type ?? ""]));
	return changes.some((c) => {
		if (c.type === "add") return !DISPLAY_ONLY_KINDS.has(c.item.type ?? "");
		if (c.type === "remove")
			return !DISPLAY_ONLY_KINDS.has(kindOf.get(c.id) ?? "");
		return false;
	});
}

/**
 * The canvas node or edge under a drop point, by id: a node's box, a thermal
 * wall's label box, or an edge's path. null over empty canvas.
 */
export function dropTargetId(el: Element | null): string | null {
	if (!el) return null;
	const node = el.closest(".react-flow__node")?.getAttribute("data-id");
	if (node) return node;
	const wallBox = el
		.closest("[data-film-wall]")
		?.getAttribute("data-film-wall");
	if (wallBox) return wallBox;
	const edge = el.closest(".react-flow__edge")?.getAttribute("data-testid");
	return edge?.startsWith("rf__edge-") ? edge.slice("rf__edge-".length) : null;
}

// --- Geometry: where the link ends ---

export interface Point {
	x: number;
	y: number;
}

export interface Rect {
	x: number;
	y: number;
	width: number;
	height: number;
}

/** A node's box in flow coordinates, as React Flow measured it. */
export function nodeRect(n: Node): Rect {
	const p = n.positionAbsolute ?? n.position;
	return { x: p.x, y: p.y, width: n.width ?? 140, height: n.height ?? 48 };
}

/** Where the segment from the box's centre toward `toward` leaves the box
 *  (the centre itself when `toward` is inside it). */
export function borderPoint(r: Rect, toward: Point): Point {
	const cx = r.x + r.width / 2;
	const cy = r.y + r.height / 2;
	const dx = toward.x - cx;
	const dy = toward.y - cy;
	const tx = dx === 0 ? Infinity : r.width / 2 / Math.abs(dx);
	const ty = dy === 0 ? Infinity : r.height / 2 / Math.abs(dy);
	const t = Math.min(tx, ty);
	if (!Number.isFinite(t) || t >= 1) return { x: cx, y: cy };
	return { x: cx + t * dx, y: cy + t * dy };
}

interface HandleBound {
	id?: string | null;
	position: Position;
	x: number;
	y: number;
	width: number;
	height: number;
}

// React Flow's own handle lookup and anchor point (getHandle and
// getHandlePosition in @reactflow/core), so the link meets the wall's label
// exactly where the wall's edge component put it.
function handleAnchor(
	n: Node,
	kind: "source" | "target",
	handleId: string | null,
): { point: Point; position: Position } | null {
	const bounds = (n as any)[internalsSymbol]?.handleBounds as
		| { source?: HandleBound[] | null; target?: HandleBound[] | null }
		| undefined;
	const list = bounds?.[kind];
	if (!list || list.length === 0) return null;
	const h =
		list.length === 1 || !handleId
			? list[0]
			: list.find((b) => b.id === handleId);
	if (!h) return null;
	const r = nodeRect(n);
	const x = r.x + h.x;
	const y = r.y + h.y;
	switch (h.position) {
		case Position.Top:
			return { point: { x: x + h.width / 2, y }, position: h.position };
		case Position.Right:
			return {
				point: { x: x + h.width, y: y + h.height / 2 },
				position: h.position,
			};
		case Position.Bottom:
			return {
				point: { x: x + h.width / 2, y: y + h.height },
				position: h.position,
			};
		default:
			return { point: { x, y: y + h.height / 2 }, position: Position.Left };
	}
}

/**
 * A wall's label point: the midpoint of its bezier (where FlowEdge and
 * ThermalEdge put the label) plus the box's saved offset. Without measured
 * handles, the midpoint of the two node centres.
 */
export function edgeLabelPoint(
	sourceNode: Node,
	targetNode: Node,
	sourceHandle: string | null,
	targetHandle: string | null,
	labelOffset?: { dx: number; dy: number },
): Point {
	const s = handleAnchor(sourceNode, "source", sourceHandle);
	const t = handleAnchor(targetNode, "target", targetHandle);
	let x: number;
	let y: number;
	if (s && t) {
		[, x, y] = getBezierPath({
			sourceX: s.point.x,
			sourceY: s.point.y,
			sourcePosition: s.position,
			targetX: t.point.x,
			targetY: t.point.y,
			targetPosition: t.position,
		});
	} else {
		const a = nodeRect(sourceNode);
		const b = nodeRect(targetNode);
		x = (a.x + a.width / 2 + b.x + b.width / 2) / 2;
		y = (a.y + a.height / 2 + b.y + b.height / 2) / 2;
	}
	const [lx, ly] = labelPosition(x, y, labelOffset);
	return { x: lx, y: ly };
}
