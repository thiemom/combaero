import type { EdgeProps } from "reactflow";
import { BaseEdge, useStore as useFlowStore } from "reactflow";
import {
	borderPoint,
	edgeLabelPoint,
	nodeRect,
	type Point,
	PROBE_LINK_COLOR,
	type ProbeLinkData,
} from "../utils/probeLink";

const SELECTED_COLOR = "#57534e";

/**
 * A probe's link to its target (#856), derived from the probe's target_id
 * (utils/probeLink.ts): thin, dotted and muted, so it never reads as flow or
 * as a heat source's dashed link. It runs border to border between the probe
 * and a target node, or to a thermal wall's label box. Selecting it and
 * pressing Delete clears the probe's target.
 */
export default function ProbeLinkEdge({
	id,
	source,
	data,
	selected,
}: EdgeProps<ProbeLinkData>) {
	const probe = useFlowStore((s) => s.nodeInternals.get(source));
	// A node target, or a wall's two ends.
	const a = useFlowStore((s) =>
		data
			? s.nodeInternals.get(data.kind === "edge" ? data.source : data.targetId)
			: undefined,
	);
	const b = useFlowStore((s) =>
		data?.kind === "edge" ? s.nodeInternals.get(data.target) : undefined,
	);
	if (!data || !probe || !a) return null;

	const probeBox = nodeRect(probe);
	const probeCentre: Point = {
		x: probeBox.x + probeBox.width / 2,
		y: probeBox.y + probeBox.height / 2,
	};
	let end: Point;
	if (data.kind === "edge") {
		if (!b) return null;
		end = edgeLabelPoint(
			a,
			b,
			data.sourceHandle,
			data.targetHandle,
			data.labelOffset,
		);
	} else {
		end = borderPoint(nodeRect(a), probeCentre);
	}
	const start = borderPoint(probeBox, end);
	const color = selected ? SELECTED_COLOR : PROBE_LINK_COLOR;

	return (
		<>
			<BaseEdge
				id={id}
				path={`M ${start.x},${start.y} L ${end.x},${end.y}`}
				interactionWidth={12}
				style={{
					stroke: color,
					strokeWidth: selected ? 1.5 : 1,
					strokeDasharray: "1 3",
					strokeLinecap: "round",
				}}
			/>
			<circle
				cx={end.x}
				cy={end.y}
				r={2}
				fill={color}
				style={{ pointerEvents: "none" }}
			/>
		</>
	);
}
