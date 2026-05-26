import type { EdgeProps } from "reactflow";
import { BaseEdge, EdgeLabelRenderer, getBezierPath } from "reactflow";

export default function FlowEdge({
	id,
	sourceX,
	sourceY,
	targetX,
	targetY,
	sourcePosition,
	targetPosition,
	style = {},
	markerEnd,
	data,
}: EdgeProps) {
	const [edgePath, labelX, labelY] = getBezierPath({
		sourceX,
		sourceY,
		sourcePosition,
		targetX,
		targetY,
		targetPosition,
	});

	const showLabel = data?.show_label && data?.label;

	return (
		<>
			<BaseEdge id={id} path={edgePath} markerEnd={markerEnd} style={style} />
			{showLabel && (
				<EdgeLabelRenderer>
					<div
						style={{
							position: "absolute",
							transform: `translate(-50%, -50%) translate(${labelX}px,${labelY}px)`,
							background: "#fff",
							padding: "2px 6px",
							borderRadius: "4px",
							fontSize: "11px",
							color: "#475569",
							border: "1.5px solid #cbd5e1",
							pointerEvents: "none",
							fontWeight: 600,
						}}
						className="nodrag nopan shadow-sm"
					>
						{data.label}
					</div>
				</EdgeLabelRenderer>
			)}
		</>
	);
}
