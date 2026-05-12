import type { EdgeProps } from "reactflow";
import {
	BaseEdge,
	EdgeLabelRenderer,
	getBezierPath,
	Position,
	useNodes,
} from "reactflow";

function rotatePosition(pos: Position, degrees: number): Position {
	const cycle: Position[] = [
		Position.Right,
		Position.Bottom,
		Position.Left,
		Position.Top,
	];
	const steps = Math.round((((degrees % 360) + 360) % 360) / 90) % 4;
	const idx = cycle.indexOf(pos);
	return idx === -1 ? pos : cycle[(idx + steps) % 4];
}

export default function FlowEdge({
	id,
	source,
	target,
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
	const nodes = useNodes();
	const sourceRotation =
		(nodes.find((n) => n.id === source)?.data as { rotation?: number })
			?.rotation ?? 0;
	const targetRotation =
		(nodes.find((n) => n.id === target)?.data as { rotation?: number })
			?.rotation ?? 0;

	const adjSourcePos = rotatePosition(sourcePosition, sourceRotation);
	const adjTargetPos = rotatePosition(targetPosition, targetRotation);

	const [edgePath, labelX, labelY] = getBezierPath({
		sourceX,
		sourceY,
		sourcePosition: adjSourcePos,
		targetX,
		targetY,
		targetPosition: adjTargetPos,
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
