import type { EdgeProps } from "reactflow";
import { EdgeLabelRenderer, getBezierPath } from "reactflow";

export default function ThermalEdge({
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

	let labelText = "Thermal";
	let profileText = "";
	if (data?.result?.Q !== undefined) {
		const Q = data.result.Q;
		const absQ = Math.abs(Q);
		const formattedValue =
			absQ >= 1000 ? `${(Q / 1000).toFixed(1)} kW` : `${Q.toFixed(1)} W`;
		labelText = formattedValue;
	}
	if (data?.result?.T_interface) {
		const temps = data.result.T_interface as number[];
		profileText = temps.map((t) => `${Math.round(t)}K`).join(" → ");
	}

	return (
		<>
			{/* Transparent path to carry the solid marker without inheriting the dash array */}
			<path
				style={{
					stroke: "transparent",
					strokeWidth: 4,
					fill: "none",
					pointerEvents: "none",
				}}
				className="react-flow__edge-path"
				d={edgePath}
				markerEnd={markerEnd}
			/>
			{/* Visible dashed path for the connection line */}
			<path
				id={id}
				style={{
					...style,
					strokeWidth: 4,
					stroke: "#ff9800",
					strokeDasharray: "5,5",
					opacity: 0.8,
					fill: "none",
				}}
				className="react-flow__edge-path"
				d={edgePath}
			/>
			<EdgeLabelRenderer>
				<div
					style={{
						position: "absolute",
						transform: `translate(-50%, -50%) translate(${labelX}px,${labelY}px)`,
						background: "#fff",
						padding: "2px 4px",
						borderRadius: "4px",
						fontSize: "10px",
						fontWeight: 700,
						color: "#ff9800",
						border: "1px solid #ff9800",
						pointerEvents: "all",
					}}
					className="nodrag nopan flex flex-col items-center"
				>
					<span>{labelText}</span>
					{profileText && (
						<span className="text-[8px] text-stone-500 font-mono tracking-tighter mt-0.5">
							{profileText}
						</span>
					)}
				</div>
			</EdgeLabelRenderer>
		</>
	);
}
