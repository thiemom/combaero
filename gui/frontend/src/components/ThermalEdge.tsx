import type { EdgeProps } from "reactflow";
import { EdgeLabelRenderer, getBezierPath } from "reactflow";
import useStore from "../store/useStore";

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
	const { unitPreferences } = useStore();
	const unit = unitPreferences.length;
	const scale = unit === "mm" ? 1000 : 1;

	const [edgePath, labelX, labelY] = getBezierPath({
		sourceX,
		sourceY,
		sourcePosition,
		targetX,
		targetY,
		targetPosition,
	});

	let labelText = "Thermal";
	let probeText = "";
	if (data?.result?.Q !== undefined) {
		const Q = data.result.Q;
		const absQ = Math.abs(Q);
		const formattedValue =
			absQ >= 1000 ? `${(Q / 1000).toFixed(1)} kW` : `${Q.toFixed(1)} W`;
		labelText = formattedValue;
	}
	if (data?.result?.T_interface && data.layers) {
		const temps = data.result.T_interface as number[];
		const targetX = data.probe_depth ?? 0;
		let currentX = 0;
		let pTemp: number | null = null;

		if (targetX <= 0) {
			pTemp = temps[0];
		} else {
			let found = false;
			for (let i = 0; i < data.layers.length; i++) {
				const t = data.layers[i].thickness;
				const nextX = currentX + t;
				if (targetX <= nextX) {
					const frac = (targetX - currentX) / t;
					pTemp = temps[i] + frac * (temps[i + 1] - temps[i]);
					found = true;
					break;
				}
				currentX = nextX;
			}
			if (!found) pTemp = temps[temps.length - 1];
		}

		if (pTemp !== null) {
			probeText = `@${(targetX * scale).toFixed(unit === "mm" ? 1 : 3)}${unit}: ${pTemp.toFixed(1)}K`;
		}
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

					{probeText && (
						<span className="text-[9px] text-orange-600 border border-orange-200 bg-orange-50 font-bold px-1 rounded tracking-tighter mt-0.5">
							{probeText}
						</span>
					)}
				</div>
			</EdgeLabelRenderer>
		</>
	);
}
