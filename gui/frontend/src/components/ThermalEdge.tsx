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
	if (data?.result?.T_interface) {
		const temps = data.result.T_interface as number[];
		const targetX = data.probe_depth ?? 0;
		let currentX = 0;
		let probeTemp: number | null = null;

		// Fallback layers for legacy/initial elements
		const layers = data.layers || [
			{
				thickness: data.thickness || 0.003,
				conductivity: data.conductivity || 20.0,
			},
		];

		if (targetX <= 0) {
			probeTemp = temps[0];
		} else {
			let found = false;
			for (let i = 0; i < layers.length; i++) {
				const t = layers[i].thickness;
				const nextX = currentX + t;
				if (targetX <= nextX) {
					const frac = (targetX - currentX) / t;
					probeTemp = temps[i] + frac * (temps[i + 1] - temps[i]);
					found = true;
					break;
				}
				currentX = nextX;
			}
			if (!found) probeTemp = temps[temps.length - 1];
		}

		if (probeTemp !== null) {
			probeText = `t=${(targetX * scale).toFixed(unit === "mm" ? 1 : 3)}${unit}: ${probeTemp.toFixed(1)}K`;
		}
	}

	return (
		<>
			{/* Transparent path to carry the solid marker without inheriting the dash array */}
			<path
				style={{
					stroke: "transparent",
					strokeWidth: 3,
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
					strokeWidth: 3,
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
						padding: "2px 6px",
						borderRadius: "4px",
						fontSize: "12px",
						color: "#ff9800",
						border: "2px solid #ff9800",
						pointerEvents: "all",
					}}
					className="nodrag nopan flex flex-col items-center shadow-sm"
				>
					<span style={{ fontWeight: 700 }}>{labelText}</span>

					{probeText && (
						<span
							style={{ fontWeight: 400 }}
							className="text-[10px] text-gray-500 font-mono whitespace-nowrap"
						>
							{probeText}
						</span>
					)}
				</div>
			</EdgeLabelRenderer>
		</>
	);
}
