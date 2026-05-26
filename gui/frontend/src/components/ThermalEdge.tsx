import type { EdgeProps } from "reactflow";
import { EdgeLabelRenderer, getBezierPath } from "reactflow";
import useStore from "../store/useStore";

const ARROW_COLOR = "#ff9800";
const MARKER_SIZE = 6;

// Hot end: existing orange. Cool end: sky blue.
// t=0 → hot, t=1 → cold.
function thermalColor(t: number): string {
	const r = Math.round(255 + (41 - 255) * t);
	const g = Math.round(152 + (182 - 152) * t);
	const b = Math.round(0 + (246 - 0) * t);
	return `rgb(${r},${g},${b})`;
}

export default function ThermalEdge({
	id,
	sourceX,
	sourceY,
	targetX,
	targetY,
	sourcePosition,
	targetPosition,
	style = {},
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

	const Q: number | undefined = data?.result?.Q;
	const flowsForward = Q === undefined || Q >= 0;

	const topoLayers: Array<{ thickness?: number }> = data?.layers || [
		{ thickness: data?.thickness || 0.003 },
	];
	const nLayers = topoLayers.length;
	const totalT =
		topoLayers.reduce((s, l) => s + (l.thickness || 0.003), 0) || 0.003;

	const markerId = `thermal-arrow-${id}`;
	const markerIdReverse = `thermal-arrow-reverse-${id}`;

	let labelText = "Thermal";
	let probeText = "";
	if (Q !== undefined) {
		const absQ = Math.abs(Q);
		const formattedValue =
			absQ >= 1000 ? `${(Q / 1000).toFixed(1)} kW` : `${Q.toFixed(1)} W`;
		labelText = formattedValue;
	}
	if (data?.result?.T_interface) {
		const temps = data.result.T_interface as number[];
		const layers = data.layers || [
			{
				thickness: data.thickness || 0.003,
				conductivity: data.conductivity || 20.0,
			},
		];

		let probeDepth = data.probe_depth ?? 0;
		let displayLabel = `d=${(probeDepth * scale).toFixed(unit === "mm" ? 1 : 3)}${unit}`;
		let hotSideShortcut = false;

		if (data.probe_mode === "preset" && data.probe_preset) {
			const { type, index } = data.probe_preset;
			const safeIndex = Math.min(index, layers.length - 1);
			const layer = layers[safeIndex];

			let runningX = 0;
			for (let i = 0; i < safeIndex; i++) {
				runningX += layers[i].thickness;
			}

			const L = layer?.thickness || 0;
			if (type === "hot") {
				probeDepth = flowsForward ? runningX : runningX + L;
				displayLabel = `L${safeIndex + 1} ${flowsForward ? "hot" : "cold"}`;
			} else if (type === "avg") {
				probeDepth = runningX + L / 2;
				displayLabel = `L${safeIndex + 1} avg`;
			} else {
				probeDepth = flowsForward ? runningX + L : runningX;
				displayLabel = `L${safeIndex + 1} ${flowsForward ? "cold" : "hot"}`;
			}
		} else if (!data.probe_mode) {
			hotSideShortcut = true;
			displayLabel = "hot side";
		}

		if (temps.length !== layers.length + 1) {
			probeText = `${displayLabel}: (solve needed)`;
		} else {
			let probeTemp: number | null = null;
			if (hotSideShortcut) {
				probeTemp = flowsForward ? temps[0] : temps[temps.length - 1];
			} else if (probeDepth <= 0) {
				probeTemp = temps[0];
			} else {
				let found = false;
				let iterX = 0;
				for (let i = 0; i < layers.length; i++) {
					const t = layers[i].thickness;
					const nextX = iterX + t;
					if (probeDepth <= nextX + 1e-9) {
						const frac = t > 0 ? (probeDepth - iterX) / t : 0;
						probeTemp = temps[i] + frac * (temps[i + 1] - temps[i]);
						found = true;
						break;
					}
					iterX = nextX;
				}
				if (!found) probeTemp = temps[temps.length - 1];
			}

			if (probeTemp !== null) {
				probeText = `${displayLabel}: ${probeTemp.toFixed(1)}K`;
			}
		}
	}

	return (
		<>
			<defs>
				<marker
					id={markerId}
					markerWidth={MARKER_SIZE}
					markerHeight={MARKER_SIZE}
					refX={MARKER_SIZE}
					refY={MARKER_SIZE / 2}
					orient="auto"
				>
					<path
						d={`M 0 0 L ${MARKER_SIZE} ${MARKER_SIZE / 2} L 0 ${MARKER_SIZE} z`}
						fill={ARROW_COLOR}
					/>
				</marker>
				<marker
					id={markerIdReverse}
					markerWidth={MARKER_SIZE}
					markerHeight={MARKER_SIZE}
					refX={0}
					refY={MARKER_SIZE / 2}
					orient="auto-start-reverse"
				>
					<path
						d={`M 0 0 L ${MARKER_SIZE} ${MARKER_SIZE / 2} L 0 ${MARKER_SIZE} z`}
						fill={ARROW_COLOR}
					/>
				</marker>
			</defs>
			{/* Transparent wider path for hit-testing */}
			<path
				style={{ stroke: "transparent", strokeWidth: 10, fill: "none" }}
				className="react-flow__edge-path"
				d={edgePath}
			/>
			{/* Visible dashed path with data-driven arrowhead */}
			<path
				id={id}
				style={{
					...style,
					strokeWidth: 3,
					stroke: ARROW_COLOR,
					strokeDasharray: "5,5",
					opacity: 0.8,
					fill: "none",
				}}
				className="react-flow__edge-path"
				d={edgePath}
				markerEnd={flowsForward ? `url(#${markerId})` : undefined}
				markerStart={!flowsForward ? `url(#${markerIdReverse})` : undefined}
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
						color: ARROW_COLOR,
						border: `2px solid ${ARROW_COLOR}`,
						pointerEvents: "all",
					}}
					className="nodrag nopan flex flex-col items-center shadow-sm"
				>
					{data?.label && data?.show_label && (
						<span className="text-[9px] text-gray-400 font-bold uppercase tracking-wider leading-none mb-1">
							{data.label}
						</span>
					)}
					<span style={{ fontWeight: 700 }}>{labelText}</span>
					<div
						style={{
							display: "flex",
							alignItems: "center",
							gap: 3,
							marginTop: 2,
							width: "100%",
						}}
					>
						<span
							style={{
								fontSize: 8,
								fontWeight: 700,
								color: ARROW_COLOR,
								lineHeight: 1,
								flexShrink: 0,
							}}
						>
							L1
						</span>
						<div
							style={{
								display: "flex",
								gap: 1,
								height: 6,
								flex: 1,
								minWidth: 32,
								borderRadius: 2,
								overflow: "hidden",
							}}
						>
							{topoLayers.map((layer, i) => {
								const topo = nLayers <= 1 ? 0 : i / (nLayers - 1);
								const t = flowsForward ? topo : 1 - topo;
								return (
									<div
										key={i}
										style={{
											flex: (layer.thickness || 0.003) / totalT,
											minWidth: 3,
											background: thermalColor(t),
										}}
									/>
								);
							})}
						</div>
						<span
							style={{
								fontSize: 8,
								fontWeight: 700,
								color: ARROW_COLOR,
								lineHeight: 1,
								flexShrink: 0,
							}}
						>
							L{nLayers}
						</span>
					</div>
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
