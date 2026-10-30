import type { PointerEvent as ReactPointerEvent } from "react";
import { useRef, useState } from "react";
import type { EdgeProps } from "reactflow";
import {
	EdgeLabelRenderer,
	getBezierPath,
	useStore as useFlowStore,
} from "reactflow";
import useStore from "../store/useStore";
import { FILM_COLOR, filmDuctId } from "../utils/filmLink";
import {
	draggedOffset,
	type LabelOffset,
	labelPosition,
	pathThroughBox,
} from "../utils/labelOffset";
import {
	orientationGlyph,
	orientationOf,
	solvedArrangement,
} from "../utils/wallArrangement";

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
	source,
	target,
	data,
	selected,
}: EdgeProps) {
	const { unitPreferences, updateEdgeData } = useStore();
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

	// Drag the box off whatever auto placement put it under; the offset is
	// saved with the wall (data.label_offset), double-click resets it.
	const zoom = useFlowStore((s) => s.transform[2]);
	const [dragOffset, setDragOffset] = useState<LabelOffset | null>(null);
	const press = useRef<{
		x: number;
		y: number;
		start: LabelOffset;
		last: LabelOffset | null;
	} | null>(null);
	const offset: LabelOffset | undefined = dragOffset ?? data?.label_offset;
	const [boxX, boxY] = labelPosition(labelX, labelY, offset);
	// A moved box keeps the arrow attached: node -> box -> node. Unmoved, the
	// wall is React Flow's own bezier, as before.
	const wallPath = offset
		? pathThroughBox(
				sourceX,
				sourceY,
				sourcePosition,
				targetX,
				targetY,
				targetPosition,
				boxX,
				boxY,
			)
		: edgePath;
	const onBoxPointerDown = (e: ReactPointerEvent<HTMLDivElement>) => {
		if (e.button !== 0) return;
		e.stopPropagation();
		e.currentTarget.setPointerCapture(e.pointerId);
		press.current = {
			x: e.clientX,
			y: e.clientY,
			start: data?.label_offset ?? { dx: 0, dy: 0 },
			last: null,
		};
	};
	const onBoxPointerMove = (e: ReactPointerEvent<HTMLDivElement>) => {
		const p = press.current;
		if (!p) return;
		const next = draggedOffset(p.start, e.clientX - p.x, e.clientY - p.y, zoom);
		if (!next) return;
		p.last = next;
		setDragOffset(next);
	};
	const onBoxPointerUp = (e: ReactPointerEvent<HTMLDivElement>) => {
		const p = press.current;
		press.current = null;
		if (e.currentTarget.hasPointerCapture(e.pointerId)) {
			e.currentTarget.releasePointerCapture(e.pointerId);
		}
		if (p?.last) {
			// A move is layout, not physics: keep the solved results.
			updateEdgeData(id, { label_offset: p.last }, true);
		}
		setDragOffset(null);
	};
	const onBoxDoubleClick = () => {
		if (data?.label_offset)
			updateEdgeData(id, { label_offset: undefined }, true);
	};

	// The film this wall holds (#386): a link from the holes to the wall's
	// label box, drawn only while the wall or the holes are selected so a
	// network with many films stays readable.
	const filmNode = useFlowStore((s) =>
		data?.film ? s.nodeInternals.get(data.film) : undefined,
	);
	const filmDuct = useFlowStore((s) =>
		data?.film
			? s.nodeInternals.get(filmDuctId({ source, target }, data.film_side))
			: undefined,
	);
	let filmLink: { d: string; x0: number; y0: number } | null = null;
	if (filmNode && (selected || filmNode.selected)) {
		// Positions are card centres (#865).
		const p = filmNode.positionAbsolute ?? filmNode.position;
		const x0 = p.x;
		const y0 = p.y;
		// It ends at the label point and runs under the label box, so it stops
		// visually at the box edge; a gentle bow keeps it off flow edges.
		const mx = (x0 + boxX) / 2 - (boxY - y0) * 0.2;
		const my = (y0 + boxY) / 2 + (boxX - x0) * 0.2;
		filmLink = { d: `M ${x0} ${y0} Q ${mx} ${my} ${boxX} ${boxY}`, x0, y0 };
	}

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
			{filmLink && (
				<g style={{ pointerEvents: "none" }} opacity={0.9}>
					<path
						d={filmLink.d}
						fill="none"
						stroke={FILM_COLOR}
						strokeWidth={1.25}
						strokeDasharray="2,3"
						strokeLinecap="round"
					/>
					<circle cx={filmLink.x0} cy={filmLink.y0} r={2.5} fill={FILM_COLOR} />
				</g>
			)}
			{/* Transparent wider path for hit-testing */}
			<path
				style={{ stroke: "transparent", strokeWidth: 10, fill: "none" }}
				className="react-flow__edge-path"
				d={wallPath}
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
				d={wallPath}
				markerEnd={flowsForward ? `url(#${markerId})` : undefined}
				markerStart={!flowsForward ? `url(#${markerIdReverse})` : undefined}
			/>
			<EdgeLabelRenderer>
				{/* biome-ignore lint/a11y/noStaticElementInteractions: a pointer drag of a canvas label, like a node drag; layout has no keyboard path in React Flow either */}
				<div
					style={{
						position: "absolute",
						transform: `translate(-50%, -50%) translate(${boxX}px,${boxY}px)`,
						background: "#fff",
						padding: "2px 6px",
						borderRadius: "4px",
						fontSize: "12px",
						color: ARROW_COLOR,
						border: `2px solid ${ARROW_COLOR}`,
						pointerEvents: "all",
						cursor: dragOffset ? "grabbing" : "grab",
						touchAction: "none",
					}}
					className="nodrag nopan flex flex-col items-center shadow-sm"
					data-film-wall={id}
					title="Drag to move; double-click to reset"
					onPointerDown={onBoxPointerDown}
					onPointerMove={onBoxPointerMove}
					onPointerUp={onBoxPointerUp}
					onPointerCancel={onBoxPointerUp}
					onDoubleClick={onBoxDoubleClick}
				>
					{data?.label && data?.show_label && (
						<span className="text-[9px] text-gray-400 font-bold uppercase tracking-wider leading-none mb-1">
							{data.label}
						</span>
					)}
					<span style={{ fontWeight: 700 }}>
						<span
							title={`${orientationOf(data)} wall${
								solvedArrangement(data?.result)
									? `, solved ${solvedArrangement(data?.result)} flow`
									: ""
							}`}
							style={{ marginRight: 4 }}
						>
							{orientationGlyph(orientationOf(data))}
						</span>
						{labelText}
					</span>
					{data?.film && filmNode && (
						<span
							className="text-[9px] leading-none font-bold"
							style={{ color: FILM_COLOR }}
							title={`Film from ${filmNode.data?.label || data.film} covers ${
								filmDuct?.data?.label || filmDuct?.id || "?"
							}`}
						>
							film
						</span>
					)}
					{solvedArrangement(data?.result) && (
						<span className="text-[9px] text-gray-500 leading-none">
							{solvedArrangement(data?.result)}
						</span>
					)}
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
