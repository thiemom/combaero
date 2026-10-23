// A thermal wall's box: auto-placed at its curve's midpoint, movable by the
// user. The offset is stored in flow coordinates (data.label_offset), so the
// box follows the wall when its nodes move.

export interface LabelOffset {
	dx: number;
	dy: number;
}

// Screen pixels the pointer must travel before a press on the box becomes a
// move, so a plain click stays a click.
export const DRAG_THRESHOLD_PX = 3;

export function labelPosition(
	labelX: number,
	labelY: number,
	offset: LabelOffset | undefined,
): [number, number] {
	return [labelX + (offset?.dx ?? 0), labelY + (offset?.dy ?? 0)];
}

// The offset after the pointer moved (sx, sy) screen pixels from where the
// press started, or null while still under the click threshold. Screen
// pixels become flow units through the zoom.
export function draggedOffset(
	start: LabelOffset,
	sx: number,
	sy: number,
	zoom: number,
): LabelOffset | null {
	if (Math.hypot(sx, sy) < DRAG_THRESHOLD_PX) return null;
	return { dx: start.dx + sx / zoom, dy: start.dy + sy / zoom };
}

type Side = "left" | "right" | "top" | "bottom";

// React Flow's bezier control offset (calculateControlOffset, curvature
// 0.25): half the distance when the handle faces the other point, a short
// loop out of the handle when it faces away.
function controlOffset(distance: number): number {
	return distance >= 0 ? 0.5 * distance : 0.25 * 25 * Math.sqrt(-distance);
}

function handleControl(
	x: number,
	y: number,
	side: Side,
	towardX: number,
	towardY: number,
): [number, number] {
	switch (side) {
		case "left":
			return [x - controlOffset(x - towardX), y];
		case "right":
			return [x + controlOffset(towardX - x), y];
		case "top":
			return [x, y - controlOffset(y - towardY)];
		default:
			return [x, y + controlOffset(towardY - y)];
	}
}

// A wall whose box was moved: one smooth curve from the source handle
// through the box to the target handle, so the arrow stays attached to the
// box. Each handle leaves the way React Flow's own bezier does; at the box
// the curve runs parallel to source -> target.
export function pathThroughBox(
	sx: number,
	sy: number,
	sourceSide: Side,
	tx: number,
	ty: number,
	targetSide: Side,
	bx: number,
	by: number,
): string {
	const [c1x, c1y] = handleControl(sx, sy, sourceSide, bx, by);
	const [c4x, c4y] = handleControl(tx, ty, targetSide, bx, by);
	const kx = (tx - sx) / 4;
	const ky = (ty - sy) / 4;
	return (
		`M ${sx},${sy} C ${c1x},${c1y} ${bx - kx},${by - ky} ${bx},${by}` +
		` C ${bx + kx},${by + ky} ${c4x},${c4y} ${tx},${ty}`
	);
}
