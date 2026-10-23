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
