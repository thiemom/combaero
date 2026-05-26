import type React from "react";
import { Position } from "reactflow";

// Clockwise rotation cycle matching CSS transform: rotate(Ndeg).
// Left→Top→Right→Bottom→Left for each 90° clockwise step.
const CYCLE = [Position.Left, Position.Top, Position.Right, Position.Bottom];

/**
 * Maps a handle's logical Position through a CSS clockwise rotation so the
 * edge tangent direction matches the handle's visual location on the node.
 *
 * Without this, ReactFlow re-measures the handle's screen coordinates via
 * updateNodeInternals (so the connection point is correct) but still uses
 * the original position prop for the spline tangent, causing edges to leave
 * at the wrong angle on rotated nodes.
 */
export const rotPos = (pos: Position, rotation: number): Position => {
	const steps = Math.round(rotation / 90) % 4;
	const idx = CYCLE.indexOf(pos);
	return CYCLE[(idx + steps + 4) % 4];
};

// Base positioning for each handle edge. Used by handleStyle below.
// The left/top/right/bottom values pin the handle to the correct physical
// edge within the node div's local (pre-rotation) coordinate system.
const HANDLE_BASE: Record<Position, React.CSSProperties> = {
	[Position.Left]: {
		left: "-4px",
		right: "auto",
		top: "50%",
		bottom: "auto",
	},
	[Position.Right]: {
		left: "auto",
		right: "-4px",
		top: "50%",
		bottom: "auto",
	},
	[Position.Top]: {
		left: "50%",
		right: "auto",
		top: "-4px",
		bottom: "auto",
	},
	[Position.Bottom]: {
		left: "50%",
		right: "auto",
		top: "auto",
		bottom: "-4px",
	},
};

// Centering offset for left/right vs top/bottom handles.
const HANDLE_CENTER: Record<Position, string> = {
	[Position.Left]: "translate(0, -50%)",
	[Position.Right]: "translate(0, -50%)",
	[Position.Top]: "translate(-50%, 0)",
	[Position.Bottom]: "translate(-50%, 0)",
};

/**
 * Returns the inline style for a Handle based on its pre-rotation base
 * position and the node's current rotation.
 *
 * Two effects:
 *  1. Overrides ReactFlow's position CSS classes (which would otherwise
 *     place the handle in the rotated coordinate system, landing it on the
 *     wrong visual edge).
 *  2. Counter-rotates the handle div by −rotation so that the ::after
 *     triangle drawn by index.css stays perpendicular to the edge it sits on,
 *     rather than inheriting the node's CSS rotation.
 */
export const handleStyle = (
	basePos: Position,
	rotation: number,
): React.CSSProperties => ({
	...HANDLE_BASE[basePos],
	transform: `${HANDLE_CENTER[basePos]} rotate(${-rotation}deg)`,
});
