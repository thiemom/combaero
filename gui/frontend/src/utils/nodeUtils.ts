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

// Half the flow-handle size (12 px) so we can center without translate.
// Using calc() in top/left instead of transform:translate avoids the
// "translate in rotated axes" pitfall: if the transform only contains
// rotate(), rotating around 50%/50% is a pure in-place spin and
// getBoundingClientRect() returns the correct screen center both before
// and after rotation.
const H = "6px";

const HANDLE_BASE: Record<Position, React.CSSProperties> = {
	[Position.Left]: {
		left: "-4px",
		right: "auto",
		top: `calc(50% - ${H})`,
		bottom: "auto",
	},
	[Position.Right]: {
		left: "auto",
		right: "-4px",
		top: `calc(50% - ${H})`,
		bottom: "auto",
	},
	[Position.Top]: {
		left: `calc(50% - ${H})`,
		right: "auto",
		top: "-4px",
		bottom: "auto",
	},
	[Position.Bottom]: {
		left: `calc(50% - ${H})`,
		right: "auto",
		top: "auto",
		bottom: "-4px",
	},
};

/**
 * Inline style for a Handle.
 *
 * 1. Pins it to the correct physical edge (overriding ReactFlow's position
 *    CSS classes, which apply in the rotated coordinate system and would
 *    land the handle on the wrong visual edge).
 * 2. Counter-rotates the handle div by −rotation.  Because the node div has
 *    transform:rotate(+rotation), the net rotation of the handle and its
 *    ::after triangle is 0° relative to the screen — keeping the triangle
 *    perpendicular to the edge it sits on.
 *    The rotation is a pure spin (no translate), so getBoundingClientRect()
 *    returns the correct screen center at every rotation step.
 */
export const handleStyle = (
	basePos: Position,
	rotation: number,
): React.CSSProperties => ({
	...HANDLE_BASE[basePos],
	transform: `rotate(${-rotation}deg)`,
});
