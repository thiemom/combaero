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

// Inline styles that pin each Handle to its correct physical edge within the
// rotated node div. ReactFlow positions handles via CSS classes
// (.react-flow__handle-left etc.) but since the node div has a CSS transform
// those classes apply in the rotated coordinate system — causing the handle
// to appear on the wrong visual edge. Inline styles (higher specificity) keep
// the handle at the intended pre-rotation location regardless of rotation.
export const HANDLE_CSS: Record<Position, React.CSSProperties> = {
	[Position.Left]: {
		left: "-4px",
		right: "auto",
		top: "50%",
		bottom: "auto",
		transform: "translate(0, -50%)",
	},
	[Position.Right]: {
		left: "auto",
		right: "-4px",
		top: "50%",
		bottom: "auto",
		transform: "translate(0, -50%)",
	},
	[Position.Top]: {
		left: "50%",
		right: "auto",
		top: "-4px",
		bottom: "auto",
		transform: "translate(-50%, 0)",
	},
	[Position.Bottom]: {
		left: "50%",
		right: "auto",
		top: "auto",
		bottom: "-4px",
		transform: "translate(-50%, 0)",
	},
};
