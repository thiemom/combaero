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
