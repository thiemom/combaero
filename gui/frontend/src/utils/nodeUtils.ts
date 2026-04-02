import { Position } from "reactflow";

/**
 * Maps a base handle position through a clockwise rotation.
 *
 * At 0° (unrotated), the canonical positions are:
 *   Left = Inlet, Right = Outlet, Top/Bottom = Thermal
 *
 * Rotating 90° CW maps:  Left→Top, Top→Right, Right→Bottom, Bottom→Left
 */
export const getHandlePosition = (
	basePosition: Position,
	rotation: number,
): Position => {
	const normRotation = ((rotation % 360) + 360) % 360;
	if (normRotation === 0) return basePosition;

	const positions = [
		Position.Top,
		Position.Right,
		Position.Bottom,
		Position.Left,
	];
	const currentIndex = positions.indexOf(basePosition);
	const steps = normRotation / 90;
	const newIndex = (currentIndex + steps) % 4;

	return positions[newIndex];
};

/**
 * Returns the wrapper square size needed for a node with given width/height.
 * The square side = max(w, h) so the rotated inner box always fits.
 */
export const getWrapperSize = (w: number, h: number): number => Math.max(w, h);
