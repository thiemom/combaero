import { Position } from "reactflow";

/**
 * Returns the actual React Flow Position based on the node's rotation.
 * This ensures that a handle defined as "Left" stays on the visually left side
 * even when the content inside is rotated.
 *
 * However, the user wants the handles to ROTATE WITH the content.
 * So if we rotate 90 deg, the "Left" handle (target) should move to "Top".
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

	// Rotate clockwise: 90 deg = +1 index
	const steps = normRotation / 90;
	const newIndex = (currentIndex + steps) % 4;

	return positions[newIndex];
};
