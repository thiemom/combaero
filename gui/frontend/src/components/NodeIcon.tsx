import type React from "react";
import {
	type AppearanceType,
	iconClass,
	NODE_APPEARANCE,
} from "../utils/nodeAppearance";

/** A node type's glyph in its colour, the same in the palette and on the canvas (#321). */
const NodeIcon: React.FC<{
	type: AppearanceType;
	size?: number;
	className?: string;
}> = ({ type, size = 16, className = "" }) => {
	const Icon = NODE_APPEARANCE[type].icon;
	return <Icon size={size} className={iconClass(type, className)} />;
};

export default NodeIcon;
