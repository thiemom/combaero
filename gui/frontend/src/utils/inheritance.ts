import type { PreviewItem } from "../api";

// The authoritative geometry-inheritance view, from the backend /resolve_preview
// endpoint (#559), keyed by canvas node id. The GUI reads this instead of
// re-deriving inheritance in its own graph walk, so the placeholder it shows and
// the value the solver uses are the same number.
export type PreviewItems = Record<string, PreviewItem>;

export interface ResolvedField {
	value: number | null;
	// "user" | "inherited from <id>" | "channel '<id>'" | "default" | ""
	source: string;
}

/**
 * The resolved value and provenance of one geometry field of a node, from the
 * preview. `undefined` when the preview has no entry yet (not fetched, or the
 * node is not a geometry item).
 */
export function inheritedField(
	preview: PreviewItems | null | undefined,
	nodeId: string,
	field: string,
): ResolvedField | undefined {
	const f = preview?.[nodeId]?.fields?.[field];
	return f ? { value: f.value, source: f.source } : undefined;
}

/**
 * When the solver cannot resolve a node's geometry (it would refuse), the reason
 * -- so the GUI can flag it before Solve rather than as a bare failure.
 * `undefined` when the node is resolved or not in the preview.
 */
export function unresolvedReason(
	preview: PreviewItems | null | undefined,
	nodeId: string,
): string | undefined {
	const item = preview?.[nodeId];
	return item && item.status === "unresolved" ? item.reason : undefined;
}

/**
 * The note shown under an inherited value, as a complete phrase. Sources come
 * in three shapes -- "channel 'x'" (or "<Element> 'x'", "node 'x'"),
 * "inherited from x", and "derived from ..." -- so the phrase is built per
 * shape rather than by a fixed "Inherited from {source}" template, which would
 * read "Inherited from inherited from x".
 *
 * Empty for a user-set or defaulted value (styled differently), or no source.
 */
export function inheritedNote(source: string | undefined): string {
	if (!source || source === "user" || source === "default") return "";
	if (
		source.startsWith("inherited from ") ||
		source.startsWith("derived from ")
	) {
		return source.charAt(0).toUpperCase() + source.slice(1);
	}
	return `Inherited from ${source}`;
}
