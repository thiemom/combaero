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
 * A human label for an inherited-from note. Empty for a user-set or defaulted
 * value (which the UI styles differently), or when there is no source.
 */
export function inheritedSourceLabel(source: string | undefined): string {
	if (!source || source === "user" || source === "default") return "";
	return source; // e.g. "channel 'ch_u'" / "inherited from ch_u"
}
