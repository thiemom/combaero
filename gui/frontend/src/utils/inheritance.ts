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
 * The note an in-pipe orifice's Inspector shows about its bore against its
 * pipes (#788, #1039):
 * - "refused": the resolve pass refused the orifice, with the solver's own
 *   reason, whichever pipe caused it -- the upstream one or the downstream
 *   one that sets beta in reverse. A refusal against a momentum chamber
 *   (#1016) is left out: the orifice-family note shows it.
 * - "equal": the bore equals the upstream pipe (the solver's pass-through).
 * - "larger": the bore exceeds the upstream pipe, before the preview has
 *   said so (it then says "refused").
 * - null: nothing to say.
 */
export type OrificePipeNote =
	| { kind: "refused"; reason: string }
	| { kind: "equal" }
	| { kind: "larger" }
	| null;

export function orificePipeNote(args: {
	preview: PreviewItems | null | undefined;
	nodeId: string;
	inPipe: boolean;
	/** The upstream pipe's diameter from the preview, null when none. */
	pipeD: number | null;
	/** The pipe's pass-through verdict (not a chamber's). */
	passThrough: boolean;
	boreDiameter: number;
}): OrificePipeNote {
	const reason = unresolvedReason(args.preview, args.nodeId);
	if (reason && !reason.includes("momentum chamber"))
		return { kind: "refused", reason };
	if (!args.inPipe || args.pipeD == null) return null;
	if (args.passThrough) return { kind: "equal" };
	if (args.boreDiameter > args.pipeD) return { kind: "larger" };
	return null;
}

/**
 * A warning the solver will issue about a node's geometry without refusing it
 * (#1023: a momentum chamber narrower than an adjoining pipe). `undefined`
 * when there is none or the node is not in the preview.
 */
export function previewWarning(
	preview: PreviewItems | null | undefined,
	nodeId: string,
): string | undefined {
	const w = preview?.[nodeId]?.warning;
	return w ? w : undefined;
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
