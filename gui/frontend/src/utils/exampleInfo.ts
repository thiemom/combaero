// Worked examples from papers (#793). A network file shipped under
// gui/examples/ carries an `example` key next to nodes/edges; the panel
// above the canvas shows it once the example loads.

export type ValidationLabel = "fidelity" | "cross-source";

export interface ExampleInfo {
	title: string;
	source: string;
	reproduced: string;
	whereToLook: string;
	label: ValidationLabel;
	labelNote: string;
	script: string | null;
}

export interface ExampleListing {
	id: string;
	title: string;
	label: string | null;
}

const text = (v: unknown): string => (typeof v === "string" ? v : "");

/** The description of a loaded file, or null when it carries none. Fields the
 *  file lacks come back empty rather than failing the load: the description is
 *  advisory, the network is what was asked for. An unknown label is null --
 *  the panel never shows a validation claim the file did not make. */
export const parseExampleInfo = (raw: unknown): ExampleInfo | null => {
	if (!raw || typeof raw !== "object") return null;
	const r = raw as Record<string, unknown>;
	const title = text(r.title);
	if (!title) return null;
	const label =
		r.label === "fidelity" || r.label === "cross-source" ? r.label : null;
	if (!label) return null;
	return {
		title,
		source: text(r.source),
		reproduced: text(r.reproduced),
		whereToLook: text(r.where_to_look),
		label,
		labelNote: text(r.label_note),
		script: text(r.script) || null,
	};
};

/** What the badge says, per docs/VALIDATION_POLICY.md: fidelity is judged on
 *  the author's own data, accuracy only cross-source. */
export const labelText = (label: ValidationLabel): string =>
	label === "fidelity"
		? "Fidelity: the source's own data"
		: "Cross-source: accuracy";
