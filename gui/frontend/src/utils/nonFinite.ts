/**
 * Non-finite numbers in a network file (#1038).
 *
 * JSON has no Infinity or NaN, but `JSON.parse` turns an overflowing literal
 * such as `1e999` into Infinity, and `JSON.stringify` writes Infinity and NaN
 * as `null`. A non-finite value that reaches the canvas is therefore saved as
 * null, which the backend refuses for a number ("Input should be a valid
 * number"). Load refuses such a file and Save refuses to write one, both
 * naming the offending fields.
 *
 * Solve results (`data.result`) are skipped: load discards them and the
 * backend never reads them.
 */

/** One non-finite value: where it is and what it was. */
export interface NonFiniteField {
	/** "node 'comp'", "edge 'e1'" or "solver settings". */
	owner: string;
	/** Dotted path inside the owner's data, e.g. "pressure_ratio". */
	field: string;
	value: number;
}

function walk(
	value: unknown,
	path: string,
	owner: string,
	out: NonFiniteField[],
): void {
	if (typeof value === "number") {
		if (!Number.isFinite(value)) out.push({ owner, field: path, value });
		return;
	}
	if (Array.isArray(value)) {
		value.forEach((v, i) => {
			walk(v, `${path}[${i}]`, owner, out);
		});
		return;
	}
	if (value !== null && typeof value === "object") {
		for (const [k, v] of Object.entries(value)) {
			if (k === "result") continue;
			walk(v, path ? `${path}.${k}` : k, owner, out);
		}
	}
}

/** Every non-finite number in the nodes', edges' and solver settings' data. */
export function findNonFinite(network: {
	nodes?: unknown;
	edges?: unknown;
	solverSettings?: unknown;
}): NonFiniteField[] {
	const out: NonFiniteField[] = [];
	const items = (list: unknown, kind: string) => {
		if (!Array.isArray(list)) return;
		for (const item of list) {
			if (item === null || typeof item !== "object") continue;
			const { id, data } = item as { id?: unknown; data?: unknown };
			walk(data, "", `${kind} '${String(id)}'`, out);
		}
	};
	items(network.nodes, "node");
	items(network.edges, "edge");
	walk(network.solverSettings, "", "solver settings", out);
	return out;
}

/** "node 'comp' pressure_ratio = Infinity; ..." for a notice. */
export function describeNonFinite(fields: NonFiniteField[]): string {
	return fields
		.map((f) => `${f.owner} ${f.field || "(value)"} = ${String(f.value)}`)
		.join("; ");
}
