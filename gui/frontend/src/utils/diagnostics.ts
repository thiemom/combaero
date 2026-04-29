/**
 * Network-wide diagnostic field discovery.
 *
 * Scans a solve result for every finite numeric scalar across all node,
 * element, and thermal-wall results. The output is topology-independent —
 * no field names are hardcoded; everything comes from the live solve data.
 */

export interface DiscoveredField {
	key: string;
	/** True when at least one node result carries a finite value for this key. */
	inNodes: boolean;
	/** True when at least one element result carries a finite value for this key. */
	inElements: boolean;
	/** True when at least one thermal-wall result carries a finite value for this key. */
	inEdges: boolean;
}

function isFiniteNum(v: unknown): v is number {
	return typeof v === "number" && Number.isFinite(v);
}

/**
 * Resolve a diagnostic value from a result object produced by the solver.
 *
 * Node results store state fields inside a `state` sub-object; element and
 * wall results keep everything at the top level. This helper tries both paths
 * so callers do not need to know which kind of result they hold.
 */
export function resolveField(result: unknown, key: string): number | undefined {
	if (!result || typeof result !== "object") return undefined;
	const r = result as Record<string, unknown>;
	const direct = r[key];
	if (isFiniteNum(direct)) return direct;
	const state = r.state as Record<string, unknown> | undefined;
	const nested = state?.[key];
	if (isFiniteNum(nested)) return nested;
	return undefined;
}

/**
 * Discover all finite scalar diagnostic fields present in a solve result.
 *
 * Returns fields sorted alphabetically. Fields that appear only as
 * arrays (Y, X) or booleans (success, is_correlation) are excluded.
 */
export function discoverFields(solveResults: unknown): DiscoveredField[] {
	if (!solveResults || typeof solveResults !== "object") return [];

	const r = solveResults as Record<string, unknown>;
	const nodeKeys = new Set<string>();
	const elemKeys = new Set<string>();
	const edgeKeys = new Set<string>();

	const SKIP = new Set(["Y", "X", "success", "is_correlation"]);

	for (const nr of Object.values(
		(r.node_results as Record<string, unknown>) ?? {},
	)) {
		const nrObj = nr as Record<string, unknown>;
		const state = nrObj.state as Record<string, unknown> | undefined;
		for (const [k, v] of Object.entries(state ?? {})) {
			if (!SKIP.has(k) && isFiniteNum(v)) nodeKeys.add(k);
		}
		// Also scan top-level node fields (e.g. phi, theta on combustor results
		// which are not part of the thermodynamic state sub-object).
		for (const [k, v] of Object.entries(nrObj)) {
			if (k !== "state" && !SKIP.has(k) && isFiniteNum(v)) nodeKeys.add(k);
		}
	}

	for (const er of Object.values(
		(r.element_results as Record<string, unknown>) ?? {},
	)) {
		for (const [k, v] of Object.entries(er as Record<string, unknown>)) {
			if (!SKIP.has(k) && isFiniteNum(v)) elemKeys.add(k);
		}
	}

	for (const wr of Object.values(
		(r.edge_results as Record<string, unknown>) ?? {},
	)) {
		for (const [k, v] of Object.entries(wr as Record<string, unknown>)) {
			if (!SKIP.has(k) && isFiniteNum(v)) edgeKeys.add(k);
		}
	}

	const allKeys = new Set([...nodeKeys, ...elemKeys, ...edgeKeys]);

	return [...allKeys].sort().map((key) => ({
		key,
		inNodes: nodeKeys.has(key),
		inElements: elemKeys.has(key),
		inEdges: edgeKeys.has(key),
	}));
}
