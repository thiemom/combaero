/**
 * The incompressible validity report (#588), for display.
 *
 * Every converged solve checks each element on an incompressible closure
 * against one Mach limit (0.2 by default): above it the element is listed
 * here, attributed to the canvas node that owns it. Report only -- the
 * solution and its success are unchanged.
 */

export interface IncompressibleFlag {
	element_id: string;
	owner_id: string;
	mach: number;
	limit: number;
}

export function warningsFor(
	nodeId: string,
	flags: IncompressibleFlag[] | undefined,
): IncompressibleFlag[] {
	return (flags ?? []).filter((f) => f.owner_id === nodeId);
}

/** "2 elements above Ma 0.2 (incompressible)", or null. */
export function summary(
	flags: IncompressibleFlag[] | undefined,
): string | null {
	const n = flags?.length ?? 0;
	if (n === 0) return null;
	const limit = flags?.[0]?.limit ?? 0.2;
	return `${n} element${n > 1 ? "s" : ""} above Ma ${limit} (incompressible)`;
}
