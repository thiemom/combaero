/**
 * Default node-card rows per palette kind (#775).
 *
 * A node card has room for one to three result rows. They used to be the
 * first keys of the global Display Fields list that resolved, in discovery
 * order, so a channel card showed m_dot, then Dh, then M_exit_duct. Now each
 * kind has a short curated list, most useful first: a card shows the first
 * `maxRows` of them that its result carries. The user can still change them
 * per kind (pins in the Live Telemetry, `cardRows` in the store); the Live
 * Telemetry itself keeps showing every selected Display Field.
 *
 * Keys are spelled as resultKeys.generated.json spells them (`state.<k>` for
 * a node's state); cardRows.test.ts fails when one is not emitted for its kind,
 * so a backend rename breaks the test rather than blanking the card.
 *
 * Only the kinds whose card draws result rows are listed: `CardKind` is the
 * `kind` prop of NodeDiagRows, so a new card without an entry does not compile.
 */
export const CARD_ROW_DEFAULTS = {
	// Nodes (three rows fit). Tt as the Display Fields default always had;
	// in a plenum P and Pt agree to the residual velocity.
	plenum: ["state.P", "state.Tt"],
	combustor: ["state.P", "state.Tt", "state.phi"],
	momentum_chamber: ["state.P", "state.Tt", "state.mach"],
	// Flow elements (one row fits today; the rest show once pinned first).
	channel: ["m_dot", "dP", "mach_out"],
	orifice: ["m_dot", "dP", "Cd"],
	area_change: ["m_dot", "dPt", "zeta"],
	discrete_loss: ["m_dot", "dPt", "mach_in"],
	vortex: ["m_dot", "dP", "M_theta_out"],
	impingement_array: ["m_dot", "dP", "htc_mean"],
	effusion_plate: ["m_dot", "T_wall_hot", "dP"],
	film_cooling: ["m_dot", "dP", "Cd"],
	slot_film: ["m_dot", "dP", "Cd"],
	effusion_liner: ["m_dot", "T_wall_hot_max", "bleed_fraction"],
	ejector: ["entrainment_ratio", "m_dot_primary", "m_dot_secondary"],
	mpce_tee: ["mass_flow_ratio", "K_branch", "K_straight"],
} as const satisfies Record<string, readonly string[]>;

export type CardKind = keyof typeof CARD_ROW_DEFAULTS;

/** True for a kind whose card draws result rows. */
export function isCardKind(kind: string | undefined): kind is CardKind {
	return kind !== undefined && Object.hasOwn(CARD_ROW_DEFAULTS, kind);
}

/** The keys a card of this kind shows, in order: the user's choice for the
 *  kind if they made one, else the curated default. */
export function cardRowKeys(
	kind: CardKind,
	overrides: Partial<Record<string, string[]>>,
): readonly string[] {
	return overrides[kind] ?? CARD_ROW_DEFAULTS[kind];
}

const rawKey = (k: string) => (k.startsWith("state.") ? k.slice(6) : k);

/** True when `key` is one of the card keys; `P` and `state.P` are the same
 *  field (the Live Telemetry lists a node's state keys bare). */
export function isCardRow(current: readonly string[], key: string): boolean {
	return current.some((k) => rawKey(k) === rawKey(key));
}

/** The card keys after pinning or unpinning `key`. A newly pinned key goes
 *  first, so it shows even on a one-row card. */
export function toggleCardRow(
	current: readonly string[],
	key: string,
): string[] {
	return isCardRow(current, key)
		? current.filter((k) => rawKey(k) !== rawKey(key))
		: [key, ...current];
}
