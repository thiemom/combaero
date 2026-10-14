/**
 * The second-law verdicts and the network balance (#545), for display.
 *
 * Every converged solve reports, without changing its success:
 * - per element, whether its control volume destroys entropy (a passive
 *   element creating total pressure), attributed to the canvas node that
 *   owns it -- a tee for its hidden port step, a liner for its segments;
 * - the ledger: mass and energy entering and leaving the network.
 */

import { nodeTypeLabel } from "./nodeUtils";

export interface SecondLawFlag {
	element_id: string;
	owner_id: string;
	sigma_norm: number;
	sigma_norm_high: number;
}

export interface LedgerStream {
	node: string;
	via: string;
	via_owner: string;
	/** wired node-to-node: via_owner is the node at the edge's far end */
	direct?: boolean;
	m_dot_in: number;
	H_in: number;
}

export interface Ledger {
	streams: LedgerStream[];
	m_dot_in: number;
	m_dot_out: number;
	mass_closure_rel: number;
	H_in: number;
	H_out: number;
	Q_given: number;
	Q_wall_out: number;
	Q_withheld: number;
	energy_closure: number;
	energy_closure_rel: number;
}

/**
 * What the canvas calls a node: its label, else its kind with the end of its
 * id, so two unlabelled boundaries stay apart without a 20-digit id.
 */
export function displayName(
	id: string,
	nodes: { id: string; type?: string; data?: { label?: string } }[],
): string {
	const n = nodes.find((x) => x.id === id);
	if (!n) return id;
	if (n.data?.label) return n.data.label;
	const kind = nodeTypeLabel(n.type);
	return kind ? `${kind} ${id.slice(-4)}` : id;
}

/** Relative closure above which the balance is reported as open. */
export const LEDGER_TOL = 1e-6;

export function failuresFor(
	nodeId: string,
	failures: SecondLawFlag[] | undefined,
): SecondLawFlag[] {
	return (failures ?? []).filter((f) => f.owner_id === nodeId);
}

export function balanceIsOpen(ledger: Ledger | undefined | null): boolean {
	if (!ledger || ledger.m_dot_in === undefined) return false;
	return (
		Math.abs(ledger.mass_closure_rel) > LEDGER_TOL ||
		Math.abs(ledger.energy_closure_rel) > LEDGER_TOL
	);
}

/** "the area step at its straight port" for a hidden element, else its id. */
export function describeElement(flag: SecondLawFlag): string {
	const step = flag.element_id.match(
		/^__tee_step__.+_(common|straight|branch)$/,
	);
	if (step) return `the area step at its ${step[1]} port`;
	if (flag.element_id === flag.owner_id) return "this element";
	return flag.element_id;
}
