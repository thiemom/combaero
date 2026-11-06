/**
 * Discrete loss with a fraction correlation (#1004): the "use the solved
 * flow as ṁ_ref" button. A fraction loss without ṁ_ref holds ΔPt = ξ·Pt at
 * every flow but zero; ṁ_ref anchors it (below it the loss fades with the
 * dynamic head). At ṁ = ṁ_ref the loss is exactly ξ, so setting ṁ_ref to
 * the solved flow leaves the solved point unchanged.
 */

/**
 * [kg/s] the flow over which a directional loss blends its forward and
 * reverse forms: FLOW_REVERSAL_M_SCALE in python/combaero/network/components.py.
 * Below it the element sits on that blend and carries no flow to copy.
 */
export const FLOW_REVERSAL_M_SCALE = 1e-4;

interface FractionLossData {
	m_ref?: number | null;
	result?: { success?: boolean; m_dot?: number | null } | null;
}

/**
 * The ṁ_ref the button sets, |ṁ| of the last converged solve, or null when
 * the button is off: no result, a failed solve, a flow on the zero-flow
 * blend, or a set ṁ_ref above the solved flow (the loss was faded there, so
 * replacing ṁ_ref would move the point).
 */
export function solvedFlowAsMref(data: FractionLossData): number | null {
	const r = data.result;
	if (r?.success !== true) return null;
	const m = Math.abs(Number(r.m_dot));
	if (!Number.isFinite(m) || m < FLOW_REVERSAL_M_SCALE) return null;
	const mRef = data.m_ref;
	if (mRef != null && mRef > 0 && m < mRef) return null;
	return m;
}
