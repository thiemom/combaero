// Display mirror of ConstantKTeeElement.bassett_design_K
// (python/combaero/network/mpce_element.py), which owns the formula and is
// what the solve uses. Bassett 2001's analytical pair at a design split:
// dividing K5/K6, joining K11/K12 (Eq 33 angle-corrected). Both sides pin
// the same values in their tests, so a drift shows up as a failure.

export function bassettDesignK(
	flowDirection: "branch" | "merge",
	qBranch: number,
	psi: number,
	thetaDeg: number,
): [number, number] {
	const q = qBranch;
	const qs = 1 - q;
	const c = Math.cos(0.75 * ((thetaDeg * Math.PI) / 180));
	if (flowDirection === "branch") {
		const k5 = qs * qs - 1.5 * qs + 0.5;
		const k6 = 1 + q * q * psi * psi - 2 * q * psi * c;
		return [k5, k6];
	}
	const d = psi + 0.5 * c;
	const n = 1 - qs * qs - q * q * psi * c;
	return [
		((2 * psi) / d) * n + qs * qs - 1,
		((2 * psi) / d) * n + q * q * psi * psi - 1,
	];
}

// Mirror of ConstantKTeeElement.design_net_loss: the pair's mass-weighted
// total-pressure loss at its design share, per common-leg head. Below zero
// the tee creates total pressure and the backend refuses the pair (#779).
export function designNetLoss(
	qBranch: number,
	kStraight: number,
	kBranch: number,
): number {
	return (1 - qBranch) * kStraight + qBranch * kBranch;
}
