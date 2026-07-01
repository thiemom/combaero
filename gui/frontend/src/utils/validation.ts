/**
 * Shared pre-solve network validation.
 *
 * Returns a list of human-readable error strings, one per invalid node/edge
 * field. The solve button shows these via `alert()` and the inspector panel
 * surfaces them in a warning box, so both entry points stay consistent.
 */

import type { Node } from "reactflow";

export function validateNetwork(nodes: Node[]): string[] {
	const errors: string[] = [];

	if (nodes.length === 0) {
		errors.push("Network is empty. Add nodes and elements before solving.");
		return errors;
	}

	for (const node of nodes) {
		if (node.type === "pressure_boundary" || node.type === "mass_boundary") {
			if ((node.data.Pt || 0) <= 0 && node.type === "pressure_boundary") {
				errors.push(`${node.id}: Total Pressure must be > 0`);
			}
			if ((node.data.Tt || 0) <= 0) {
				errors.push(`${node.id}: Temperature must be > 0`);
			}
			if (node.data.composition?.source === "custom") {
				const values = Object.values(
					node.data.composition.custom_fractions || {},
				) as number[];
				const sum = values.reduce((a: number, b: number) => a + (b || 0), 0);
				if (Math.abs(sum - 1.0) > 1e-3) {
					errors.push(
						`${node.id}: Custom composition sum is ${sum.toFixed(3)} (expected 1.0)`,
					);
				}
			}
		}
		if (node.type === "channel") {
			// null/undefined D means "inherit from graph" — only flag explicit non-positive values
			if (node.data.D != null && node.data.D <= 0)
				errors.push(`${node.id}: Channel diameter must be > 0`);
			if ((node.data.L || 0) <= 0)
				errors.push(`${node.id}: Channel length must be > 0`);
		}
		if (node.type === "orifice") {
			// null diameter means "use default 0.08 m" — only flag explicit non-positive values
			if (node.data.diameter != null && node.data.diameter <= 0)
				errors.push(`${node.id}: Orifice diameter must be > 0`);
		}
		if (node.type === "momentum_chamber") {
			// null area means "derive from Dh" — only flag explicit non-positive values
			if (node.data.area != null && node.data.area <= 0)
				errors.push(`${node.id}: Momentum Chamber area must be > 0`);
			if (node.data.Dh != null && node.data.Dh <= 0)
				errors.push(
					`${node.id}: Momentum Chamber Dh must be > 0 (or left blank to inherit)`,
				);
		}
		if (node.type === "combustor") {
			// null area means "derive from Dh" — only flag explicit non-positive values
			if (node.data.area != null && node.data.area <= 0)
				errors.push(`${node.id}: Combustor area must be > 0`);
			if (node.data.Dh != null && node.data.Dh <= 0)
				errors.push(
					`${node.id}: Combustor Dh must be > 0 (or left blank to inherit)`,
				);
		}
		if (node.type === "discrete_loss") {
			// `area` is optional: when undefined/null it is inferred from the
			// upstream node in the backend graph builder. We only reject
			// *explicit* non-positive values, which would otherwise silently
			// fall back to a 0.1 m² default and produce wrong HTC/Nu results.
			const area = node.data.area;
			if (area !== undefined && area !== null && area <= 0) {
				errors.push(
					`${node.id}: Discrete loss area must be > 0 (or left blank to auto-infer)`,
				);
			}
			// Head-type correlations use zeta * 0.5 * rho * v^2 / P and require
			// a non-zero zeta coefficient; zero would trivially zero the loss.
			const corr = node.data.correlation_type;
			if (
				(corr === "constant_head" || corr === "linear_theta_head") &&
				(node.data.zeta ?? 0) <= 0 &&
				(node.data.zeta0 ?? 0) <= 0
			) {
				errors.push(`${node.id}: Discrete loss zeta must be > 0`);
			}
		}
		if (node.type === "mpce_tee") {
			const thetaDeg = node.data.theta_deg ?? 90;
			if (Math.abs(thetaDeg) > 90) {
				errors.push(
					`${node.id}: Tee branch angle |theta| must be <= 90 deg, got ${thetaDeg} deg`,
				);
			}
			// null F_C means "inherit from connected channels" — only flag explicit non-positive values
			if (node.data.F_C != null && node.data.F_C <= 0)
				errors.push(`${node.id}: Tee common arm area F_C must be > 0`);
			if ((node.data.psi ?? 1) <= 0)
				errors.push(`${node.id}: Tee area ratio psi must be > 0`);
		}
		if (node.type === "area_change") {
			// null F0/F1 means "inherit from connected channels" — only flag explicit non-positive values
			if (node.data.F0 != null && node.data.F0 <= 0)
				errors.push(`${node.id}: AreaChange Upstream Area (F0) must be > 0`);
			if (node.data.F1 != null && node.data.F1 <= 0)
				errors.push(`${node.id}: AreaChange Downstream Area (F1) must be > 0`);
			if (
				node.data.D_h !== undefined &&
				node.data.D_h !== null &&
				node.data.D_h < 0
			)
				errors.push(
					`${node.id}: AreaChange Dh must be >= 0 (or left blank for auto)`,
				);
		}
	}

	return errors;
}
