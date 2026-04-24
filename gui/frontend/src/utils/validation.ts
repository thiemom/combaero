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
			if ((node.data.D || 0) <= 0)
				errors.push(`${node.id}: Channel diameter must be > 0`);
			if ((node.data.L || 0) <= 0)
				errors.push(`${node.id}: Channel length must be > 0`);
		}
		if (node.type === "orifice") {
			if ((node.data.area || 0) <= 0 && (node.data.diameter || 0) <= 0)
				errors.push(`${node.id}: Orifice area or diameter must be > 0`);
		}
		if (node.type === "momentum_chamber") {
			if ((node.data.area ?? 0) <= 0)
				errors.push(`${node.id}: Momentum Chamber area must be > 0`);
			if (node.data.Dh !== undefined && node.data.Dh < 0)
				errors.push(
					`${node.id}: Momentum Chamber Dh must be > 0 (or 0 for auto)`,
				);
		}
		if (node.type === "combustor") {
			if ((node.data.area ?? 0.1) <= 0)
				errors.push(`${node.id}: Combustor area must be > 0`);
			if (node.data.Dh !== undefined && node.data.Dh < 0)
				errors.push(`${node.id}: Combustor Dh must be > 0 (or 0 for auto)`);
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
		if (node.type === "area_change") {
			if ((node.data.F0 ?? 0.01) <= 0)
				errors.push(`${node.id}: AreaChange Upstream Area (F0) must be > 0`);
			if ((node.data.F1 ?? 0.01) <= 0)
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
