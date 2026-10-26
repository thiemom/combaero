/**
 * Which drawn arm a junction's `port_<i>_*` result key belongs to.
 *
 * The backend numbers a junction's ports inlets first, then outlets
 * (gui/backend/graph_builder.py):
 *   mpce_tee, flow_direction "branch": [common, straight, branch]
 *   mpce_tee, flow_direction "merge":  [straight, branch, common]
 *   ejector:                           [primary, secondary, outlet]
 * The letters are the ones the node cards draw at the handles.
 * python/tests/test_gui_port_arms.py solves both tee directions and the
 * ejector and checks each port's flow against the arm wired to it, reading
 * PORT_ARMS from this file.
 */

export const PORT_ARMS = {
	mpce_tee_branch: ["C", "S", "B"],
	mpce_tee_merge: ["S", "B", "C"],
	ejector: ["P", "S", "O"],
} as const;

/** The arm letters of a node's ports, in port order, or null when the node
 *  is not a multi-port junction. */
export function portArms(
	type: string | undefined,
	data: { flow_direction?: string } | undefined,
): readonly string[] | null {
	if (type === "mpce_tee") {
		return (data?.flow_direction ?? "branch") === "merge"
			? PORT_ARMS.mpce_tee_merge
			: PORT_ARMS.mpce_tee_branch;
	}
	if (type === "ejector") return PORT_ARMS.ejector;
	return null;
}

/** `label (arm)` for a `port_<i>_*` key of a junction; `label` otherwise. */
export function withPortArm(
	label: string,
	key: string,
	type: string | undefined,
	data: { flow_direction?: string } | undefined,
): string {
	const arms = portArms(type, data);
	const m = /^(?:state\.)?port_(\d+)_/.exec(key);
	if (!arms || !m) return label;
	const arm = arms[Number(m[1])];
	return arm ? `${label} (${arm})` : label;
}
