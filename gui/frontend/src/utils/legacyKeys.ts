// Renamed input keys in saved networks (#773). A file written before a
// rename still loads: the old key moves to the new one unless the new one is
// already set. The backend accepts the old names too (pydantic aliases in
// gui/backend/schemas.py), so a headless payload needs no migration.

/** omega_rpm became shaft_speed_rpm: omega is 2 pi f, not a speed in rpm. */
const RENAMED: Record<string, string> = { omega_rpm: "shaft_speed_rpm" };

const renameKeys = (obj: any): any => {
	if (!obj || typeof obj !== "object") return obj;
	let out = obj;
	for (const [oldKey, newKey] of Object.entries(RENAMED)) {
		if (!(oldKey in out)) continue;
		const { [oldKey]: value, ...rest } = out;
		out =
			newKey in rest && rest[newKey] != null
				? rest
				: { ...rest, [newKey]: value };
	}
	return out;
};

/** Node data (the vortex's shaft speed) under its current key names. */
export const migrateNodeData = (data: any): any => renameKeys(data);

/** Solver settings (the global shaft speed) under their current key names. */
export const migrateSolverSettings = (settings: any): any =>
	renameKeys(settings);
