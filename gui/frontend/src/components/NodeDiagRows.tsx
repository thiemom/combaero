/**
 * Shared micro-component for node cards.
 *
 * Renders up to `maxRows` diagnostic values from a solve result, using
 * whichever keys the user has selected in `displaySettings`. Only fields
 * that resolve to a finite numeric value in this specific result are shown â
 * stale or absent keys are silently skipped, so no hardcoded field lists are
 * needed here or in any consuming node component.
 */

import { useMemo } from "react";
import useStore from "../store/useStore";
import { resolveField } from "../utils/diagnostics";
import { quantityMeta } from "../utils/quantities";

interface Props {
	result: unknown;
	maxRows?: number;
}

export function NodeDiagRows({ result, maxRows = 3 }: Props) {
	const { displaySettings } = useStore();

	const rows = useMemo(() => {
		if (!result) return [];
		return displaySettings
			.map((key) => {
				const val = resolveField(result, key);
				if (val === undefined) return null;
				const meta = quantityMeta(key);
				return {
					key,
					label: meta.label,
					formatted: meta.format(val),
					unit: meta.unit === "—" ? "" : meta.unit,
				};
			})
			.filter(Boolean)
			.slice(0, maxRows) as {
			key: string;
			label: string;
			formatted: string;
			unit: string;
		}[];
	}, [result, displaySettings, maxRows]);

	if (rows.length === 0) return null;

	return (
		<>
			{rows.map((row) => (
				<div
					key={row.key}
					className="text-[9px] text-gray-500 font-mono whitespace-nowrap truncate max-w-full"
				>
					{row.label}: {row.formatted}{" "}
					<span className="text-gray-400">{row.unit}</span>
				</div>
			))}
		</>
	);
}
