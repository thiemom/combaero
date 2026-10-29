/**
 * Shared micro-component for node cards.
 *
 * Renders up to `maxRows` diagnostic values from a solve result: the card
 * rows of this kind (utils/cardRows.ts), which are the curated per-kind
 * default unless the user pinned others in the Live Telemetry (#775). Only
 * keys that resolve to a finite number in this result are shown, so a row a
 * variant does not carry is skipped and the next one takes its place.
 */

import { useMemo } from "react";
import useStore from "../store/useStore";
import { type CardKind, cardRowKeys } from "../utils/cardRows";
import { resolveKey } from "../utils/diagnostics";
import { quantityMeta } from "../utils/quantities";

interface Props {
	kind: CardKind;
	result: unknown;
	maxRows?: number;
}

export function NodeDiagRows({ kind, result, maxRows = 3 }: Props) {
	const cardRows = useStore((s) => s.cardRows);

	const rows = useMemo(() => {
		if (!result) return [];
		return cardRowKeys(kind, cardRows)
			.map((key) => {
				const val = resolveKey(result, key);
				if (val === undefined) return null;
				const meta = quantityMeta(key);
				return {
					key,
					label: meta.cardLabel ?? meta.label,
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
	}, [kind, result, cardRows, maxRows]);

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
