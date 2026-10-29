/**
 * Live Telemetry controls for the node-card rows of a kind (#775): a pin per
 * field, and a reset to the curated default (utils/cardRows.ts). The choice
 * applies to every card of the kind, like the Display Fields apply to every
 * Live Telemetry.
 */

import { Pin } from "lucide-react";
import useStore from "../store/useStore";
import {
	type CardKind,
	cardRowKeys,
	isCardRow,
	toggleCardRow,
} from "../utils/cardRows";
import { nodeTypeLabel } from "../utils/nodeUtils";

export function CardRowPin({
	kind,
	fieldKey,
}: {
	kind: CardKind;
	fieldKey: string;
}) {
	const cardRows = useStore((s) => s.cardRows);
	const setCardRows = useStore((s) => s.setCardRows);
	const current = cardRowKeys(kind, cardRows);
	const pinned = isCardRow(current, fieldKey);
	const what = `${nodeTypeLabel(kind)} cards`;
	return (
		<button
			type="button"
			onClick={() => setCardRows(kind, toggleCardRow(current, fieldKey))}
			title={
				pinned
					? `Shown on ${what}: click to remove`
					: `Show on ${what} (first, the card shows the top rows that fit)`
			}
			aria-label={pinned ? `Unpin from ${what}` : `Pin to ${what}`}
			aria-pressed={pinned}
			className={`shrink-0 ${pinned ? "text-blue-500" : "text-stone-300 hover:text-stone-500"}`}
		>
			<Pin size={9} fill={pinned ? "currentColor" : "none"} />
		</button>
	);
}

export function CardRowsReset({ kind }: { kind: CardKind }) {
	const overridden = useStore((s) => s.cardRows[kind] !== undefined);
	const setCardRows = useStore((s) => s.setCardRows);
	if (!overridden) return null;
	return (
		<button
			type="button"
			onClick={() => setCardRows(kind, null)}
			title={`Restore the default rows of ${nodeTypeLabel(kind)} cards`}
			className="text-[9px] font-bold uppercase tracking-tight px-1.5 py-0.5 rounded border border-stone-200 hover:bg-white text-stone-400"
		>
			Reset card rows
		</button>
	);
}
