import type React from "react";
import { geometryView, storedFrom } from "../utils/geometry";
import { inheritedNote, type ResolvedField } from "../utils/inheritance";
import AreaInput from "./AreaInput";
import LengthInput from "./LengthInput";

// One flow face's cross-section -- D | Dh | A -- rendered identically on every
// element (#559). By default the face inherits: the solver's resolved values
// are shown greyed and locked. Unticking "Inherit" opts out: D and A become
// inputs (two views of one stored value; entering either computes the other)
// and Dh can be overridden for a non-circular duct. Ticking it again clears
// the override.

interface GeometryFieldsProps {
	idPrefix: string;
	title?: string;
	stores: "D" | "A";
	userValue: number | null | undefined;
	inherited: ResolvedField | undefined;
	onChange: (stored: number | null) => void;
	// Present when the element stores its own Dh (else Dh = D, computed).
	dh?: {
		userValue: number | null | undefined;
		inherited?: ResolvedField;
		onChange: (v: number | null) => void;
	};
	inheritable?: boolean;
	// Amber hint when nothing can be inherited and no value is set.
	missingHint?: string;
	// Replaces the generic note when the inherited value is a nominal default.
	defaultHint?: string;
}

const NOTE = "text-[9px] italic text-gray-400";

const GeometryFields: React.FC<GeometryFieldsProps> = ({
	idPrefix,
	title,
	stores,
	userValue,
	inherited,
	onChange,
	dh,
	inheritable,
	missingHint,
	defaultHint,
}) => {
	const v = geometryView({
		stores,
		userValue,
		inherited,
		userDh: dh?.userValue,
		inheritedDh: dh?.inherited,
		inheritable,
	});

	const setInherit = (inherit: boolean) => {
		if (inherit) {
			onChange(null);
			dh?.onChange(null);
		} else if (inherited?.value != null) {
			// Opting out pins the current inherited value as the user's own, so
			// the face stops following its neighbours from here on.
			onChange(inherited.value);
		}
	};

	const valueClass = (user: boolean) =>
		v.locked
			? "text-gray-400 cursor-not-allowed"
			: user
				? "font-semibold"
				: "text-gray-400";
	const note = v.userSet ? "" : inheritedNote(v.source);
	const circular =
		v.D !== null && v.Dh !== null && Math.abs(v.Dh - v.D) <= 1e-12 * v.D;

	return (
		<div className="flex flex-col gap-2 rounded border border-stone-200 p-2">
			<div className="flex items-center justify-between">
				<span className="text-xs font-bold uppercase text-gray-500">
					{title ?? "Cross-section"}
				</span>
				{v.canInherit && (
					<label className="flex items-center gap-1 text-[10px] text-gray-500">
						<input
							type="checkbox"
							checked={v.locked}
							onChange={(e) => setInherit(e.target.checked)}
						/>
						Inherit
					</label>
				)}
			</div>
			<LengthInput
				id={`${idPrefix}_D`}
				label="Diameter D"
				value={v.valueD}
				placeholder={v.D ?? (stores === "D" ? "set D" : "set D or A")}
				readOnly={v.locked}
				onChange={(val) => onChange(storedFrom(stores, "D", val))}
				onClear={() => onChange(null)}
				inputClassName={valueClass(v.userSet && stores === "D")}
			/>
			{dh ? (
				<LengthInput
					id={`${idPrefix}_Dh`}
					label="Hydraulic diameter Dh"
					value={v.userSetDh ? v.Dh : null}
					placeholder={v.Dh ?? "= D"}
					readOnly={v.locked}
					onChange={(val) => dh.onChange(val)}
					onClear={() => dh.onChange(null)}
					inputClassName={valueClass(v.userSetDh)}
				/>
			) : (
				<LengthInput
					id={`${idPrefix}_Dh`}
					label="Hydraulic diameter Dh"
					value={null}
					placeholder={v.D ?? "= D"}
					readOnly
					onChange={() => {}}
					inputClassName="text-gray-400 cursor-not-allowed"
				/>
			)}
			<AreaInput
				id={`${idPrefix}_A`}
				label="Area A"
				value={v.valueA}
				placeholder={v.A ?? (stores === "A" ? "set A" : "set D or A")}
				readOnly={v.locked}
				onChange={(val) => onChange(storedFrom(stores, "A", val))}
				onClear={() => onChange(null)}
				inputClassName={valueClass(v.userSet && stores === "A")}
			/>
			{note && <p className={NOTE}>{note}.</p>}
			{!v.userSet && v.source === "default" && (
				<p className="text-[9px] text-amber-700">
					{defaultHint ??
						(inheritable === false
							? "Not set — a nominal default is used. Enter D or A to set it."
							: "Nothing to inherit from — a nominal default is used. Untick Inherit to set it.")}
				</p>
			)}
			{v.userSet && (
				<p className={NOTE}>
					{stores === "D" ? "D" : "A"} is set; the other is computed (A =
					πD²/4).
				</p>
			)}
			{dh && !v.userSetDh && !circular && v.Dh !== null && (
				<p className={NOTE}>Dh: {inheritedNote(v.dhSource) || v.dhSource}.</p>
			)}
			{(dh ? !v.userSetDh && circular : true) && v.D !== null && (
				<p className={NOTE}>
					Dh = D (circular).
					{dh && !v.locked ? " Set Dh for a non-circular duct." : ""}
				</p>
			)}
			{v.D === null && missingHint && (
				<p className="text-[9px] text-amber-700">{missingHint}</p>
			)}
		</div>
	);
};

export default GeometryFields;
