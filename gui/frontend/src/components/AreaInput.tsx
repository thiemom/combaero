import type React from "react";
import useStore from "../store/useStore";
import NumericInput from "./NumericInput";

interface AreaInputProps {
	label: string;
	// null/undefined = no value (e.g. an inherited area the solver cannot
	// resolve yet): shows the placeholder rather than a coerced 0.
	value?: number | null;
	onChange: (val: number) => void;
	id?: string;
	// A number is an area in m^2 (e.g. an inherited value shown greyed until
	// the user overrides it) and is scaled to the display unit.
	placeholder?: string | number;
	onClear?: () => void;
	inputClassName?: string;
	readOnly?: boolean;
}

const AreaInput: React.FC<AreaInputProps> = ({
	label,
	value,
	onChange,
	id,
	onClear,
	placeholder = "0.00",
	inputClassName,
	readOnly,
}) => {
	const { unitPreferences, setAreaUnit } = useStore();
	const unit = unitPreferences.area;

	// Scale factor to convert internal (m2) to display (m2 or mm2)
	// 1 m2 = 1,000,000 mm2
	const displayScale = unit === "mm²" ? 1000000 : 1;

	const scaledPlaceholder =
		typeof placeholder === "number"
			? Number((placeholder * displayScale).toPrecision(4)).toString()
			: placeholder;

	return (
		<div className="flex flex-col gap-1">
			<label htmlFor={id} className="text-xs font-medium text-stone-500">
				{label}
			</label>
			<div className="flex gap-1">
				<NumericInput
					id={id}
					value={
						value !== undefined && value !== null
							? value * displayScale
							: undefined
					}
					onChange={(v) => onChange(v / displayScale)}
					onClear={onClear}
					readOnly={readOnly}
					className={`min-w-0 w-full p-1.5 border rounded text-xs bg-stone-50 focus:bg-white focus:ring-1 focus:ring-orange-500 outline-none${inputClassName ? ` ${inputClassName}` : ""}`}
					placeholder={scaledPlaceholder}
				/>
				<select
					value={unit}
					onChange={(e) => setAreaUnit(e.target.value as "m²" | "mm²")}
					className="w-14 h-[30px] border rounded text-[9px] bg-white outline-none shrink-0"
				>
					<option value="m²">m²</option>
					<option value="mm²">mm²</option>
				</select>
			</div>
		</div>
	);
};

export default AreaInput;
