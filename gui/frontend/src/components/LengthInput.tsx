import type React from "react";
import useStore from "../store/useStore";
import NumericInput from "./NumericInput";

interface LengthInputProps {
	label: string;
	value: number;
	onChange: (val: number) => void;
	id?: string;
	placeholder?: string;
}

const LengthInput: React.FC<LengthInputProps> = ({
	label,
	value,
	onChange,
	id,
	placeholder = "0.00",
}) => {
	const { unitPreferences, setLengthUnit } = useStore();
	const unit = unitPreferences.length;

	// Scale factor to convert internal (m) to display (m or mm)
	const displayScale = unit === "mm" ? 1000 : 1;

	return (
		<div className="flex flex-col gap-1">
			<label htmlFor={id} className="text-xs font-medium text-stone-500">
				{label}
			</label>
			<div className="flex gap-1">
				<NumericInput
					id={id}
					value={value * displayScale}
					onChange={(v) => onChange(v / displayScale)}
					className="min-w-0 w-full p-1.5 border rounded text-xs bg-stone-50 focus:bg-white focus:ring-1 focus:ring-orange-500 outline-none"
					placeholder={placeholder}
				/>
				<select
					value={unit}
					onChange={(e) => setLengthUnit(e.target.value as "m" | "mm")}
					className="w-14 h-[30px] border rounded text-[9px] bg-white outline-none shrink-0"
				>
					<option value="m">m</option>
					<option value="mm">mm</option>
				</select>
			</div>
		</div>
	);
};

export default LengthInput;
