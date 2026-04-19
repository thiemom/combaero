import type React from "react";
import useStore from "../store/useStore";
import NumericInput from "./NumericInput";

interface UnitInputProps {
	label: string;
	value: number;
	onChange: (val: number) => void;
	id?: string;
}

const UnitInput: React.FC<UnitInputProps> = ({
	label,
	value,
	onChange,
	id,
}) => {
	const { unitPreferences, setPressureUnit } = useStore();
	const unit = unitPreferences.pressure;

	// Scale factor based on unit
	const scale = unit === "kPa" ? 1e3 : unit === "MPa" ? 1e6 : 1;

	return (
		<div className="flex flex-col gap-1">
			<label htmlFor={id} className="text-xs font-medium text-stone-500">
				{label}
			</label>
			<div className="flex gap-1">
				<NumericInput
					id={id}
					value={value / scale}
					onChange={(v) => onChange(v * scale)}
					className="flex-grow p-1.5 border rounded text-sm bg-stone-50 focus:bg-white focus:ring-1 focus:ring-orange-500 outline-none"
					placeholder="0.00"
				/>
				<select
					value={unit}
					onChange={(e) =>
						setPressureUnit(e.target.value as "Pa" | "kPa" | "MPa")
					}
					className="w-16 p-1.5 border rounded text-xs bg-white outline-none"
				>
					<option value="Pa">Pa</option>
					<option value="kPa">kPa</option>
					<option value="MPa">MPa</option>
				</select>
			</div>
		</div>
	);
};

export default UnitInput;
