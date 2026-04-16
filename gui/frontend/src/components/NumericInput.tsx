import type React from "react";
import { useEffect, useState } from "react";

interface NumericInputProps {
	value?: number | null;
	onChange: (val: number) => void;
	id?: string;
	className?: string;
	placeholder?: string;
	step?: string;
}

const NumericInput: React.FC<NumericInputProps> = ({
	value,
	onChange,
	id,
	className,
	placeholder,
}) => {
	const [localValue, setLocalValue] = useState(
		value !== undefined && value !== null ? value.toString() : "",
	);
	const [isFocused, setIsFocused] = useState(false);

	// Update local state if the external value changes, but ONLY when not focused
	useEffect(() => {
		if (isFocused) return;
		if (value === undefined || value === null) {
			setLocalValue("");
			return;
		}
		const safeValue = value;
		const currentLocal = Number.parseFloat(localValue);
		if (
			Number.isNaN(currentLocal) ||
			Math.abs(currentLocal - safeValue) > 1e-10
		) {
			setLocalValue(safeValue.toString());
		}
	}, [value, isFocused, localValue]);

	const handleChange = (e: React.ChangeEvent<HTMLInputElement>) => {
		const raw = e.target.value;
		// Allow typing decimal points and negative signs
		if (raw === "" || raw === "-" || raw === "." || raw === "-.") {
			setLocalValue(raw);
			return;
		}

		setLocalValue(raw);

		const parsed = Number.parseFloat(raw);
		if (!Number.isNaN(parsed) && /^-?\d*\.?\d*$/.test(raw)) {
			onChange(parsed);
		}
	};

	const handleBlur = () => {
		setIsFocused(false);
		const parsed = Number.parseFloat(localValue);
		if (Number.isNaN(parsed)) {
			setLocalValue(value.toString());
		} else {
			onChange(parsed);
			setLocalValue(parsed.toString()); // Force clean format on blur
		}
	};

	return (
		<input
			id={id}
			type="text"
			value={localValue}
			onChange={handleChange}
			onFocus={() => setIsFocused(true)}
			onBlur={handleBlur}
			className={className}
			placeholder={placeholder}
		/>
	);
};

export default NumericInput;
