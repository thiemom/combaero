import type React from "react";
import { useEffect, useState } from "react";

interface NumericInputProps {
	value?: number | null;
	onChange: (val: number) => void;
	/**
	 * Optional callback invoked when the user leaves the field empty on blur.
	 * When provided, the cleared state is propagated to the parent instead of
	 * silently snapping back to the last ``value`` prop. Use this for inputs
	 * where "empty" is a meaningful state (e.g. auto-inferred from topology).
	 */
	onClear?: () => void;
	id?: string;
	className?: string;
	placeholder?: string;
	step?: string;
	min?: number;
}

const NumericInput: React.FC<NumericInputProps> = ({
	value,
	onChange,
	onClear,
	id,
	className,
	placeholder,
	min,
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
		// Allow typing decimal points; only allow negative sign if no min constraint
		if (
			raw === "" ||
			raw === "." ||
			(min === undefined && (raw === "-" || raw === "-."))
		) {
			setLocalValue(raw);
			return;
		}

		setLocalValue(raw);

		const parsed = Number.parseFloat(raw);
		if (!Number.isNaN(parsed) && /^-?\d*\.?\d*$/.test(raw)) {
			if (min !== undefined && parsed < min) return;
			onChange(parsed);
		}
	};

	const handleBlur = () => {
		setIsFocused(false);
		// Intentional clear: empty input + parent accepts cleared state.
		if (localValue === "" && onClear) {
			onClear();
			return;
		}
		const parsed = Number.parseFloat(localValue);
		if (Number.isNaN(parsed)) {
			setLocalValue(value != null ? value.toString() : "");
		} else {
			const clamped = min !== undefined && parsed < min ? min : parsed;
			onChange(clamped);
			setLocalValue(clamped.toString());
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
