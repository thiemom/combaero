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
	max?: number;
	disabled?: boolean;
	// Shown but not editable (an inherited value locked until the user opts out).
	readOnly?: boolean;
	/**
	 * Display-only format (e.g. a value copied from a solve, shown at 6
	 * significant figures). The stored value stays exact: leaving the field
	 * without typing in it never commits the rounded display, and what the
	 * user types replaces it as usual. Without it the field shows
	 * `value.toString()`.
	 */
	format?: (v: number) => string;
}

const NumericInput: React.FC<NumericInputProps> = ({
	value,
	onChange,
	onClear,
	id,
	className,
	placeholder,
	min,
	max,
	disabled,
	readOnly,
	format,
}) => {
	const [localValue, setLocalValue] = useState(
		value !== undefined && value !== null ? (format ?? String)(value) : "",
	);
	const [isFocused, setIsFocused] = useState(false);
	// Typed in since focus: only then may a formatted display be committed.
	const [dirty, setDirty] = useState(false);

	// Update local state if the external value changes, but ONLY when not focused
	useEffect(() => {
		if (isFocused) return;
		if (value === undefined || value === null) {
			setLocalValue("");
			return;
		}
		if (format) {
			const shown = format(value);
			if (localValue !== shown) setLocalValue(shown);
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
	}, [value, isFocused, localValue, format]);

	const handleChange = (e: React.ChangeEvent<HTMLInputElement>) => {
		const raw = e.target.value;
		setDirty(true);
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
			if (max !== undefined && parsed > max) return;
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
		if (format && !dirty) {
			// Clicked in and out: the display is rounded, the value is not.
			setLocalValue(value != null ? format(value) : "");
			return;
		}
		const parsed = Number.parseFloat(localValue);
		// parseFloat reads "Infinity" and "1e999" as Infinity: not a value a
		// field can hold or a file can save (#1038). Snap back, as for text.
		if (!Number.isFinite(parsed)) {
			setLocalValue(value != null ? (format ?? String)(value) : "");
		} else {
			const clamped =
				min !== undefined && parsed < min
					? min
					: max !== undefined && parsed > max
						? max
						: parsed;
			// Only propagate if the value actually changed; clicking into a field
			// that shows an inherited/default value and clicking away without typing
			// should not mark it as user-set.
			if (value == null || Math.abs(clamped - value) > 1e-10) {
				onChange(clamped);
			}
			setLocalValue((format ?? String)(clamped));
		}
	};

	const handleFocus = () => {
		setIsFocused(true);
		setDirty(false);
	};

	return (
		<input
			id={id}
			type="text"
			value={localValue}
			onChange={handleChange}
			onFocus={handleFocus}
			onBlur={handleBlur}
			className={className}
			placeholder={placeholder}
			disabled={disabled}
			readOnly={readOnly}
		/>
	);
};

export default NumericInput;
