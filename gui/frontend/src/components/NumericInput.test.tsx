/** @vitest-environment jsdom */
import { act } from "react";
import { createRoot, type Root } from "react-dom/client";
import { afterEach, beforeEach, describe, expect, it, vi } from "vitest";
import NumericInput from "./NumericInput";

// #1038: parseFloat reads "Infinity" and "1e999" as Infinity. The field
// filtered them while typing but committed them on blur, and an infinite
// value was then saved as null.

(
	globalThis as { IS_REACT_ACT_ENVIRONMENT?: boolean }
).IS_REACT_ACT_ENVIRONMENT = true;

let host: HTMLDivElement;
let root: Root;

beforeEach(() => {
	host = document.createElement("div");
	document.body.appendChild(host);
	root = createRoot(host);
});

afterEach(() => {
	act(() => root.unmount());
	host.remove();
});

function typeAndBlur(text: string, onChange: (v: number) => void): void {
	act(() => {
		root.render(<NumericInput value={4} onChange={onChange} />);
	});
	const input = host.querySelector("input") as HTMLInputElement;
	const setValue = Object.getOwnPropertyDescriptor(
		HTMLInputElement.prototype,
		"value",
	)?.set as (v: string) => void;
	act(() => {
		input.focus();
	});
	act(() => {
		setValue.call(input, text);
		input.dispatchEvent(new Event("input", { bubbles: true }));
	});
	act(() => {
		input.blur();
	});
}

describe("NumericInput refuses a non-finite value on blur", () => {
	it.each(["Infinity", "1e999", "-Infinity"])("%s snaps back", (text) => {
		const onChange = vi.fn();
		typeAndBlur(text, onChange);
		for (const [v] of onChange.mock.calls)
			expect(Number.isFinite(v)).toBe(true);
	});

	it("commits a finite value", () => {
		const onChange = vi.fn();
		typeAndBlur("6", onChange);
		expect(onChange).toHaveBeenLastCalledWith(6);
	});
});
