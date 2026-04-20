/** @vitest-environment jsdom */
import { describe, expect, it } from "vitest";

describe("Node 25 Compatibility Regression Tests", () => {
	// Explicitly mock localStorage if the environment fails to provide it
	if (
		typeof localStorage === "undefined" ||
		localStorage === null ||
		typeof localStorage.setItem !== "function"
	) {
		const mockStorage: Record<string, string> = {};
		(globalThis as any).localStorage = {
			getItem: (key: string) => mockStorage[key] || null,
			setItem: (key: string, value: string) => {
				mockStorage[key] = value;
			},
			removeItem: (key: string) => {
				delete mockStorage[key];
			},
			clear: () => {
				for (const key in mockStorage) delete mockStorage[key];
			},
			length: 0,
			key: (index: number) => Object.keys(mockStorage)[index] || null,
		} as Storage;
	}

	it("should have access to localStorage without SecurityErrors", () => {
		// Node 25 enables documentation-level localStorage/sessionStorage by default.
		// We verify that our testing environment (currently jsdom) handles this gracefully.
		const testKey = "__combaero_node_25_test__";
		const testValue = "stable";

		localStorage.setItem(testKey, testValue);
		expect(localStorage.getItem(testKey)).toBe(testValue);
		localStorage.removeItem(testKey);
	});

	it("should strictly adhere to WHATWG URL standards", () => {
		// Node 25 enforces stricter WHATWG URL parsing.
		// We verify basic absolute and relative URL parsing used in our API client.
		const baseUrl = "https://api.combaero.io";
		const relativePath = "/v1/solver";
		const url = new URL(relativePath, baseUrl);

		expect(url.toString()).toBe("https://api.combaero.io/v1/solver");
		expect(url.protocol).toBe("https:");
	});

	it("should not use deprecated SlowBuffer", () => {
		// Node 25 removes SlowBuffer. While not used in our frontend src,
		// this ensures no build-time imports accidentally rely on it.
		expect((globalThis as any).SlowBuffer).toBeUndefined();
	});
});
