import { afterEach, describe, expect, it, vi } from "vitest";
import {
	AUTOSAVE_KEY,
	clearAutosave,
	readAutosave,
	restoredNote,
	serialize,
	startAutosave,
	writeAutosave,
} from "./autosave";

class MemoryStorage {
	data = new Map<string, string>();
	getItem(k: string) {
		return this.data.has(k) ? (this.data.get(k) as string) : null;
	}
	setItem(k: string, v: string) {
		this.data.set(k, v);
	}
	removeItem(k: string) {
		this.data.delete(k);
	}
}

const throwing = {
	getItem: () => {
		throw new Error("blocked");
	},
	setItem: () => {
		throw new Error("quota");
	},
	removeItem: () => {
		throw new Error("blocked");
	},
};

const snap = () => ({
	nodes: [{ id: "a", data: { label: "A", result: { P: 1 } } }],
	edges: [
		{ id: "e", source: "a", target: "b", data: { result: { m_dot: 1 } } },
	],
	solverSettings: { global_regime: "compressible" },
});

describe("autosave", () => {
	afterEach(() => vi.useRealTimers());

	it("saves what Save saves, without solve results", () => {
		const d = serialize(snap(), new Date("2026-10-17T14:32:00Z"));
		expect(d.version).toBe("1.1");
		expect(d.nodes[0].data).toEqual({ label: "A" });
		expect(d.edges[0].data).toEqual({});
		expect(d.solverSettings).toEqual({ global_regime: "compressible" });
		expect(d.savedAt).toBe("2026-10-17T14:32:00.000Z");
	});

	it("round-trips through storage and clears", () => {
		const st = new MemoryStorage();
		expect(writeAutosave(st, snap())).toBe(true);
		const back = readAutosave(st);
		expect(back?.nodes[0].id).toBe("a");
		clearAutosave(st);
		expect(readAutosave(st)).toBeNull();
	});

	it("restores nothing from an empty canvas or a corrupt copy", () => {
		const st = new MemoryStorage();
		writeAutosave(st, { nodes: [], edges: [], solverSettings: {} });
		expect(readAutosave(st)).toBeNull();
		st.setItem(AUTOSAVE_KEY, "{not json");
		expect(readAutosave(st)).toBeNull();
		st.setItem(AUTOSAVE_KEY, JSON.stringify({ nodes: "x" }));
		expect(readAutosave(st)).toBeNull();
	});

	it("never throws when storage is blocked or full", () => {
		expect(writeAutosave(throwing, snap())).toBe(false);
		expect(readAutosave(throwing)).toBeNull();
		expect(() => clearAutosave(throwing)).not.toThrow();
	});

	it("writes once per burst of edits, and only for the saved parts", () => {
		vi.useFakeTimers();
		let state = { ...snap(), other: 0 };
		const listeners: ((s: any, p: any) => void)[] = [];
		const store = {
			getState: () => state,
			subscribe: (l: (s: any, p: any) => void) => {
				listeners.push(l);
				return () => listeners.splice(listeners.indexOf(l), 1);
			},
		};
		const set = (patch: any) => {
			const prev = state;
			state = { ...state, ...patch };
			for (const l of listeners) l(state, prev);
		};
		const st = new MemoryStorage();
		const writes = vi.spyOn(st, "setItem");
		const stop = startAutosave(store, st, 1000);
		set({ other: 1 }); // not a saved part
		vi.advanceTimersByTime(2000);
		expect(writes).not.toHaveBeenCalled();
		for (let i = 0; i < 5; i++) {
			set({ nodes: [...state.nodes] });
			vi.advanceTimersByTime(200);
		}
		expect(writes).not.toHaveBeenCalled();
		vi.advanceTimersByTime(1000);
		expect(writes).toHaveBeenCalledTimes(1);
		stop();
		set({ nodes: [] });
		vi.advanceTimersByTime(2000);
		expect(writes).toHaveBeenCalledTimes(1);
	});

	it("tells the user the canvas was restored, with the time", () => {
		expect(restoredNote("2026-10-17T14:32:00Z")).toMatch(
			/^Restored your last session \(autosaved .+\)\.$/,
		);
		expect(restoredNote("garbage")).toBe("Restored your last session.");
	});
});
