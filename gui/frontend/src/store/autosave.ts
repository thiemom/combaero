// Auto-save the network to the browser's localStorage (#596).
//
// The copy is the same JSON the Save button writes ({nodes, edges,
// solverSettings, version}) without solve results, which are stale after any
// edit. A reload restores it through loadNetwork, so a restore behaves like
// opening a saved file. One slot shared by every tab: the last tab edited wins.

export const AUTOSAVE_KEY = "combaero.autosave.v1";
export const AUTOSAVE_DELAY_MS = 1000;

export interface AutosaveData {
	nodes: any[];
	edges: any[];
	solverSettings: Record<string, unknown>;
	version: string;
	savedAt: string;
}

interface StorageLike {
	getItem(key: string): string | null;
	setItem(key: string, value: string): void;
	removeItem(key: string): void;
}

interface Snapshot {
	nodes: any[];
	edges: any[];
	solverSettings: Record<string, unknown>;
}

const withoutResult = (items: any[]): any[] =>
	items.map((it) => {
		if (!it?.data || !("result" in it.data)) return it;
		const { result: _r, ...rest } = it.data;
		return { ...it, data: rest };
	});

export function serialize(s: Snapshot, now: Date = new Date()): AutosaveData {
	return {
		nodes: withoutResult(s.nodes),
		edges: withoutResult(s.edges),
		solverSettings: s.solverSettings,
		version: "1.1",
		savedAt: now.toISOString(),
	};
}

// A storage that throws (private mode, quota, blocked site data) must never
// break the editor: every access is guarded and a failure is a no-op.
export function writeAutosave(storage: StorageLike, s: Snapshot): boolean {
	try {
		storage.setItem(AUTOSAVE_KEY, JSON.stringify(serialize(s)));
		return true;
	} catch {
		return false;
	}
}

export function readAutosave(storage: StorageLike): AutosaveData | null {
	try {
		const raw = storage.getItem(AUTOSAVE_KEY);
		if (!raw) return null;
		const data = JSON.parse(raw);
		if (!Array.isArray(data?.nodes) || !Array.isArray(data?.edges)) return null;
		// An empty canvas is nothing to restore.
		if (data.nodes.length === 0 && data.edges.length === 0) return null;
		return data as AutosaveData;
	} catch {
		return null;
	}
}

export function clearAutosave(storage: StorageLike): void {
	try {
		storage.removeItem(AUTOSAVE_KEY);
	} catch {
		// nothing to do
	}
}

interface Subscribable<S> {
	getState(): S;
	subscribe(listener: (state: S, prev: S) => void): () => void;
}

// Write on every change of the saved parts, debounced: one write per burst of
// edits (a drag moves nodes on every frame).
export function startAutosave<S extends Snapshot>(
	store: Subscribable<S>,
	storage: StorageLike,
	delayMs: number = AUTOSAVE_DELAY_MS,
): () => void {
	let timer: ReturnType<typeof setTimeout> | null = null;
	const unsubscribe = store.subscribe((state, prev) => {
		if (
			state.nodes === prev.nodes &&
			state.edges === prev.edges &&
			state.solverSettings === prev.solverSettings
		) {
			return;
		}
		if (timer !== null) clearTimeout(timer);
		timer = setTimeout(() => {
			timer = null;
			writeAutosave(storage, store.getState());
		}, delayMs);
	});
	return () => {
		if (timer !== null) clearTimeout(timer);
		unsubscribe();
	};
}

export const restoredNote = (savedAt: string): string => {
	const t = new Date(savedAt);
	const hhmm = Number.isNaN(t.getTime())
		? ""
		: ` (autosaved ${t.toLocaleTimeString([], { hour: "2-digit", minute: "2-digit" })})`;
	return `Restored your last session${hhmm}.`;
};
