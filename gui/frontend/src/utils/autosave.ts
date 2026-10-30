/**
 * Autosave of the network to the browser's localStorage (#596).
 *
 * The record is what Save writes ({nodes, edges, solverSettings, version}),
 * without solve results (stale after any edit, and re-solving is cheap), plus
 * the time it was written. One slot per site, shared by every tab: the last
 * tab edited wins.
 *
 * Cross-browser: every storage access is guarded. Reading
 * `window.localStorage` itself throws a SecurityError when site data is
 * blocked (Safari, Firefox, Chrome and Edge all can), setItem throws
 * QuotaExceededError when full, and a private window may offer a storage that
 * is wiped on close. None of that may break the editor, so a failed write
 * returns false and a failed read returns null. The flush on leaving uses
 * `pagehide` and `visibilitychange`: Safari, and every mobile browser, do not
 * reliably fire `beforeunload`.
 */
import type { Edge, Node } from "reactflow";
import { stripResults } from "../store/useStore";
import { NETWORK_FILE_VERSION } from "./nodeOrigin";

export const AUTOSAVE_KEY = "combaero.autosave.v1";
export const AUTOSAVE_DELAY_MS = 1000;

export interface AutosaveRecord {
	nodes: Node[];
	edges: Edge[];
	solverSettings: Record<string, unknown>;
	version: string;
	savedAt: string;
}

export interface AutosaveSnapshot {
	nodes: Node[];
	edges: Edge[];
	solverSettings: object;
}

type StorageGetter = () => Storage | null;

const defaultStorage: StorageGetter = () => {
	try {
		return typeof window !== "undefined" ? window.localStorage : null;
	} catch {
		return null;
	}
};

/** Write the snapshot; an empty network removes the record instead. */
export const writeAutosave = (
	snap: AutosaveSnapshot,
	storage: StorageGetter = defaultStorage,
	now: Date = new Date(),
): boolean => {
	const s = storage();
	if (!s) return false;
	try {
		if (snap.nodes.length === 0 && snap.edges.length === 0) {
			s.removeItem(AUTOSAVE_KEY);
			return true;
		}
		const record: AutosaveRecord = {
			nodes: stripResults(snap.nodes),
			edges: stripResults(snap.edges),
			solverSettings: { ...snap.solverSettings } as Record<string, unknown>,
			version: NETWORK_FILE_VERSION,
			savedAt: now.toISOString(),
		};
		s.setItem(AUTOSAVE_KEY, JSON.stringify(record));
		return true;
	} catch {
		return false;
	}
};

/** The saved record, or null when there is none or it is unreadable. */
export const readAutosave = (
	storage: StorageGetter = defaultStorage,
): AutosaveRecord | null => {
	const s = storage();
	if (!s) return null;
	try {
		const raw = s.getItem(AUTOSAVE_KEY);
		if (!raw) return null;
		const rec = JSON.parse(raw);
		if (
			!rec ||
			!Array.isArray(rec.nodes) ||
			!Array.isArray(rec.edges) ||
			typeof rec.savedAt !== "string"
		) {
			return null;
		}
		return rec as AutosaveRecord;
	} catch {
		return null;
	}
};

export const clearAutosave = (
	storage: StorageGetter = defaultStorage,
): void => {
	const s = storage();
	if (!s) return;
	try {
		s.removeItem(AUTOSAVE_KEY);
	} catch {
		// Blocked storage: nothing was saved there either.
	}
};

interface Subscribable {
	getState: () => AutosaveSnapshot;
	subscribe: (
		listener: (state: AutosaveSnapshot, prev: AutosaveSnapshot) => void,
	) => () => void;
}

/**
 * Write the store to storage AUTOSAVE_DELAY_MS after the last change to its
 * nodes, edges or solver settings, and at once when the page is hidden or
 * left. Returns a stop function.
 */
export const startAutosave = (
	store: Subscribable,
	storage: StorageGetter = defaultStorage,
	delayMs: number = AUTOSAVE_DELAY_MS,
): (() => void) => {
	let timer: ReturnType<typeof setTimeout> | null = null;
	let pending = false;

	const flush = () => {
		if (timer !== null) {
			clearTimeout(timer);
			timer = null;
		}
		if (!pending) return;
		pending = false;
		const { nodes, edges, solverSettings } = store.getState();
		writeAutosave({ nodes, edges, solverSettings }, storage);
	};

	const unsubscribe = store.subscribe((state, prev) => {
		if (
			state.nodes === prev.nodes &&
			state.edges === prev.edges &&
			state.solverSettings === prev.solverSettings
		) {
			return;
		}
		pending = true;
		if (timer !== null) clearTimeout(timer);
		timer = setTimeout(flush, delayMs);
	});

	const onHidden = () => {
		if (document.visibilityState === "hidden") flush();
	};
	const hasWindow = typeof window !== "undefined";
	if (hasWindow) {
		window.addEventListener("pagehide", flush);
		document.addEventListener("visibilitychange", onHidden);
	}

	return () => {
		flush();
		unsubscribe();
		if (hasWindow) {
			window.removeEventListener("pagehide", flush);
			document.removeEventListener("visibilitychange", onHidden);
		}
	};
};
