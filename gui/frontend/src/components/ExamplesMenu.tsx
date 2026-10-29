import { BookOpen, ChevronDown } from "lucide-react";
import { useEffect, useRef, useState } from "react";
import { fetchExample, fetchExamples } from "../api";
import useStore from "../store/useStore";
import type { ExampleListing } from "../utils/exampleInfo";

// Worked examples from papers (#793): networks shipped with the package,
// listed by the backend's /examples. Loading one replaces the canvas like
// Load does, so it asks first when the canvas holds something.
const ExamplesMenu = () => {
	const loadNetwork = useStore((s) => s.loadNetwork);
	const hasWork = useStore((s) => s.nodes.length > 0);
	const [open, setOpen] = useState(false);
	const [items, setItems] = useState<ExampleListing[] | null>(null);
	const [error, setError] = useState<string | null>(null);
	const ref = useRef<HTMLDivElement>(null);

	useEffect(() => {
		if (!open || items !== null) return;
		fetchExamples()
			.then(setItems)
			.catch(() =>
				setError("Could not list the examples: is the backend running?"),
			);
	}, [open, items]);

	useEffect(() => {
		if (!open) return;
		const onDown = (e: MouseEvent) => {
			if (ref.current && !ref.current.contains(e.target as globalThis.Node)) {
				setOpen(false);
			}
		};
		document.addEventListener("mousedown", onDown);
		return () => document.removeEventListener("mousedown", onDown);
	}, [open]);

	const choose = async (ex: ExampleListing) => {
		if (
			hasWork &&
			!window.confirm(`Replace the canvas with the example "${ex.title}"?`)
		) {
			return;
		}
		setOpen(false);
		try {
			loadNetwork(await fetchExample(ex.id));
		} catch {
			alert(`Could not load the example "${ex.title}".`);
		}
	};

	return (
		<div className="relative" ref={ref}>
			<button
				type="button"
				onClick={() => setOpen((o) => !o)}
				aria-haspopup="menu"
				aria-expanded={open}
				className="flex items-center gap-2 px-4 py-1.5 border rounded-md hover:bg-slate-50 transition-colors text-sm font-medium"
			>
				<BookOpen size={16} /> Examples <ChevronDown size={14} />
			</button>
			{open && (
				<div
					role="menu"
					className="absolute right-0 mt-1 w-80 bg-white border rounded-md shadow-lg z-20 py-1"
				>
					<div className="px-3 py-1 text-[10px] font-bold text-gray-500 uppercase tracking-wider">
						Worked examples from papers
					</div>
					{error && (
						<div className="px-3 py-2 text-xs text-red-600">{error}</div>
					)}
					{!error && items === null && (
						<div className="px-3 py-2 text-xs text-gray-500">Loading…</div>
					)}
					{items?.length === 0 && (
						<div className="px-3 py-2 text-xs text-gray-500">
							No examples shipped.
						</div>
					)}
					{items?.map((ex) => (
						<button
							key={ex.id}
							type="button"
							role="menuitem"
							onClick={() => void choose(ex)}
							className="w-full text-left px-3 py-2 text-sm hover:bg-slate-50 flex items-center justify-between gap-2"
						>
							<span>{ex.title}</span>
							{ex.label && (
								<span className="text-[10px] uppercase tracking-wider text-emerald-700">
									{ex.label}
								</span>
							)}
						</button>
					))}
				</div>
			)}
		</div>
	);
};

export default ExamplesMenu;
