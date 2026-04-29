/**
 * DiagnosticsPanel — sidebar widget for network-wide diagnostic discovery.
 *
 * After a solve the panel scans every finite numeric scalar across all node,
 * element, and thermal-wall results and presents a searchable, checkboxed
 * list. Checking a field adds it to `displaySettings`; unchecking removes it.
 *
 * "Filter and forget": the search input resets to empty whenever the solve
 * results change, so the list is always presented clean for the new topology.
 */

import { ChevronDown, ChevronRight, Search, X } from "lucide-react";
import { useEffect, useMemo, useRef, useState } from "react";
import useStore from "../store/useStore";
import { discoverFields } from "../utils/diagnostics";
import { QUANTITY_CATALOGUE } from "../utils/quantities";

function SourceBadge({
	label,
	title,
	color,
}: {
	label: string;
	title: string;
	color: string;
}) {
	return (
		<span
			title={title}
			className={`text-[8px] font-bold leading-none ${color} select-none`}
		>
			{label}
		</span>
	);
}

export function DiagnosticsPanel() {
	const { solveResults, displaySettings, setDisplaySettings } = useStore();
	const [search, setSearch] = useState("");
	const [collapsed, setCollapsed] = useState(false);
	const prevResultsRef = useRef<unknown>(null);

	// "Filter and forget": search text resets on every new solve result.
	useEffect(() => {
		if (solveResults !== prevResultsRef.current) {
			prevResultsRef.current = solveResults;
			setSearch("");
		}
	}, [solveResults]);

	const allFields = useMemo(() => discoverFields(solveResults), [solveResults]);

	const visible = useMemo(() => {
		const q = search.trim().toLowerCase();
		if (!q) return allFields;
		return allFields.filter((f) => {
			if (f.key.toLowerCase().includes(q)) return true;
			const cat = QUANTITY_CATALOGUE[f.key];
			return (
				cat?.label.toLowerCase().includes(q) ||
				cat?.unit.toLowerCase().includes(q)
			);
		});
	}, [allFields, search]);

	const selectedSet = useMemo(
		() => new Set(displaySettings),
		[displaySettings],
	);

	const toggle = (key: string) => {
		setDisplaySettings(
			selectedSet.has(key)
				? displaySettings.filter((k) => k !== key)
				: [...displaySettings, key],
		);
	};

	const selectVisible = () =>
		setDisplaySettings([
			...new Set([...displaySettings, ...visible.map((f) => f.key)]),
		]);

	const clearVisible = () => {
		const toRemove = new Set(visible.map((f) => f.key));
		setDisplaySettings(displaySettings.filter((k) => !toRemove.has(k)));
	};

	// Count selected keys that actually exist in current results.
	const activeCount = useMemo(
		() =>
			displaySettings.filter((k) => allFields.some((f) => f.key === k)).length,
		[displaySettings, allFields],
	);

	return (
		<div className="flex flex-col gap-2">
			{/* Section header — collapsible */}
			<button
				type="button"
				onClick={() => setCollapsed((c) => !c)}
				className="flex items-center justify-between w-full text-left group"
			>
				<span className="text-sm font-bold text-gray-500 uppercase tracking-wider">
					Display Fields
				</span>
				<span className="text-stone-400 group-hover:text-stone-600">
					{collapsed ? <ChevronRight size={14} /> : <ChevronDown size={14} />}
				</span>
			</button>

			{!collapsed &&
				(allFields.length === 0 ? (
					<p className="text-[10px] text-stone-400 italic text-center py-1">
						Run a solve to discover available fields
					</p>
				) : (
					<>
						{/* Count + bulk action buttons */}
						<div className="flex items-center justify-between">
							<span className="text-[10px] text-stone-400">
								{activeCount} / {allFields.length} shown
							</span>
							<div className="flex gap-1">
								<button
									type="button"
									onClick={clearVisible}
									className="text-[9px] font-bold uppercase tracking-tight px-1.5 py-0.5 rounded border border-stone-200 hover:bg-stone-50 text-stone-400"
								>
									None
								</button>
								<button
									type="button"
									onClick={selectVisible}
									className="text-[9px] font-bold uppercase tracking-tight px-1.5 py-0.5 rounded border border-stone-200 hover:bg-stone-50 text-stone-400"
								>
									All
								</button>
							</div>
						</div>

						{/* Search input */}
						<div className="relative">
							<Search
								size={11}
								className="absolute left-2 top-1/2 -translate-y-1/2 text-stone-300 pointer-events-none"
							/>
							<input
								type="text"
								value={search}
								onChange={(e) => setSearch(e.target.value)}
								placeholder="Filter fields…"
								className="w-full pl-6 pr-6 py-1 text-xs border rounded bg-white outline-none focus:ring-1 focus:ring-stone-200 placeholder:text-stone-300"
							/>
							{search && (
								<button
									type="button"
									onClick={() => setSearch("")}
									className="absolute right-2 top-1/2 -translate-y-1/2 text-stone-300 hover:text-stone-500"
								>
									<X size={11} />
								</button>
							)}
						</div>

						{/* Field list */}
						<div className="flex flex-col gap-px max-h-52 overflow-y-auto -mx-1 px-1">
							{visible.length === 0 ? (
								<p className="text-[10px] text-stone-400 italic text-center py-2">
									No fields match "{search}"
								</p>
							) : (
								visible.map((field) => {
									const cat = QUANTITY_CATALOGUE[field.key];
									const isChecked = selectedSet.has(field.key);
									return (
										<label
											key={field.key}
											className="flex items-center gap-1.5 cursor-pointer py-0.5 px-1 rounded hover:bg-stone-50"
										>
											<input
												type="checkbox"
												checked={isChecked}
												onChange={() => toggle(field.key)}
												className="accent-orange-500 w-3 h-3 flex-shrink-0"
											/>
											{/* Label / key */}
											<span className="font-mono text-[11px] text-stone-700 w-[4.5rem] flex-shrink-0 truncate">
												{cat?.label ?? field.key}
											</span>
											{/* Raw key + unit for non-catalogue fields */}
											<span className="text-[10px] text-stone-400 flex-1 min-w-0 truncate">
												{cat ? field.key : ""}
											</span>
											<span className="text-[9px] text-stone-300 flex-shrink-0 w-10 text-right truncate">
												{cat?.unit ?? ""}
											</span>
											{/* Source badges: N = node, E = element, W = wall */}
											<div className="flex gap-px flex-shrink-0">
												{field.inNodes && (
													<SourceBadge
														label="N"
														title="Present in node results"
														color="text-indigo-300"
													/>
												)}
												{field.inElements && (
													<SourceBadge
														label="E"
														title="Present in element results"
														color="text-amber-400"
													/>
												)}
												{field.inEdges && (
													<SourceBadge
														label="W"
														title="Present in thermal wall results"
														color="text-orange-400"
													/>
												)}
											</div>
										</label>
									);
								})
							)}
						</div>
					</>
				))}
		</div>
	);
}
