import { BookOpen, X } from "lucide-react";
import useStore from "../store/useStore";
import { labelText } from "../utils/exampleInfo";

// The description of a loaded worked example (#793): the paper, what is
// reproduced, where to look, and the validation label. Dismissible; a Clear
// or another load replaces it.
const ExamplePanel = () => {
	const info = useStore((s) => s.exampleInfo);
	const dismiss = useStore((s) => s.dismissExampleInfo);
	if (!info) return null;
	const fidelity = info.label === "fidelity";
	return (
		<section
			aria-label="Example description"
			className="px-6 py-3 bg-sky-50 border-b border-sky-200 text-xs text-slate-700 flex gap-4"
		>
			<BookOpen size={18} className="text-sky-600 shrink-0 mt-0.5" />
			<div className="flex-1 grid gap-1.5 max-w-5xl">
				<div className="flex items-center gap-3 flex-wrap">
					<span className="text-sm font-bold text-slate-900">{info.title}</span>
					<span
						className={`px-2 py-0.5 rounded-full text-[10px] font-bold uppercase tracking-wider ${
							fidelity
								? "bg-emerald-100 text-emerald-800"
								: "bg-amber-100 text-amber-800"
						}`}
					>
						{labelText(info.label)}
					</span>
				</div>
				{info.source && (
					<p>
						<span className="font-semibold">Source. </span>
						{info.source}
					</p>
				)}
				{info.reproduced && (
					<p>
						<span className="font-semibold">What is reproduced. </span>
						{info.reproduced}
					</p>
				)}
				{info.whereToLook && (
					<p>
						<span className="font-semibold">Where to look. </span>
						{info.whereToLook}
					</p>
				)}
				{info.labelNote && (
					<p>
						<span className="font-semibold">Validation. </span>
						{info.labelNote}
					</p>
				)}
				{info.script && (
					<p className="text-slate-500">
						The full comparison against the published values:{" "}
						<code className="font-mono">{info.script}</code>
					</p>
				)}
			</div>
			<button
				type="button"
				onClick={dismiss}
				className="self-start p-1 rounded hover:bg-sky-100"
				aria-label="Dismiss example description"
			>
				<X size={14} />
			</button>
		</section>
	);
};

export default ExamplePanel;
