import {
	ArrowRight,
	ChevronRight,
	Database,
	LogIn,
	LogOut,
} from "lucide-react";
import type React from "react";

const Sidebar = () => {
	const onDragStart = (event: React.DragEvent, nodeType: string) => {
		event.dataTransfer.setData("application/reactflow", nodeType);
		event.dataTransfer.effectAllowed = "move";
	};

	return (
		<aside className="w-64 border-r bg-white p-4 flex flex-col gap-4">
			<div className="text-sm font-bold text-gray-500 uppercase tracking-wider">
				Nodes
			</div>

			<button
				type="button"
				className="flex items-center gap-2 p-2 border rounded cursor-grab hover:bg-stone-50 transition-colors w-full text-left bg-white"
				onDragStart={(event) => onDragStart(event, "mass_boundary")}
				draggable
			>
				<LogIn size={18} className="text-orange-500" />
				<span>Mass Boundary</span>
			</button>

			<button
				type="button"
				className="flex items-center gap-2 p-2 border rounded cursor-grab hover:bg-stone-50 transition-colors w-full text-left bg-white"
				onDragStart={(event) => onDragStart(event, "plenum")}
				draggable
			>
				<Database size={18} className="text-stone-500" />
				<span>Plenum</span>
			</button>

			<button
				type="button"
				className="flex items-center gap-2 p-2 border rounded cursor-grab hover:bg-stone-50 transition-colors w-full text-left bg-white"
				onDragStart={(event) => onDragStart(event, "pressure_boundary")}
				draggable
			>
				<LogOut size={18} className="text-blue-500" />
				<span>Pressure Boundary</span>
			</button>

			<div className="mt-4 text-sm font-bold text-gray-500 uppercase tracking-wider">
				Elements
			</div>

			<button
				type="button"
				className="flex items-center gap-2 p-2 border rounded cursor-grab hover:bg-stone-50 transition-colors w-full text-left bg-white"
				onDragStart={(event) => onDragStart(event, "pipe")}
				draggable
			>
				<ArrowRight size={18} className="text-stone-500" />
				<span>Pipe</span>
			</button>

			<button
				type="button"
				className="flex items-center gap-2 p-2 border rounded cursor-grab hover:bg-stone-50 transition-colors w-full text-left bg-white"
				onDragStart={(event) => onDragStart(event, "orifice")}
				draggable
			>
				<ChevronRight size={18} className="text-orange-400" />
				<span>Orifice</span>
			</button>
		</aside>
	);
};

export default Sidebar;
