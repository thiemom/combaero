import useStore from "../store/useStore";

const Inspector = () => {
	const { nodes, edges, updateNodeData } = useStore();

	// Get selected node
	const selectedNode = nodes.find((n) => n.selected);
	const selectedEdge = edges.find((e) => e.selected);

	if (!selectedNode && !selectedEdge) {
		return (
			<aside className="w-80 border-l bg-white p-4 flex flex-col gap-4 italic text-gray-400">
				Select a node or edge to edit properties
			</aside>
		);
	}

	if (selectedNode) {
		return (
			<aside className="w-80 border-l bg-white p-4 flex flex-col gap-4 overflow-y-auto">
				<h2 className="text-lg font-bold border-b pb-2 capitalize">
					{selectedNode.type} Node
				</h2>

				<div className="flex flex-col gap-2">
					<label
						htmlFor={`node_id_${selectedNode.id}`}
						className="text-xs font-bold text-gray-500 uppercase"
					>
						Node ID
					</label>
					<input
						id={`node_id_${selectedNode.id}`}
						type="text"
						value={selectedNode.id}
						className="p-2 border rounded bg-gray-50 cursor-not-allowed"
						disabled
					/>
				</div>

				{selectedNode.type === "mass_boundary" && (
					<>
						<div className="flex flex-col gap-2">
							<label
								htmlFor={`m_dot_${selectedNode.id}`}
								className="text-xs font-bold text-gray-500 uppercase"
							>
								Mass Flow (kg/s)
							</label>
							<input
								id={`m_dot_${selectedNode.id}`}
								type="number"
								value={selectedNode.data.m_dot || 1.0}
								onChange={(e) =>
									updateNodeData(selectedNode.id, {
										m_dot: parseFloat(e.target.value),
									})
								}
								className="p-2 border rounded"
							/>
						</div>
						<div className="flex flex-col gap-2">
							<label
								htmlFor={`T_total_${selectedNode.id}`}
								className="text-xs font-bold text-gray-500 uppercase"
							>
								T_total (K)
							</label>
							<input
								id={`T_total_${selectedNode.id}`}
								type="number"
								value={selectedNode.data.T_total || 300}
								onChange={(e) =>
									updateNodeData(selectedNode.id, {
										T_total: parseFloat(e.target.value),
									})
								}
								className="p-2 border rounded"
							/>
						</div>
					</>
				)}

				{selectedNode.type === "pressure_boundary" && (
					<div className="flex flex-col gap-2">
						<label
							htmlFor={`P_total_${selectedNode.id}`}
							className="text-xs font-bold text-gray-500 uppercase"
						>
							P_total (Pa)
						</label>
						<input
							id={`P_total_${selectedNode.id}`}
							type="number"
							value={selectedNode.data.P_total || 101325}
							onChange={(e) =>
								updateNodeData(selectedNode.id, {
									P_total: parseFloat(e.target.value),
								})
							}
							className="p-2 border rounded"
						/>
					</div>
				)}

				{selectedNode.data.result && (
					<div className="mt-4 p-3 bg-blue-50 border border-blue-200 rounded">
						<h3 className="text-sm font-bold text-blue-700 mb-2">
							Live Telemetry
						</h3>
						<div className="text-xs grid grid-cols-2 gap-y-1">
							<span className="text-blue-600">Pressure:</span>
							<span className="font-mono">
								{(selectedNode.data.result.P / 1e5).toFixed(4)} bar
							</span>
							<span className="text-blue-600">Temperature:</span>
							<span className="font-mono">
								{selectedNode.data.result.T.toFixed(2)} K
							</span>
						</div>
					</div>
				)}
			</aside>
		);
	}

	// Edge Inspector placeholder
	return (
		<aside className="w-80 border-l bg-white p-4 flex flex-col gap-4 italic text-gray-400">
			Edge editing (Pipe/Orifice) implemented in next step
		</aside>
	);
};

export default Inspector;
