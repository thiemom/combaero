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

				{selectedNode.type === "pipe" && (
					<>
						<div className="flex flex-col gap-2">
							<label
								htmlFor={`L_${selectedNode.id}`}
								className="text-xs font-bold text-gray-500 uppercase"
							>
								Length (m)
							</label>
							<input
								id={`L_${selectedNode.id}`}
								type="number"
								step="0.1"
								value={selectedNode.data.L || 1.0}
								onChange={(e) =>
									updateNodeData(selectedNode.id, {
										L: parseFloat(e.target.value),
									})
								}
								className="p-2 border rounded"
							/>
						</div>
						<div className="flex flex-col gap-2">
							<label
								htmlFor={`D_${selectedNode.id}`}
								className="text-xs font-bold text-gray-500 uppercase"
							>
								Diameter (m)
							</label>
							<input
								id={`D_${selectedNode.id}`}
								type="number"
								step="0.01"
								value={selectedNode.data.D || 0.1}
								onChange={(e) =>
									updateNodeData(selectedNode.id, {
										D: parseFloat(e.target.value),
									})
								}
								className="p-2 border rounded"
							/>
						</div>
						<div className="flex flex-col gap-2">
							<label
								htmlFor={`roughness_${selectedNode.id}`}
								className="text-xs font-bold text-gray-500 uppercase"
							>
								Roughness (m)
							</label>
							<input
								id={`roughness_${selectedNode.id}`}
								type="number"
								step="1e-6"
								value={selectedNode.data.roughness || 1e-5}
								onChange={(e) =>
									updateNodeData(selectedNode.id, {
										roughness: parseFloat(e.target.value),
									})
								}
								className="p-2 border rounded"
							/>
						</div>
					</>
				)}

				{selectedNode.type === "orifice" && (
					<>
						<div className="flex flex-col gap-2">
							<label
								htmlFor={`area_${selectedNode.id}`}
								className="text-xs font-bold text-gray-500 uppercase"
							>
								Area (m²)
							</label>
							<input
								id={`area_${selectedNode.id}`}
								type="number"
								step="0.001"
								value={selectedNode.data.area || 0.01}
								onChange={(e) =>
									updateNodeData(selectedNode.id, {
										area: parseFloat(e.target.value),
									})
								}
								className="p-2 border rounded"
							/>
						</div>
						<div className="flex flex-col gap-2">
							<label
								htmlFor={`Cd_${selectedNode.id}`}
								className="text-xs font-bold text-gray-500 uppercase"
							>
								Discharge Coeff (Cd)
							</label>
							<input
								id={`Cd_${selectedNode.id}`}
								type="number"
								step="0.01"
								value={selectedNode.data.Cd || 0.6}
								onChange={(e) =>
									updateNodeData(selectedNode.id, {
										Cd: parseFloat(e.target.value),
									})
								}
								className="p-2 border rounded"
							/>
						</div>
					</>
				)}

				{selectedNode.data.result && (
					<div className="mt-4 p-3 bg-blue-50 border border-blue-200 rounded">
						<h3 className="text-sm font-bold text-blue-700 mb-2">
							Live Telemetry
						</h3>
						{selectedNode.data.result.P !== undefined ? (
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
						) : (
							<div className="text-xs grid grid-cols-2 gap-y-1">
								<span className="text-blue-600">Mass Flow:</span>
								<span className="font-mono">
									{selectedNode.data.result.m_dot?.toFixed(4)} kg/s
								</span>
							</div>
						)}
					</div>
				)}
			</aside>
		);
	}

	// Edge Inspector
	return (
		<aside className="w-80 border-l bg-white p-4 flex flex-col gap-4 italic text-gray-400">
			Connection is structural. Topographic link only.
		</aside>
	);
};

export default Inspector;
