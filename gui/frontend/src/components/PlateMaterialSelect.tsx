import type React from "react";
import { useEffect, useState } from "react";
import { fetchMaterials } from "../api";

interface PlateMaterialSelectProps {
	value: string | null | undefined;
	onChange: (material: string | null) => void;
}

// An effusion plate's alloy, for its wall temperature limit only (#998): the
// hot end of its k(T) envelope (#882). Unset, the limit is the generic metal
// cap. The plate's conductivity stays its own Wall Conductivity.
const PlateMaterialSelect: React.FC<PlateMaterialSelectProps> = ({
	value,
	onChange,
}) => {
	const [materials, setMaterials] = useState<string[]>([]);
	useEffect(() => {
		fetchMaterials().then(setMaterials).catch(console.error);
	}, []);
	return (
		<div className="flex flex-col gap-1">
			<label className="text-xs font-bold text-gray-500 uppercase">
				Plate Material (T limit)
			</label>
			<select
				className="p-2 border rounded bg-white text-xs border-stone-200"
				value={value ?? ""}
				onChange={(e) => onChange(e.target.value || null)}
			>
				<option value="">Generic metal (1533 K cap)</option>
				{materials.map((m) => (
					<option key={m} value={m}>
						{m}
					</option>
				))}
			</select>
		</div>
	);
};

export default PlateMaterialSelect;
