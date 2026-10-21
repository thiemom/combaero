import { Handle, Position } from "reactflow";
import { FILM_COLOR, FILM_HANDLE } from "../../utils/filmLink";
import { handleStyle, rotPos } from "../../utils/nodeUtils";

// The "f" handle every film provider carries at 6 o'clock (#386): dragged
// onto a thermal wall's label box it names the node as that wall's film
// (NetworkCanvas); it never makes an edge.
export const FilmHandle = ({ rotation }: { rotation: number }) => (
	<>
		<div
			className="absolute -bottom-3 left-1/2 -translate-x-1/2 text-[7px] font-extrabold leading-none select-none pointer-events-none bg-white/70 rounded-sm px-0.5"
			style={{ transform: `rotate(${-rotation}deg)`, color: FILM_COLOR }}
		>
			f
		</div>
		<Handle
			type="source"
			position={rotPos(Position.Bottom, rotation)}
			style={handleStyle(Position.Bottom, rotation)}
			id={FILM_HANDLE}
			title="Drag onto a thermal wall's label to cover it with this film"
		/>
	</>
);
