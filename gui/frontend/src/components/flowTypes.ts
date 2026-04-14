import ChannelNode from "./nodes/ChannelNode";
import CombustorNode from "./nodes/CombustorNode.tsx";
import LosslessNode from "./nodes/LosslessNode";
import MassBoundaryNode from "./nodes/MassBoundaryNode";
import MomentumChamberNode from "./nodes/MomentumChamberNode.tsx";
import OrificeNode from "./nodes/OrificeNode";
import PlenumNode from "./nodes/PlenumNode";
import PressureBoundaryNode from "./nodes/PressureBoundaryNode";
import ProbeNode from "./nodes/ProbeNode";
import ThermalEdge from "./ThermalEdge";

export const nodeTypes = {
	plenum: PlenumNode,
	mass_boundary: MassBoundaryNode,
	pressure_boundary: PressureBoundaryNode,
	channel: ChannelNode,
	orifice: OrificeNode,
	combustor: CombustorNode,
	momentum_chamber: MomentumChamberNode,
	lossless_connection: LosslessNode,
	probe: ProbeNode,
};

export const edgeTypes = {
	thermal: ThermalEdge,
};
