import FlowEdge from "./FlowEdge";
import AreaChangeNode from "./nodes/AreaChangeNode";
import ChannelNode from "./nodes/ChannelNode";
import CombustorNode from "./nodes/CombustorNode.tsx";
import DiscreteLossNode from "./nodes/DiscreteLossNode";
import EffusionPlateNode from "./nodes/EffusionPlateNode";
import EjectorNode from "./nodes/EjectorNode";
import ImpingementArrayNode from "./nodes/ImpingementArrayNode";
import LosslessNode from "./nodes/LosslessNode";
import MassBoundaryNode from "./nodes/MassBoundaryNode";
import MomentumChamberNode from "./nodes/MomentumChamberNode.tsx";
import MPCETeeNode from "./nodes/MPCETeeNode";
import OrificeNode from "./nodes/OrificeNode";
import PlenumNode from "./nodes/PlenumNode";
import PressureBoundaryNode from "./nodes/PressureBoundaryNode";
import ProbeNode from "./nodes/ProbeNode";
import VortexNode from "./nodes/VortexNode";
import WallNode from "./nodes/WallNode";
import ThermalEdge from "./ThermalEdge";

export const nodeTypes = {
	plenum: PlenumNode,
	mass_boundary: MassBoundaryNode,
	pressure_boundary: PressureBoundaryNode,
	channel: ChannelNode,
	orifice: OrificeNode,
	impingement_array: ImpingementArrayNode,
	effusion_plate: EffusionPlateNode,
	combustor: CombustorNode,
	momentum_chamber: MomentumChamberNode,
	discrete_loss: DiscreteLossNode,
	lossless_connection: LosslessNode,
	probe: ProbeNode,
	area_change: AreaChangeNode,
	mpce_tee: MPCETeeNode,
	ejector: EjectorNode,
	vortex: VortexNode,
	wall: WallNode,
};

export const edgeTypes = {
	thermal: ThermalEdge,
	default: FlowEdge,
};
