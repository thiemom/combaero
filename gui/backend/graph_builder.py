from combaero.network import (
    FlowNetwork,
    LosslessConnectionElement,
    MassFlowBoundary,
    MomentumChamberNode,
    OrificeElement,
    PipeElement,
    PlenumNode,
    PressureBoundary,
)

from .schemas import (
    LosslessConnectionData,
    MassBoundaryData,
    MomentumChamberData,
    NetworkGraphSchema,
    OrificeData,
    PipeData,
    PlenumData,
    PressureBoundaryData,
)


def build_network_from_schema(schema: NetworkGraphSchema) -> FlowNetwork:
    net = FlowNetwork()
    nodes_map = {}

    # 1. Add Nodes
    for node_schema in schema.nodes:
        node_id = node_schema.id
        node_data = node_schema.data

        if isinstance(node_data, PlenumData):
            node = PlenumNode(node_id)
        elif isinstance(node_data, MassBoundaryData):
            node = MassFlowBoundary(
                node_id, m_dot=node_data.m_dot, T_total=node_data.T_total, Y=node_data.Y
            )
        elif isinstance(node_data, PressureBoundaryData):
            node = PressureBoundary(node_id, P_total=node_data.P_total)
        elif isinstance(node_data, MomentumChamberData):
            node = MomentumChamberNode(node_id)
        else:
            raise ValueError(f"Unknown node type: {node_schema.type}")

        net.add_node(node)
        nodes_map[node_id] = node

    # 2. Add Elements (Edges)
    for edge_schema in schema.edges:
        edge_id = edge_schema.id
        source_id = edge_schema.source
        target_id = edge_schema.target
        edge_data = edge_schema.data

        if isinstance(edge_data, PipeData):
            elem = PipeElement(
                edge_id,
                from_node=nodes_map[source_id],
                to_node=nodes_map[target_id],
                L=edge_data.L,
                D=edge_data.D,
                roughness=edge_data.roughness,
            )
        elif isinstance(edge_data, OrificeData):
            elem = OrificeElement(
                edge_id,
                from_node=nodes_map[source_id],
                to_node=nodes_map[target_id],
                area=edge_data.area,
                Cd=edge_data.Cd,
            )
        elif isinstance(edge_data, LosslessConnectionData):
            elem = LosslessConnectionElement(
                edge_id, from_node=nodes_map[source_id], to_node=nodes_map[target_id]
            )
        else:
            raise ValueError(f"Unknown edge type: {edge_schema.type}")

        net.add_element(elem)

    return net
