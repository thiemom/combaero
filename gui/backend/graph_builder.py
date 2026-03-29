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
    nodes_map = {}  # ID -> Physical Node
    element_nodes = []  # List of node schemas that are actually elements

    # 1. First Pass: Create Physical Nodes
    for node_schema in schema.nodes:
        node_id = node_schema.id
        node_data = node_schema.data

        node_type = node_schema.type

        if node_type == "plenum":
            node = PlenumNode(node_id)
            net.add_node(node)
            nodes_map[node_id] = node
        elif node_type == "mass_boundary":
            data = MassBoundaryData(**node_data)
            node = MassFlowBoundary(
                node_id, m_dot=data.m_dot, T_total=data.T_total, Y=data.Y
            )
            net.add_node(node)
            nodes_map[node_id] = node
        elif node_type == "pressure_boundary":
            data = PressureBoundaryData(**node_data)
            node = PressureBoundary(node_id, P_total=data.P_total)
            net.add_node(node)
            nodes_map[node_id] = node
        elif node_type == "momentum_chamber":
            node = MomentumChamberNode(node_id)
            net.add_node(node)
            nodes_map[node_id] = node
        else:
            # This is an element node (Pipe, Orifice, etc.)
            element_nodes.append(node_schema)

    element_ids = {node.id for node in element_nodes}

    # 1b. Create implicit junction nodes for element -> element visual links.
    # This allows series chains like:
    # mass_boundary -> pipe -> orifice -> pressure_boundary
    # by inserting an internal physical node between adjacent elements.
    edge_junction_map: dict[tuple[str, str], str] = {}
    for edge in schema.edges:
        if edge.source in element_ids and edge.target in element_ids:
            junction_id = f"__junction__{edge.source}__{edge.target}"
            if junction_id not in nodes_map:
                junction_node = MomentumChamberNode(junction_id)
                net.add_node(junction_node)
                nodes_map[junction_id] = junction_node
            edge_junction_map[(edge.source, edge.target)] = junction_id

    # 2. Second Pass: Create Elements by tracing visual edges
    # We need to find connections for each element node
    for elem_node_schema in element_nodes:
        elem_id = elem_node_schema.id
        elem_data = elem_node_schema.data

        # Find visual edges connected to this element node
        # Incoming edge: Source Node -> Element Node
        # Outgoing edge: Element Node -> Target Node
        source_id = None
        target_id = None
        incoming_count = 0
        outgoing_count = 0

        for edge in schema.edges:
            if edge.target == elem_id:
                incoming_count += 1
                if incoming_count > 1:
                    raise ValueError(
                        f"Element '{elem_id}' has multiple upstream links. "
                        "Each element must have exactly one upstream connection."
                    )

                if edge.source in nodes_map:
                    source_id = edge.source
                elif edge.source in element_ids:
                    source_id = edge_junction_map.get((edge.source, elem_id))
                else:
                    source_id = None

            if edge.source == elem_id:
                outgoing_count += 1
                if outgoing_count > 1:
                    raise ValueError(
                        f"Element '{elem_id}' has multiple downstream links. "
                        "Each element must have exactly one downstream connection."
                    )

                if edge.target in nodes_map:
                    target_id = edge.target
                elif edge.target in element_ids:
                    target_id = edge_junction_map.get((elem_id, edge.target))
                else:
                    target_id = None

        if not source_id or not target_id:
            raise ValueError(
                f"Element '{elem_id}' is not fully connected. "
                "Each element must have exactly one upstream and one downstream node."
            )

        if source_id not in nodes_map or target_id not in nodes_map:
            raise ValueError(
                f"Element '{elem_id}' must connect between physical nodes. "
                f"Got source='{source_id}', target='{target_id}'."
            )

        elem_type = elem_node_schema.type

        if elem_type == "pipe":
            data = PipeData(**elem_data)
            elem = PipeElement(
                elem_id,
                from_node=source_id,
                to_node=target_id,
                length=data.L,
                diameter=data.D,
                roughness=data.roughness,
            )
            net.add_element(elem)
        elif elem_type == "orifice":
            data = OrificeData(**elem_data)
            elem = OrificeElement(
                elem_id,
                from_node=source_id,
                to_node=target_id,
                area=data.area,
                Cd=data.Cd,
            )
            net.add_element(elem)
        elif elem_type == "lossless_connection":
            elem = LosslessConnectionElement(elem_id, from_node=source_id, to_node=target_id)
            net.add_element(elem)

    return net
