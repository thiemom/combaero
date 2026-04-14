from combaero.network import (
    CombustorNode,
    ConvectiveSurface,
    FlowNetwork,
    LosslessConnectionElement,
    MassFlowBoundary,
    MomentumChamberNode,
    OrificeElement,
    NetworkElement,
    PipeElement,
    PlenumNode,
    PressureBoundary,
    ThermalWall,
    WallConnection,
    WallLayer,
)

from .schemas import (
    CombustorData,
    CompositionData,
    MassBoundaryData,
    MomentumChamberData,
    NetworkGraphSchema,
    OrificeData,
    PipeData,
    PlenumData,
    PressureBoundaryData,
    ThermalWallData,
)


def resolve_composition(comp: CompositionData, T: float, P: float) -> list[float]:
    """
    Converts UI CompositionData into a mass fraction list (Y)
    using combaero.species utilities.
    """
    import combaero as cb

    if comp.source == "dry_air":
        Y = cb.species.dry_air_mass()
    elif comp.source == "humid_air":
        # humid_air_mass takes ambient T, P, RH
        Y = cb.species.humid_air_mass(comp.ambient_T, comp.ambient_P, comp.relative_humidity)
    elif comp.source == "fuel":
        # Default fuel to CH4 for now, or use a specific fuel if provided
        Y = cb.species.pure_species("CH4")
        Y = cb.species.to_mass(Y)
    elif comp.source == "custom" and comp.custom_fractions:
        vec = cb.species.from_mapping(comp.custom_fractions)
        Y = cb.species.to_mass(vec) if comp.mode == "mole" else vec
    else:
        Y = cb.species.dry_air_mass()

    return [float(y) for y in Y]


def build_network_from_schema(schema: NetworkGraphSchema) -> FlowNetwork:
    net = FlowNetwork()
    nodes_map = {}  # ID -> Physical Node
    plenum_node_ids: set[str] = set()
    element_nodes = []  # List of node schemas that are actually elements

    # 1. First Pass: Create Physical Nodes
    for node_schema in schema.nodes:
        node_id = node_schema.id
        node_data = node_schema.data
        node_type = node_schema.type

        if node_type == "plenum":
            data = PlenumData(**node_data)
            node = PlenumNode(node_id)
            node.initial_guess = data.initial_guess
            net.add_node(node)
            nodes_map[node_id] = node
            plenum_node_ids.add(node_id)
        elif node_type == "mass_boundary":
            data = MassBoundaryData(**node_data)
            Y = resolve_composition(data.composition, data.T_total, 101325.0)
            node = MassFlowBoundary(node_id, m_dot=data.m_dot, T_total=data.T_total, Y=Y)
            node.initial_guess = data.initial_guess
            net.add_node(node)
            nodes_map[node_id] = node
        elif node_type == "pressure_boundary":
            data = PressureBoundaryData(**node_data)
            Y = resolve_composition(data.composition, data.T_total, data.P_total)
            node = PressureBoundary(node_id, P_total=data.P_total, T_total=data.T_total, Y=Y)
            node.initial_guess = data.initial_guess
            net.add_node(node)
            nodes_map[node_id] = node
        elif node_type == "combustor":
            data = CombustorData(**node_data)
            node = CombustorNode(node_id, method=data.method)
            node.initial_guess = data.initial_guess
            net.add_node(node)
            nodes_map[node_id] = node
        elif node_type == "momentum_chamber":
            data = MomentumChamberData(**node_data)
            node = MomentumChamberNode(
                node_id,
                area=data.area,
                surface=ConvectiveSurface(
                    area=data.area,
                    Nu_multiplier=data.Nu_multiplier,
                    # MomentumChamberNode is modeled as lossless for pressure-flow,
                    # so friction multiplier is fixed to 1.0 by design.
                    f_multiplier=1.0,
                ),
            )
            node.initial_guess = data.initial_guess
            net.add_node(node)
            nodes_map[node_id] = node
        elif node_type == "probe":
            pass  # Display-only diagnostic nodes -- not part of the solver graph
        else:
            # This is an element node (Pipe, Orifice, etc.)
            element_nodes.append(node_schema)

    element_ids = {node.id for node in element_nodes}

    # 1b. Create implicit junction nodes for element -> element visual links.
    edge_junction_map: dict[tuple[str, str], str] = {}
    for edge in schema.edges:
        if edge.data and edge.data.get("type") == "thermal":
            continue
        if edge.source in element_ids and edge.target in element_ids:
            junction_id = f"__junction__{edge.source}__{edge.target}"
            if junction_id not in nodes_map:
                junction_node = MomentumChamberNode(junction_id)
                net.add_node(junction_node)
                nodes_map[junction_id] = junction_node
            edge_junction_map[(edge.source, edge.target)] = junction_id

    # 2. Second Pass: Create Elements
    for elem_node_schema in element_nodes:
        elem_id = elem_node_schema.id
        elem_data = elem_node_schema.data

        source_id = None
        target_id = None
        incoming_count = 0
        outgoing_count = 0

        for edge in schema.edges:
            if edge.data and edge.data.get("type") == "thermal":
                continue
            if edge.target == elem_id:
                incoming_count += 1
                if incoming_count > 1:
                    raise ValueError(f"Element '{elem_id}' has multiple upstream links.")
                if edge.source in nodes_map:
                    source_id = edge.source
                elif edge.source in element_ids:
                    source_id = edge_junction_map.get((edge.source, elem_id))
            if edge.source == elem_id:
                outgoing_count += 1
                if outgoing_count > 1:
                    raise ValueError(f"Element '{elem_id}' has multiple downstream links.")
                if edge.target in nodes_map:
                    target_id = edge.target
                elif edge.target in element_ids:
                    target_id = edge_junction_map.get((elem_id, edge.target))

        if not source_id or not target_id:
            raise ValueError(f"Element '{elem_id}' is not fully connected.")

        elem_type = elem_node_schema.type
        if elem_type == "pipe":
            data = PipeData(**elem_data)
            conv_area = 3.1415926535 * data.D * data.L
            elem = PipeElement(
                elem_id,
                from_node=source_id,
                to_node=target_id,
                length=data.L,
                diameter=data.D,
                roughness=data.roughness,
                surface=ConvectiveSurface(
                    area=conv_area,
                    Nu_multiplier=data.Nu_multiplier,
                    f_multiplier=data.f_multiplier,
                ),
            )
            regime = (
                data.regime if data.regime != "default" else schema.solver_settings.global_regime
            )
            elem.regime = regime
            elem.initial_guess = data.initial_guess
            net.add_element(elem)
        elif elem_type == "orifice":
            data = OrificeData(**elem_data)
            elem = OrificeElement(
                elem_id,
                from_node=source_id,
                to_node=target_id,
                diameter=data.diameter,
                Cd=data.Cd,
                auto_Cd=data.auto_Cd,
                plate_thickness=data.plate_thickness,
                edge_radius=data.edge_radius,
            )
            regime = (
                data.regime if data.regime != "default" else schema.solver_settings.global_regime
            )
            elem.regime = regime
            elem.initial_guess = data.initial_guess
            net.add_element(elem)
        elif elem_type == "lossless_connection":
            elem = LosslessConnectionElement(elem_id, from_node=source_id, to_node=target_id)
            net.add_element(elem)

    # 3. Third Pass: Auto-connections
    for edge in schema.edges:
        if edge.source in nodes_map and edge.target in nodes_map and edge.data is None:
            auto_id = f"__auto_link__{edge.source}__{edge.target}"
            if auto_id not in net.elements:
                elem = LosslessConnectionElement(
                    auto_id, from_node=edge.source, to_node=edge.target
                )
                net.add_element(elem)

    # 4. Fourth Pass: Thermal Walls
    for edge in schema.edges:
        if edge.data and edge.data.get("type") == "thermal":
            if edge.source in plenum_node_ids or edge.target in plenum_node_ids:
                continue
            wall_id = edge.id
            data = ThermalWallData(**edge.data)

            # Build layers list
            layers = []
            if data.layers:
                # New multi-layer format
                layers = [
                    WallLayer(thickness=L.thickness, conductivity=L.conductivity, material=L.material)
                    for L in data.layers
                ]
            else:
                # Legacy scalar format
                t = data.thickness if data.thickness is not None else 0.003
                k = data.conductivity if data.conductivity is not None else 20.0
                layers = [WallLayer(thickness=t, conductivity=k)]

            net.add_wall(
                ThermalWall(
                    id=wall_id,
                    element_a=edge.source,
                    element_b=edge.target,
                    layers=layers,
                    contact_area=data.area,
                    R_fouling=data.R_fouling,
                )
            )
            net.thermal_coupling_enabled = True

    return net
