from combaero.network import (
    ChannelElement,
    CombustorNode,
    ConstantFractionLoss,
    ConstantHeadLoss,
    ConvectiveSurface,
    DimpledModel,
    FlowNetwork,
    ImpingementModel,
    LinearThetaFractionLoss,
    LinearThetaHeadLoss,
    LosslessConnectionElement,
    MassFlowBoundary,
    MomentumChamberNode,
    OrificeElement,
    PinFinModel,
    PlenumNode,
    PressureBoundary,
    PressureLossElement,
    RibbedModel,
    SmoothModel,
    ThermalWall,
    WallLayer,
)

from .schemas import (
    ChannelData,
    CombustorData,
    CompositionData,
    ConstantFractionLossData,
    ConstantHeadLossData,
    DimpledModelData,
    DiscreteLossData,
    ImpingementModelData,
    LinearThetaFractionLossData,
    LinearThetaHeadLossData,
    MassBoundaryData,
    MomentumChamberData,
    NetworkGraphSchema,
    OrificeData,
    PinFinModelData,
    PlenumData,
    PressureBoundaryData,
    RibbedModelData,
    SmoothModelData,
    ThermalWallData,
)


def build_pressure_loss(data, area: float = 0.1):
    """Converts PressureLossData schema to a pressure-loss callable."""
    if data is None:
        return None
    t = getattr(data, "type", None)
    if isinstance(data, ConstantFractionLossData) or t == "constant_fraction":
        return ConstantFractionLoss(xi=data.xi)
    if isinstance(data, LinearThetaFractionLossData) or t == "linear_theta_fraction":
        return LinearThetaFractionLoss(k=data.k, xi0=data.xi0)
    if isinstance(data, ConstantHeadLossData) or t == "constant_head":
        return ConstantHeadLoss(zeta=data.zeta, area=area)
    if isinstance(data, LinearThetaHeadLossData) or t == "linear_theta_head":
        return LinearThetaHeadLoss(k=data.k, zeta0=data.zeta0, area=area)
    return None


def map_surface_model(data):
    """Maps UI SurfaceModelData to combaero.network models."""
    if isinstance(data, SmoothModelData) or data.type == "smooth":
        return SmoothModel()
    elif isinstance(data, RibbedModelData) or data.type == "ribbed":
        return RibbedModel(
            e_D=data.e_D,
            pitch_to_height=data.pitch_to_height,
            alpha_deg=data.alpha_deg,
        )
    elif isinstance(data, DimpledModelData) or data.type == "dimpled":
        return DimpledModel(
            d_Dh=data.d_Dh,
            h_d=data.h_d,
            S_d=data.S_d,
        )
    elif isinstance(data, PinFinModelData) or data.type == "pin_fin":
        return PinFinModel(
            pin_diameter=data.pin_diameter,
            channel_height=data.channel_height,
            S_D=data.S_D,
            X_D=data.X_D,
            N_rows=data.N_rows,
            is_staggered=data.is_staggered,
        )
    elif isinstance(data, ImpingementModelData) or data.type == "impingement":
        return ImpingementModel(
            d_jet=data.d_jet,
            z_D=data.z_D,
            x_D=data.x_D,
            y_D=data.y_D,
            A_target=data.A_target,
            Cd_jet=data.Cd_jet,
        )
    return SmoothModel()


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


def _expand_initial_guess(short_guess: dict, prefix: str) -> dict:
    """Translate short GUI guess keys (``P``, ``T``, ``m_dot``) into the fully
    qualified solver unknown names (e.g. ``"plenum1.P"``). ``P`` and ``T`` are
    broadcast to both static and total variants because the solver exposes
    them as separate unknowns on some node types; the solver silently ignores
    entries that do not correspond to an actual unknown.
    Any key already containing a ``.`` is passed through unchanged so that
    hand-authored qualified guesses still work. Explicit short keys
    (e.g. ``P_total``) take precedence over broadcast values from ``P``/``T``.
    """
    if not short_guess:
        return {}
    result: dict = {}
    # Pass 1 — broadcast P and T to both static and total variants.
    for key, val in short_guess.items():
        if "." in key:
            continue
        if key == "P":
            result[f"{prefix}.P"] = val
            result[f"{prefix}.P_total"] = val
        elif key == "T":
            result[f"{prefix}.T"] = val
            result[f"{prefix}.T_total"] = val
    # Pass 2 — explicit short keys (P_total, T_total, m_dot, X[i], etc.)
    # are applied after the broadcast so they win on conflict.
    for key, val in short_guess.items():
        if "." in key:
            result[key] = val
            continue
        if key in ("P", "T"):
            continue  # already handled in pass 1
        result[f"{prefix}.{key}"] = val
    return result


def _build_discrete_loss_correlation(
    corr_type: str,
    xi: float,
    k: float,
    xi0: float,
    zeta: float,
    zeta0: float,
    area: float,
):
    """Build a pressure loss correlation object from discrete loss params."""
    if corr_type == "constant_fraction":
        return ConstantFractionLoss(xi=xi)
    if corr_type == "linear_theta_fraction":
        return LinearThetaFractionLoss(k=k, xi0=xi0)
    if corr_type == "constant_head":
        return ConstantHeadLoss(zeta=zeta, area=area)
    if corr_type == "linear_theta_head":
        return LinearThetaHeadLoss(k=k, zeta0=zeta0, area=area)
    return ConstantFractionLoss(xi=xi)


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
            node.initial_guess = _expand_initial_guess(data.initial_guess, node_id)
            net.add_node(node)
            nodes_map[node_id] = node
            plenum_node_ids.add(node_id)
        elif node_type == "mass_boundary":
            data = MassBoundaryData(**node_data)
            Y = resolve_composition(data.composition, data.T_total, 101325.0)
            node = MassFlowBoundary(node_id, m_dot=data.m_dot, T_total=data.T_total, Y=Y)
            node.initial_guess = _expand_initial_guess(data.initial_guess, node_id)
            net.add_node(node)
            nodes_map[node_id] = node
        elif node_type == "pressure_boundary":
            data = PressureBoundaryData(**node_data)
            Y = resolve_composition(data.composition, data.T_total, data.P_total)
            node = PressureBoundary(node_id, P_total=data.P_total, T_total=data.T_total, Y=Y)
            node.initial_guess = _expand_initial_guess(data.initial_guess, node_id)
            net.add_node(node)
            nodes_map[node_id] = node
        elif node_type == "combustor":
            data = CombustorData(**node_data)
            node = CombustorNode(
                node_id,
                method=data.method,
                area=data.area,
                Dh=data.Dh,
                surface=ConvectiveSurface(
                    area=data.area,
                    model=map_surface_model(data.surface),
                    Nu_multiplier=data.Nu_multiplier,
                    f_multiplier=1.0,
                ),
            )
            node.initial_guess = _expand_initial_guess(data.initial_guess, node_id)
            net.add_node(node)
            nodes_map[node_id] = node
        elif node_type == "momentum_chamber":
            data = MomentumChamberData(**node_data)
            node = MomentumChamberNode(
                node_id,
                area=data.area,
                surface=ConvectiveSurface(
                    area=data.area,
                    model=map_surface_model(data.surface),
                    Nu_multiplier=data.Nu_multiplier,
                    # MomentumChamberNode is modeled as lossless for pressure-flow,
                    # so friction multiplier is fixed to 1.0 by design.
                    f_multiplier=1.0,
                ),
            )
            node.initial_guess = _expand_initial_guess(data.initial_guess, node_id)
            net.add_node(node)
            nodes_map[node_id] = node
        elif node_type == "probe":
            pass  # Display-only diagnostic nodes -- not part of the solver graph
        else:
            # This is an element node (Channel, Orifice, etc.)
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
        if elem_type == "channel":
            data = ChannelData(**elem_data)
            conv_area = 3.1415926535 * data.D * data.L
            elem = ChannelElement(
                elem_id,
                from_node=source_id,
                to_node=target_id,
                length=data.L,
                diameter=data.D,
                roughness=data.roughness,
                surface=ConvectiveSurface(
                    area=conv_area,
                    model=map_surface_model(data.surface),
                    Nu_multiplier=data.Nu_multiplier,
                    f_multiplier=data.f_multiplier,
                ),
            )
            regime = (
                data.regime if data.regime != "default" else schema.solver_settings.global_regime
            )
            elem.regime = regime
            elem.initial_guess = _expand_initial_guess(data.initial_guess, elem_id)
            net.add_element(elem)
        elif elem_type == "orifice":
            data = OrificeData(**elem_data)
            elem = OrificeElement(
                elem_id,
                from_node=source_id,
                to_node=target_id,
                diameter=data.diameter,
                Cd=data.Cd,
                correlation=data.correlation,
                plate_thickness=data.plate_thickness,
                edge_radius=data.edge_radius,
            )
            regime = (
                data.regime if data.regime != "default" else schema.solver_settings.global_regime
            )
            elem.regime = regime
            elem.initial_guess = _expand_initial_guess(data.initial_guess, elem_id)
            net.add_element(elem)
        elif elem_type == "lossless_connection":
            elem = LosslessConnectionElement(elem_id, from_node=source_id, to_node=target_id)
            net.add_element(elem)
        elif elem_type == "discrete_loss":
            data = DiscreteLossData(**elem_data)
            # Infer area from upstream node if not explicitly set
            area = data.area
            if area is None and source_id in nodes_map:
                area = getattr(nodes_map[source_id], "area", 0.1)
            area = area if area and area > 0 else 0.1
            correlation = _build_discrete_loss_correlation(
                data.correlation_type, data.xi, data.k, data.xi0, data.zeta, data.zeta0, area
            )
            theta_source = data.theta_source if data.theta_source not in (None, "none") else None
            surface = ConvectiveSurface(
                area=area,
                model=map_surface_model(data.surface),
                Nu_multiplier=data.Nu_multiplier,
                f_multiplier=data.f_multiplier,
            )
            elem = PressureLossElement(
                elem_id,
                from_node=source_id,
                to_node=target_id,
                correlation=correlation,
                theta_source=theta_source,
                area=area,
                surface=surface,
            )
            elem.initial_guess = _expand_initial_guess(data.initial_guess, elem_id)
            net.add_element(elem)

    # 3. Third Pass: Auto-connections
    for edge in schema.edges:
        if (
            edge.source in nodes_map
            and edge.target in nodes_map
            and not (edge.data or {}).get("type")
        ):
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
                    WallLayer(
                        thickness=L.thickness, conductivity=L.conductivity, material=L.material
                    )
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
