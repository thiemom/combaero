import math as _math

from combaero.network import (
    AreaChangeElement,
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
    TeeJunctionElement,
    ThermalWall,
    WallLayer,
)

from .schemas import (
    AreaChangeData,
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
    TeeJunctionData,
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
    (e.g. ``Pt``) take precedence over broadcast values from ``P``/``T``.
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
            result[f"{prefix}.Pt"] = val
        elif key == "T":
            result[f"{prefix}.T"] = val
            result[f"{prefix}.Tt"] = val
    # Pass 2 — explicit short keys (Pt, Tt, m_dot, X[i], etc.)
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


def _port_from_handle(handle: str | None) -> str | None:
    """Extract port name from a tee handle string.

    'port-common-source' -> 'common', 'port-branch-target' -> 'branch', etc.
    Returns None for non-tee handles.
    """
    if (
        handle
        and handle.startswith("port-")
        and (handle.endswith("-source") or handle.endswith("-target"))
    ):
        return handle.split("-")[1]
    return None


def _find_tee_port(
    edges: list,
    tee_id: str,
    port: str,
    nodes_map: dict,
    tee_port_node_map: dict | None = None,
) -> str:
    """Return the physical-node ID connected to a named tee port.

    Checks both edge directions: tee as source (sourceHandle == port-X-source)
    and tee as target (targetHandle == port-X-target).
    Falls back to *tee_port_node_map* for auto-created junction nodes when an
    element is connected directly to a tee port without an explicit plenum.
    """
    src_handle = f"port-{port}-source"
    tgt_handle = f"port-{port}-target"
    for edge in edges:
        if edge.data and edge.data.get("type") == "thermal":
            continue
        if edge.source == tee_id and edge.sourceHandle == src_handle:
            nid = edge.target
        elif edge.target == tee_id and edge.targetHandle == tgt_handle:
            nid = edge.source
        else:
            continue
        if nid in nodes_map:
            return nid
        if tee_port_node_map and (tee_id, port) in tee_port_node_map:
            return tee_port_node_map[(tee_id, port)]
        raise ValueError(
            f"Tee '{tee_id}' port '{port}' connects to '{nid}', which is not a "
            "physical node (plenum or boundary). All tee ports must connect directly "
            "to plena or boundaries."
        )
    raise ValueError(
        f"Tee '{tee_id}' port '{port}' is not connected. "
        "All three tee ports (common, straight, branch) must be wired to a node."
    )


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
            Y = resolve_composition(data.composition, data.Tt, 101325.0)
            node = MassFlowBoundary(node_id, m_dot=data.m_dot, Tt=data.Tt, Y=Y)
            node.initial_guess = _expand_initial_guess(data.initial_guess, node_id)
            net.add_node(node)
            nodes_map[node_id] = node
        elif node_type == "pressure_boundary":
            data = PressureBoundaryData(**node_data)
            Y = resolve_composition(data.composition, data.Tt, data.Pt)
            node = PressureBoundary(node_id, Pt=data.Pt, Tt=data.Tt, Y=Y)
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
                    f_multiplier=data.f_multiplier,
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
                    f_multiplier=data.f_multiplier,
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
    tee_ids = {node.id for node in element_nodes if node.type == "tee_junction"}

    # 1b. Create implicit junction nodes for element -> element visual links.
    # For direct element<->element connections we auto-insert a MomentumChamberNode.
    # For direct element<->tee-port connections we auto-insert a PlenumNode so that
    # users don't need to manually place a plenum on every tee arm.
    edge_junction_map: dict[tuple[str, str], str] = {}
    tee_port_node_map: dict[tuple[str, str], str] = {}  # (tee_id, port) -> node_id
    for edge in schema.edges:
        if edge.data and edge.data.get("type") == "thermal":
            continue
        src_elem = edge.source in element_ids
        tgt_elem = edge.target in element_ids
        src_tee = edge.source in tee_ids
        tgt_tee = edge.target in tee_ids

        if src_elem and tgt_elem and not src_tee and not tgt_tee:
            # Standard element-to-element: auto-insert MomentumChamberNode
            junction_id = f"__junction__{edge.source}__{edge.target}"
            if junction_id not in nodes_map:
                junction_node = MomentumChamberNode(junction_id)
                net.add_node(junction_node)
                nodes_map[junction_id] = junction_node
            edge_junction_map[(edge.source, edge.target)] = junction_id

        elif src_tee and tgt_elem:
            # Tee port → element (e.g. tee.common-source → channel): auto-insert PlenumNode
            port = _port_from_handle(edge.sourceHandle)
            if port and (edge.source, port) not in tee_port_node_map:
                jid = f"__tee_jct__{edge.source}_{port}"
                if jid not in nodes_map:
                    jn = PlenumNode(jid)
                    net.add_node(jn)
                    nodes_map[jid] = jn
                tee_port_node_map[(edge.source, port)] = jid
                edge_junction_map[(edge.source, edge.target)] = jid

        elif src_elem and tgt_tee:
            # Element → tee port (e.g. channel → tee.straight-target): auto-insert PlenumNode
            port = _port_from_handle(edge.targetHandle)
            if port and (edge.target, port) not in tee_port_node_map:
                jid = f"__tee_jct__{edge.target}_{port}"
                if jid not in nodes_map:
                    jn = PlenumNode(jid)
                    net.add_node(jn)
                    nodes_map[jid] = jn
                tee_port_node_map[(edge.target, port)] = jid
                edge_junction_map[(edge.source, edge.target)] = jid

    # 2. Second Pass: Create Elements
    for elem_node_schema in element_nodes:
        elem_id = elem_node_schema.id
        elem_data = elem_node_schema.data
        elem_type = elem_node_schema.type

        # 3-port element: bypass the 2-port source/target logic entirely.
        if elem_type == "tee_junction":
            _d = TeeJunctionData(**elem_data)
            _tee = TeeJunctionElement(
                id=elem_id,
                common_node=_find_tee_port(
                    schema.edges, elem_id, "common", nodes_map, tee_port_node_map
                ),
                straight_node=_find_tee_port(
                    schema.edges, elem_id, "straight", nodes_map, tee_port_node_map
                ),
                branch_node=_find_tee_port(
                    schema.edges, elem_id, "branch", nodes_map, tee_port_node_map
                ),
                theta=_math.radians(_d.theta_deg),
                F_C=_d.F_C,
                psi=_d.psi,
                tee_type=_d.tee_type,
            )
            _tee.initial_guess = _expand_initial_guess(_d.initial_guess, elem_id)
            net.add_element(_tee)
            continue

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
                friction_model=data.friction_model,
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
        elif elem_type == "area_change":
            data = AreaChangeData(**elem_data)
            elem = AreaChangeElement(
                id=elem_id,
                from_node=source_id,
                to_node=target_id,
                F0=data.F0,
                F1=data.F1,
                model_type=data.model_type,
                length=data.length,
                D_h=data.D_h,
            )
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
