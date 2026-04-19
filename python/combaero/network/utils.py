from typing import TYPE_CHECKING

if TYPE_CHECKING:
    import pandas as pd

    from .graph import FlowNetwork


def results_to_dataframe(net: "FlowNetwork") -> "pd.DataFrame":
    """
    Extracts the current state of a solved network into a pandas DataFrame.

    Args:
        net: A solved FlowNetwork instance.

    Returns:
        pd.DataFrame containing thermodynamics and flow results for nodes and elements.
    """
    import pandas as pd

    data = []

    # 1. Process Nodes
    for node_name, node in net.nodes.items():
        state = node.get_state()
        row = {
            "type": "node",
            "id": node_name,
            "T": state.T,
            "P": state.P,
            "rho": state.rho,
            "h": state.h,
            "s": state.s,
            "M": state.M,
            "v": state.v,
        }
        # Add species fractions
        for i, y_val in enumerate(state.Y):
            row[f"Y_{i}"] = y_val

        data.append(row)

    # 2. Process Elements (Edges)
    for elem_name, elem in net.elements.items():
        row = {
            "type": "element",
            "id": elem_name,
            "m_dot": elem.m_dot if hasattr(elem, "m_dot") else 0.0,
            "dP": elem.dP if hasattr(elem, "dP") else 0.0,
        }
        data.append(row)

    return pd.DataFrame(data)
