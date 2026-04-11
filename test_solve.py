import traceback
import sys
from gui.backend.graph_builder import build_network_from_schema
from combaero.network.solver import NetworkSolver
from gui.backend.schemas import NetworkGraphSchema
import json

try:
    with open('gui/tmp/fuel_air_transport_comb_test.json', 'r') as f:
        data = json.load(f)
    graph_data = NetworkGraphSchema(**data)
    # Let's import main._solve_sync to test the exact pipeline!
    from gui.backend.main import _solve_sync
    result, node_results, element_results, edge_results, _ = _solve_sync(graph_data)

    print("Diagnostics pipeline completed!")

    plenum = node_results.get('node_1775378197276')
    if plenum:
        print("Plenum Y =", plenum.state.Y)

    combustor = node_results.get('node_1775378191205')
    if combustor:
        print("Combustor T =", combustor.state.T)
        print("Combustor phi =", combustor.state.model_extra.get('phi'))


except Exception:
    traceback.print_exc(file=sys.stdout)
