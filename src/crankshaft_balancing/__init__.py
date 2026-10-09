"""crankshaft-balancing: моменты от центробежных сил и подбор противовесов коленвала."""
from .analysis import BalanceResult, ShaftForces, analyze, compute_forces, raw_moments
from .config import EngineConfig
from .constraints import CounterweightConstraints, resolve_constraints
from .forces import (centrifugal_force_magnitude, centrifugal_force_vectors,
                     distance_from_right_end, radial_distance)
from .geometry import (CrankshaftGeometry, ShaftComponent, build_geometry,
                       derive_crankshaft_dimensions)
from .loader import (load_constraints, load_counterweight_positions,
                     load_shaft)
from .optimization import (CounterweightLP, build_counterweight_lp,
                           make_initial_coefficients, optimize_counterweights,
                           optimize_counterweights_history)
from .report import (CounterweightRow, counterweight_rows,
                     format_counterweight_table, rows_from_components,
                     save_counterweight_csv)
from .sections import (NMM_TO_NM, heaviside_ramp, section_contribution,
                       section_grid, total_distribution)
from .verification import (sweep_sensitivity, verify_kkt_conditions,
                           verify_optimum_via_duality)

__version__ = "0.1.0"

__all__ = [
    "BalanceResult", "ShaftForces", "analyze", "compute_forces", "raw_moments",
    "EngineConfig", "CrankshaftGeometry", "ShaftComponent", "build_geometry",
    "derive_crankshaft_dimensions", "load_shaft", "load_constraints",
    "load_counterweight_positions",
    "CounterweightConstraints", "resolve_constraints",
    "centrifugal_force_magnitude", "centrifugal_force_vectors",
    "distance_from_right_end", "radial_distance",
    "CounterweightLP", "build_counterweight_lp", "make_initial_coefficients",
    "optimize_counterweights", "optimize_counterweights_history",
    "CounterweightRow", "counterweight_rows", "format_counterweight_table",
    "rows_from_components", "save_counterweight_csv",
    "NMM_TO_NM", "heaviside_ramp", "section_contribution", "section_grid",
    "total_distribution",
    "sweep_sensitivity", "verify_kkt_conditions", "verify_optimum_via_duality",
]
