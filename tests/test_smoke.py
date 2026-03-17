from research.control_gap_junction.model import SimulationConfig
from research.control_gap_junction.simulator import run_simulation


def test_control_package_smoke() -> None:
    config = SimulationConfig(lattice_size=4, steps=2)
    history = run_simulation(config)
    assert len(history) == 2
