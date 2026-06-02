"""
No_flame_sim/__main__.py
------------------------
Run as:  python -m No_flame_sim
         python -m No_flame_sim --phi 0.7 1.0 1.3
         python -m No_flame_sim --coord x_raw --csv
         python -m No_flame_sim --mechanism gri30.yaml --verbose
"""
import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))


def main():
    parser = argparse.ArgumentParser(
        description="Cantera 1-D H₂/air freely-propagating flame simulations with NO emissions (AE4262 §4.3)"
    )
    parser.add_argument(
        "--phi", nargs="+", type=float, default=None,
        metavar="PHI",
        help="Equivalence ratios to simulate (default: 0.7 1.0 1.3).",
    )
    parser.add_argument(
        "--coord", choices=["x_shifted", "x_raw", "c"], default=None,
        help="Spatial coordinate for plots (default: x_shifted).",
    )
    parser.add_argument(
        "--mechanism", default=None,
        help=(
            "Cantera mechanism YAML file (default: gri30.yaml). "
            "NOTE: NO profiles require GRI-Mech 3.0. "
            "h2o2.yaml is faster but contains no nitrogen chemistry."
        ),
    )
    parser.add_argument(
        "--csv", action="store_true",
        help="Also save per-phi results as CSV files.",
    )
    parser.add_argument(
        "--verbose", action="store_true",
        help="Pass loglevel=1 to Cantera solver.",
    )
    args = parser.parse_args()

    import No_flame_sim as sim

    if args.coord is not None:
        sim.COORD = args.coord
    if args.mechanism is not None:
        sim.MECHANISM = args.mechanism

    print(f"Mechanism : {sim.MECHANISM}")
    print(f"Coordinate: {sim.COORD}")

    sim.run_all_phi(
        phi_values=args.phi,
        loglevel=1 if args.verbose else 0,
        save_csv=args.csv,
    )


if __name__ == "__main__":
    main()