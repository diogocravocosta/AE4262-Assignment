"""
ch4_flame_sim/__main__.py
-------------------------
Run as:  python -m ch4_flame_sim
         python -m ch4_flame_sim --phi 0.8 1.0 1.2
         python -m ch4_flame_sim --coord x_raw --csv
         python -m ch4_flame_sim --verbose
"""
import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Cantera 1-D CH₄/air freely-propagating flame simulation "
            "with RADIS emission spectrum synthesis (GRI-Mech 3.0)."
        )
    )
    parser.add_argument(
        "--phi", nargs="+", type=float, default=None,
        metavar="PHI",
        help="Equivalence ratios to simulate (default: 0.8 1.0 1.2).",
    )
    parser.add_argument(
        "--coord", choices=["x_shifted", "x_raw", "c"], default=None,
        help=(
            "Spatial coordinate for plots (default: x_shifted). "
            "'c' = progress variable (0=reactants → 1=products)."
        ),
    )
    parser.add_argument(
        "--mechanism", default=None,
        help="Cantera mechanism YAML file (default: gri30.yaml).",
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

    import ch4_flame_sim as sim

    print("Mechanism:", sim.MECHANISM)

    if args.coord is not None:
        sim.COORD = args.coord
    if args.mechanism is not None:
        sim.MECHANISM = args.mechanism

    sim.run_all_phi(
        phi_values=args.phi,
        loglevel=1 if args.verbose else 0,
        save_csv=args.csv,
    )


if __name__ == "__main__":
    main()