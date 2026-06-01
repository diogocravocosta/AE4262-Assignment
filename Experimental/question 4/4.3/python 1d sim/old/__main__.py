"""
h2_flame_sim/__main__.py
------------------------
Run as:  python -m h2_flame_sim
         python -m h2_flame_sim --phi 0.7 1.0 1.3   # specific phi only
         python -m h2_flame_sim --coord x_raw --csv  # change coordinate, save CSV
"""
import argparse
import sys

def main():
    parser = argparse.ArgumentParser(
        description="Cantera 1-D H₂/air freely-propagating flame simulations (AE4262 §4.3)"
    )
    parser.add_argument(
        "--phi", nargs="+", type=float, default=None,
        metavar="PHI",
        help="Equivalence ratios to simulate (default: 0.5 0.6 … 1.3)."
    )
    parser.add_argument(
        "--coord", choices=["x_shifted", "x_raw", "Z"], default=None,
        help="Spatial coordinate for plots (default: x_shifted). "
             "NOTE: 'Z' is not meaningful for premixed flames – see __init__.py docstring."
    )
    parser.add_argument(
        "--mechanism", default=None,
        help="Cantera mechanism YAML file (default: h2o2.yaml)."
    )
    parser.add_argument(
        "--csv", action="store_true",
        help="Also save per-phi results as CSV files."
    )
    parser.add_argument(
        "--verbose", action="store_true",
        help="Pass loglevel=1 to Cantera solver."
    )
    args = parser.parse_args()

    # Apply overrides to module-level constants before importing run logic
    import h2_flame_sim as sim

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
