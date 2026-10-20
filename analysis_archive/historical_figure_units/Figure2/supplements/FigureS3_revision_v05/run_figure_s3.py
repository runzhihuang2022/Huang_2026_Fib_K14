import argparse

from scripts.analyze_figure_s3_v04 import main as analyze
from scripts.build_figure_s3_v05 import build


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--reanalyze",
        action="store_true",
        help="Recompute source tables from the 11 registered H5AD sections.",
    )
    args = parser.parse_args()
    if args.reanalyze:
        analyze()
    build()
