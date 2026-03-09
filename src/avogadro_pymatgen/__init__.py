"""Avogadro plugin for building structures using pymatgen."""

import argparse
import json
import sys


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("feature")
    parser.add_argument("--lang", nargs="?", default="en")
    parser.add_argument("--debug", action="store_true")
    args = parser.parse_args()

    avo_input = json.load(sys.stdin)
    output = None

    match args.feature:
        case "slab":
            from .slab import run
            output = run(avo_input)

    if output is not None:
        print(json.dumps(output))
