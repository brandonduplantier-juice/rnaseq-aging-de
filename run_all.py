"""
Run the full RNA-seq aging DE pipeline: download, differential expression, plots.

Skips the download if data/counts.csv already exists, so reruns are fast.
Usage:  python run_all.py   (or --force-download to refetch)
"""

import os
import runpy
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "src")
DATA = os.path.join(HERE, "data", "counts.csv")


def run(stage):
    print("\n=== {} ===".format(stage))
    runpy.run_path(os.path.join(SRC, stage), run_name="__main__")


def main():
    if os.path.exists(DATA) and "--force-download" not in sys.argv:
        print("[run_all] data/counts.csv exists, skipping download "
              "(pass --force-download to refetch)")
    else:
        run("download_data.py")
    run("de_analysis.py")
    run("plots.py")
    print("\n[run_all] done. See results/ for the table, metrics, and plots.")


if __name__ == "__main__":
    main()
