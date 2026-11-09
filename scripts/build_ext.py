"""Build libalp.so with CMake and copy it into the package.

poetry-core runs this file as a script during ``pip install`` / ``poetry build``.
It can also be run by hand: ``python scripts/build_ext.py``.
"""
import multiprocessing
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
BUILD_DIR = ROOT / "build"
LIB = "libalp.so"


def main() -> None:
    cmake = shutil.which("cmake")
    if cmake is None:
        sys.exit("ERROR: cmake not found on PATH (pip install cmake).")
    subprocess.run([cmake, "-B", str(BUILD_DIR), "-S", str(ROOT),
                    "-DCMAKE_BUILD_TYPE=Release"], check=True)
    subprocess.run([cmake, "--build", str(BUILD_DIR), "-j",
                    str(multiprocessing.cpu_count())], check=True)
    built = BUILD_DIR / LIB
    if not built.exists():
        sys.exit(f"ERROR: expected {built} after build.")
    shutil.copy2(built, ROOT / "alp" / LIB)


main()
