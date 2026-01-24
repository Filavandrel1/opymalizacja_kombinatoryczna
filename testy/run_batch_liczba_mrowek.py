#!/usr/bin/env python3
"""Uruchamia tryb batch w programie C++ i wypisuje wyniki.

Wymaga, aby binarka `./main` była zbudowana w katalogu projektu.
Program C++ obsługuje tryb: `./main --batch test.txt`.

Skrypt:
- buduje program (g++ -std=c++17)
- uruchamia batch dla podanego pliku instancji
- wypisuje stdout (tabela TSV)
"""

from __future__ import annotations

import argparse
import subprocess
from pathlib import Path


def run(cmd: list[str], cwd: Path) -> None:
    subprocess.run(cmd, cwd=str(cwd), check=True)


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("instance", help="Plik instancji, np. test.txt")
    args = p.parse_args()

    root = Path(__file__).resolve().parent

    # Build
    run(["g++", "-std=c++17", "-O2", "-g", "main.cpp", "-o", "main"], cwd=root)

    # Execute batch
    proc = subprocess.run(
        ["./main", "--batch", args.instance],
        cwd=str(root),
        check=True,
        capture_output=True,
        text=True,
    )
    print(proc.stdout)
    if proc.stderr:
        print(proc.stderr, end="", file=__import__("sys").stderr)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
