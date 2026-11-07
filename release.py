#!/usr/bin/env python3

###############################################################################
#
#    Copyright (C) 2020 Ben Woodcroft
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
###############################################################################

__author__ = "Ben Woodcroft"
__copyright__ = "Copyright 2023"
__credits__ = ["Ben Woodcroft"]
__license__ = "GPL3"
__maintainer__ = "Ben Woodcroft"
__email__ = "benjwoodcroft near gmail.com"
__status__ = "Development"

import argparse
import os
import subprocess
import sys
from pathlib import Path

# No cargo update here: the tag has to build the Cargo.lock that was measured.
GATES: tuple[tuple[list[str], dict[str, str]], ...] = (
    (["cargo", "fmt", "--check"], {}),
    (
        [
            "cargo",
            "clippy",
            "--locked",
            "--all-targets",
            "--all-features",
            "--",
            "-D",
            "warnings",
        ],
        {},
    ),
    (["cargo", "doc", "--locked", "--no-deps"], {"RUSTDOCFLAGS": "-D warnings"}),
    (["cargo", "test", "--locked"], {}),
    (["cargo", "package", "--locked"], {}),
)


class Args(argparse.Namespace):
    version: str = ""


def clean() -> bool:
    return all(
        subprocess.run(["git", "diff", "--exit-code", *staged], check=False).returncode
        == 0
        for staged in ([], ["--cached"])
    )


def cargo_version() -> str:
    return Path("Cargo.toml").read_text().split('version = "')[1].split('"')[0]


def main() -> None:
    parser = argparse.ArgumentParser()
    _ = parser.add_argument("--version", required=True, help="e.g. 0.0.6")
    args = parser.parse_args(namespace=Args())

    if args.version.startswith("v"):
        sys.exit("--version shouldn't have a leading 'v'")
    if not clean():
        sys.exit(
            "There appears to be uncommited code changes. Please commit everything first."
        )
    if cargo_version() != args.version:
        sys.exit(
            f"Cargo.toml version ({cargo_version()}) doesn't match --version ({args.version})"
        )

    for command, extra in GATES:
        print(f"$ {' '.join(command)}", flush=True)
        _ = subprocess.run(command, check=True, env=os.environ | extra)

    _ = subprocess.run(["git", "tag", f"v{args.version}"], check=True)
    print("Now run cargo publish --locked && git push && git push --tags")


if __name__ == "__main__":
    main()
