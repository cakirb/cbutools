"""Command line entry point: ``python -m cbutools``.

Exists so that a consumer which is not Python -- a Nextflow config, a shell
script -- can read recipe metadata without importing the package.
"""

import sys

from .recipes import _main as _recipes_main

USAGE = """usage: python -m cbutools recipes [NAME] [--names-only] [--include-pending]

Print cbutools barcode recipe metadata as JSON.
"""


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    if argv and argv[0] == "recipes":
        return _recipes_main(argv[1:])
    if argv and argv[0] in ("-h", "--help"):
        print(USAGE)
        return 0
    print(USAGE, file=sys.stderr)
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
