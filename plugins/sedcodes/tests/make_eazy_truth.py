"""This script does not call an external package and does not write a substitute truth catalog."""
import sys


def main():
    sys.stderr.write(
        "the external package is not called; this script does not generate "
        "a substitute truth catalog.\n"
    )
    return 2


if __name__ == "__main__":
    sys.exit(main())
