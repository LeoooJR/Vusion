#!/usr/bin/python3

from cli import EntryPoint
from sys import exit, stderr

def main():

    try:

        EntryPoint().launch()

    except SystemExit as e:

        print(f"Error: {e}", file=stderr)
        exit(1)


if __name__ == "__main__":

    main()
