#!/usr/bin/python3

"""
Main script for the program.
"""

from sys import exit

from loguru import logger
from rich import box
from rich.panel import Panel

from cli import EntryPoint
from console import print_stderr


def main() -> None:
    """
    Main function for the program.
    """

    # Try to launch the program
    try:
        # Launch the program
        exit_code: int = EntryPoint().launch()
        logger.success(f"Program exited with status code {exit_code}.")
    except SystemExit as e:
        # Retrieve the original cause (parent error) if present
        cause = e.__cause__
        if cause is not None:
            message = f"{str(e)}\n[bold red]Caused by:[/bold red] {repr(cause)}"
        else:
            message = str(e)
        # Print the caught exception and its cause to standard error stream
        print_stderr(
            Panel.fit(
                message,
                box=box.ROUNDED,
                title="Execution error",
                subtitle="System exit as 1",
                highlight=True,
            ),
        )
        # Exit the program with status code 1 (error) as Unix convention
        exit(1)


# Call the main function if the script is executed directly
if __name__ == "__main__":

    main()
