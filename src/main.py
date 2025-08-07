#!/usr/bin/python3

from cli import EntryPoint
from console import stderr_console
from rich.panel import Panel
from rich import box
from sys import exit

def main():

    # Try to launch the program
    try:
        # Launch the program
        EntryPoint().launch()
        # Exit the program with status code 0 as Unix convention
        exit(0)

    except SystemExit as e:
        # Print the catched exception to standard error stream
        stderr_console.print(Panel.fit(str(e), box=box.ROUNDED, title="Execution error", subtitle="System exit as 1", highlight=True), style="error")
        # Exit the program with status code 1 as Unix convention
        exit(1)

# Call the main function if the script is executed directly
if __name__ == "__main__":

    main()
