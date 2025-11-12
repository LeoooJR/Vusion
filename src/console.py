"""
A module to manage the console outputs.
"""

from typing import Final

from rich.console import Console
from rich.style import Style
from rich.theme import Theme

# Theme for the console
THEMES: Final[Theme] = Theme(
    styles={
        "result": Style(color="green1", bold=True),
        "info": Style(color="sky_blue3", bold=True),
        "warning": Style(color="orange_red1", bold=True),
        "error": Style(color="red1", bold=True),
    }
)

# Console for standard output stream
STDOUT_CONSOLE: Final[Console] = Console(color_system="auto", theme=THEMES)

# Console for standard error stream
STDERR_CONSOLE: Final[Console] = Console(color_system="auto", stderr=True, theme=THEMES)


# Print a message to the standard output stream
def print_stdout(message: object) -> None:
    """
    Print a message to the standard output stream.
    """
    STDOUT_CONSOLE.print(message, style="result")


# Print a message to the standard error stream
def print_stderr(message: object) -> None:
    """
    Print a message to the standard error stream.
    """
    STDERR_CONSOLE.print(message, style="error")
