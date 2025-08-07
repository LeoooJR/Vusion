from rich.console import Console
from rich.style import Style
from rich.theme import Theme

# Theme for the console
theme = Theme(styles={
    "result": Style(color="green1", bold=True),
    "info": Style(color="sky_blue3", bold=True),
    "warning": Style(color="orange_red1", bold=True),
    "error": Style(color="red1", bold=True),
})

# Console for standard output stream
stdout_console = Console(color_system="auto", theme=theme)

# Console for standard error stream
stderr_console = Console(color_system="auto", stderr=True, theme=theme)