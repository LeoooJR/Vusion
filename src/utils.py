"""
Utility functions for the program.
"""

import ast
import os
import re
from datetime import datetime, timezone
from hashlib import sha256
from pathlib import Path
from typing import Callable, Final, final


class Cache:
    """
    Least Recently Used (LRU) cache implementation.
    """

    def __init__(self, func: Callable, max_size: Final[int | None] = None) -> None:
        """
        Initialize the cache.

        Args:
            func: The function to cache output from.
            max_size: The maximum size of the cache.
        """

        # Maximum size of the cache
        self.max_size: Final[int] = max_size

        # Dictionnary to implement the cache
        # The key is the hash of the arguments
        # The value is the result of the function
        self.cache: Final[dict] = {}

        # The function to cache output from.
        self.func: Final[Callable] = func

    def add(self, args: list[str], key: object):
        """
        Add the result of the function to the cache

        Args:
            args: The arguments of the function
            key: The key of the cache
        """

        # If the cache is full, clear it
        if self.max_size and len(self.cache) >= self.max_size:

            self.cache.clear()

        # Add the result of the function to the cache
        self.cache[key] = self.func(*args)

    def call(self, args: list[str], key: object):
        """
        Call the function and add the result to the cache

        Args:
            args: The arguments of the function
            key: The key of the cache
        """
        # Check if the key is hashable
        if key.__hash__:

            # If the key is not in the cache,
            # add the result of the function to the cache
            if not key in self.cache:

                self.add(args, key)

            # O(1) access to the cache
            # Return the result of the function
            # from the cache
            return self.cache[key]

        raise TypeError(f"Key is not hashable: {key}")


class PluginPythonChecker(ast.NodeVisitor):
    """
    Check if a Python file is safe.
    """

    # List of dangerous calls.
    DANGEROUS_CALL: list[str] = ["exec", "eval", "compile"]

    def __init__(self):

        # List of imports in the Python file.
        self.imports: list[str] = []

        # List of calls in the Python file.
        self.calls: list[str] = []

        # List of dangerous calls in the Python file.
        self.not_safe_calls: list[str] = []

    @final
    def visit_Import(self, node: ast.Import) -> None:
        """
        Visit an Import node.
        """

        for alias in node.names:

            self.imports.append(alias.name)

        self.generic_visit(node)

    @final
    def visit_ImportFrom(self, node: ast.ImportFrom) -> None:
        """
        Visit an ImportFrom node.
        """

        if node.module:

            self.imports.append(node.module)

        self.generic_visit(node)

    @final
    def visit_Call(self, node: ast.Call) -> None:
        """
        Visit a Call node.
        """

        if isinstance(node.func, ast.Name):

            if node.func.id in self.DANGEROUS_CALL:

                self.not_safe_calls.append(node.func.id)

        self.generic_visit(node)


# ===========================================================================================
# Basics functions on dictionary
# ===========================================================================================


def merge_collections(collections: list[object]) -> object:
    """
    Merge a list of collections.

    Args:
        collections: A list of collections to merge.

    Returns:
        A merged collection of the same type as inputed collections.
    """
    if isinstance(
        collections[0], dict
    ):  # Check if the first collection is a dictionary

        # Check if all the collections are dictionaries
        is_dict: list[bool] = list(
            map(lambda collection: isinstance(collection, dict), collections)
        )

        if all(is_dict):

            output: dict = {}

            for collection in collections:

                output.update(collection)

        else:

            raise ValueError(
                f"Not all of the collections being merged are of the same data type: {is_dict}"
            )

    else:

        raise ValueError(f"{type(collections[0])} cannot be merged.")

    return output


# ===========================================================================================
# Filsystem
# ===========================================================================================


def get_project_dir() -> str:
    """
    Get the project directory.
    """
    return os.path.dirname(os.path.abspath(__file__))


def get_or_create_config_dir() -> Path:
    """Get user configuration directory following XDG spec"""
    if os.name == "nt":  # Windows
        base = os.environ.get("APPDATA", Path.home() / "AppData" / "Roaming")
    else:  # Unix-like
        base = os.environ.get("XDG_CONFIG_HOME", Path.home() / ".config")

    # The path tho the callers configuration directory
    config: Path = Path(base) / "vusion" / "callers"

    # If the directory does not exist, create it
    if not config.exists():
        create_config_dir(config)

    # Return the path to the callers configuration directory
    return config


def create_config_dir(path: Path) -> Path:
    """
    Create a filesystem .config directory.

    Args:
        path (Path): The path to the directory to create.

    Returns:
        Path: The path to the created directory.
    """

    path.mkdir(parents=True, exist_ok=True)

    # Create an empty __init__.py file
    path.joinpath("__init__.py").touch()

    # Return the path to the callers configuration directory
    return path


def hash_file(path: Path) -> str:
    """
    Hash a file.

    Args:
        path (Path): The path to the file to hash.

    Returns:
        str: The hash of the file.
    """
    # Check if the file exists
    if path.exists():

        # Create a SHA-256 hash object
        hash_object_file: sha256 = sha256()

        # Open the file and read its content
        with open(path, mode="r", encoding="utf-8") as plugin:

            # Read the content of the file
            for line in plugin:

                # Update the hash with the content of the file
                hash_object_file.update(line.encode())

        # Return the hash of the file
        return hash_object_file.hexdigest()

    # Raise an error if the file does not exist
    raise FileNotFoundError(f"File {path} not found.")


def file_infos(path: str) -> dict:
    """
    Get file stats.

    Args:
        path (str): The path to the file to get stats from.

    Returns:
        A dictionary containing the file stats.
    """
    statinfo = os.stat(path)

    return {
        "basename": os.path.basename(path),
        "path": os.path.dirname(path),
        "size": round(statinfo.st_size / pow(1024, 2), 2),
        "mtime": datetime.fromtimestamp(statinfo.st_mtime, tz=timezone.utc),
    }


def clean(files: list[Path]):
    """
    Remove a list of files from filesystem.

    Args:
        files: A list of files to remove.
    """

    # Remove each file from the filesystem
    for file in files:

        file.unlink(missing_ok=True)


# ===========================================================================================
# Functions on variants
# ===========================================================================================


def estimate_brc_r_e(variant, pileup_line_info):
    """
    estimate <BRC/R/E> (background read counts/ratio/enrichment)
    BRC : background read counts
    BRR : background read ratio
    BRE : background read enrichment

    Parameters:
    - dic (dict): A dictionary containing variant information,
                  including the ALT Read Count Ratio (ARR).
    - variant_key (str): The key representing the specific variant in the dictionary.
    - pileup_line_info (list): List containing information from pileup

    Returns : a tuple containing the estimated BRC, BRR, and BRE.
    """
    ref_and_alt_read_count = (
        variant["sample"]["ARC-"]
        + variant["sample"]["ARC+"]
        + variant["sample"]["RRC-"]
        + variant["sample"]["RRC+"]
    )
    ins_read_counts = 0

    total_read_count = variant["sample"]["TRC"]

    if pileup_line_info[14] != "None":
        # Get number of read with ins format is A:1,0
        tmp_table_count = re.findall(r"\d+", pileup_line_info[14])
        for tmp_count in tmp_table_count:
            ins_read_counts += int(tmp_count)

    if variant["type"] == "INS":
        ins_read_counts -= variant["sample"]["ARC+"] + variant["sample"]["ARC-"]

    # Removing DEL counts if variant is at some position of a DEL.
    # Because we miss valid variants in specific case like that
    # Clintool bug where non-existing deletion is reported and start with A, C, T or G
    if (variant["type"] != "DEL") and (pileup_line_info[15] != "None"):

        # Escaping clintool bug where non-existing deletion is reported and start with A, C, T or G
        if pileup_line_info[15][0] == "*":

            del_read_count = int(
                pileup_line_info[15].strip().split(";")[0].split(":")[1]
            )
            tmp_total_read_count = total_read_count - int(del_read_count)

            variant["sample"]["BRC"] = (
                tmp_total_read_count
                + ins_read_counts
                - min([tmp_total_read_count, ref_and_alt_read_count])
            )

        else:

            variant["sample"]["BRC"] = (
                total_read_count
                + ins_read_counts
                - min([total_read_count, ref_and_alt_read_count])
            )

    elif variant["type"] != "DEL":

        variant["sample"]["BRC"] = (
            total_read_count
            + ins_read_counts
            - min([total_read_count, ref_and_alt_read_count])
        )
    else:
        # 230413 Not counting snp background for Deletion
        del_read_count = pileup_line_info[15].strip().split(";")
        del_alt_count = 0
        for del_info in del_read_count:
            tmp_del_info = del_info.split(":")
            if (tmp_del_info[0] != "*") and (
                tmp_del_info[0] != (variant["collection"]["REF"])[1:]
            ):
                del_alt_count += int(tmp_del_info[1].split(",")[0]) + int(
                    tmp_del_info[1].split(",")[1]
                )
        variant["sample"]["BRC"] = del_alt_count

    background_read_counts = variant["sample"]["BRC"]

    variant["sample"]["BRR"] = round(background_read_counts / total_read_count, 5)
    background_read_ratio = variant["sample"]["BRR"]

    alt_read_count_ratio = variant["sample"]["ARR"] / 100
    variant["sample"]["BRE"] = background_read_ratio / (
        alt_read_count_ratio + background_read_ratio
    )
    # scale ratio from [0-1] to [0-100]
    variant["sample"]["BRE"] = round(float(variant["sample"]["BRE"]) * 100, 5)

    variant["sample"]["BRR"] = round(float(variant["sample"]["BRR"]) * 100, 5)

    return (
        variant["sample"]["BRC"],
        variant["sample"]["BRR"],
        variant["sample"]["BRE"],
    )


def categorize_background_signal(variant, thresholds):
    """
    categorize background based on background read enrichment (BRE) thresholds :

    Parameters:
    - dic (dict): A dictionary containing variant information,
                  including background read enrichment (BRE).
    - variant_key (str): The key representing the specific variant in the dictionary.
    - thresholds (list): A list of thresholds passed in mandatory options

    Returns:
    str: Returns the categorized background signal,
         which could be one of the following: 'LNO', 'PNO', 'PCL', or 'LCL'.

    Note:
    - The input dictionary is expected to have 'sample' containing 'ARR' and 'BRE'.
    - update : 230413 If variant ratio is higher that 30% we consider it automatically clean

    L/P-NO Likely/Probably NOisy
    L/P-CL Likely/Probably CLean
    [ LCL [ PCL [ PNO [          LNO          ]
    |     |     |     |                       |
    0    [6]   [7]   [8]                    100 (BRE)
    """
    alt_read_count_ratio = float(variant["sample"]["ARR"])
    if alt_read_count_ratio >= 30:
        variant["sample"]["BKG"] = "PCL"
    else:
        background_read_enrichment = float(variant["sample"]["BRE"])

        if background_read_enrichment >= thresholds[8]:
            variant["sample"]["BKG"] = "LNO"
        elif background_read_enrichment >= thresholds[7]:
            variant["sample"]["BKG"] = "PNO"
        elif background_read_enrichment >= thresholds[6]:
            variant["sample"]["BKG"] = "PCL"
        else:
            variant["sample"]["BKG"] = "LCL"

    return variant["sample"]["BKG"]
