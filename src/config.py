"""
A module to parse the config file.
"""

import re
from copy import copy
from dataclasses import dataclass
from typing import List, Optional, Union

import cerberus
import yaml
from lark import Lark, Transformer, UnexpectedInput
from loguru import logger
from rich.panel import Panel
from rich.text import Text

from console import print_stderr
from exceptions import ConfigError


@dataclass
class Metadatas:
    """
    A class to store metadata.
    """

    header: Optional[str] = None
    index: Optional[int] = None
    unit: Optional[str] = None


@dataclass
class Term:
    """
    A class to store a term.
    """

    field: str
    metadata: Optional[Metadatas] = None


@dataclass
class Expression:
    """
    A class to store an expression.
    """

    terms: List[Union[Term, "Expression"]]
    operator: Optional[str] = None


class TreeToExpression(Transformer):
    """
    A class to transform a tree to an expression.
    """

    def STRING(self, token):
        """
        Convert a string token to a string.
        """
        return str(token.value)

    def INDEX(self, token):
        """
        Convert an integer token to an integer.
        """
        return int(token.value)

    def HEADER(self, token):
        """
        Convert a header token to a string.
        """
        return str(token.value)

    def UNIT(self, token):
        """
        Convert a unit token to a string.
        """
        return str(token.value)

    def OPERATOR(self, token):
        """
        Convert an operator token to a string.
        """
        return str(token.value)

    def FIELD(self, token):
        """
        Convert a field token to a string.
        """
        return str(token.value)

    def metadata(self, item):
        """
        Convert a metadata item to a metadata object.
        """
        return item[0]

    def indexing(self, items: list):
        """
        Convert a list of items to a metadata object.
        """
        header = None
        index = None
        unit = None

        for item in items:
            if isinstance(item, int):
                index = item
            elif isinstance(item, str) and item == "%":
                unit = item
            else:
                header = item

        return Metadatas(header=header, index=index, unit=unit)

    def term(self, items):
        """
        Convert a list of items to a term object.
        """
        if len(items) == 1:
            return Term(field=items[0])
        return Term(field=items[0], metadata=items[1])

    def expression(self, items):
        """
        Convert a list of items to an expression object.
        """
        if len(items) == 1:
            return items[0]

        terms = []
        current_operator = None

        for item in items:
            if isinstance(item, str) and item in ["+", "-", "*", "/"]:
                current_operator = item
            else:
                if current_operator and terms:
                    terms.append(
                        Expression(terms=[terms[-1], item], operator=current_operator)
                    )
                    current_operator = None
                else:
                    terms.append(item)

        return terms[-1] if len(terms) > 1 else terms[0]


class ExpressionVisitor:
    """Class to handle visiting of expression."""

    def visit_EXPRESSION(self, expression) -> dict:
        """
        Visit an expression object and return a dictionary.
        """
        if isinstance(expression, Term):
            if isinstance(expression.field, Expression):
                return {
                    "type": "expression",
                    "operator": expression.field.operator,
                    "terms": [
                        self.visit_EXPRESSION(term) for term in expression.field.terms
                    ],
                }
            else:
                return {
                    "type": "term",
                    "field": expression.field,
                    "metadata": (
                        {
                            "header": (
                                expression.metadata.header
                                if expression.metadata
                                else None
                            ),
                            "index": (
                                expression.metadata.index
                                if expression.metadata
                                else None
                            ),
                            "unit": (
                                expression.metadata.unit
                                if expression.metadata
                                else None
                            ),
                        }
                        if expression.metadata
                        else None
                    ),
                }

        return {
            "type": "expression",
            "operator": expression.operator,
            "terms": [self.visit_EXPRESSION(term) for term in expression.terms],
        }


class ExpressionTemplate:
    """Class to handle templating of expressions into different formats."""

    def __init__(self, format_fields: list = None, infos_fields: list = None):
        """Initialize with the list of valid format fields.

        Args:
            format_fields (list): List of valid format field names
        """
        self._format_fields = format_fields

        self._infos_fields = infos_fields

    @property
    def format(self):
        """
        Get the format fields.
        """
        return getattr(self, "_format_fields", None)

    @format.setter
    def format(self, value):
        """
        Set the format fields.
        """
        self._format_fields = value

    @property
    def infos(self):
        """
        Get the infos fields.
        """
        return getattr(self, "_infos_field", None)

    @infos.setter
    def infos(self, value):
        """
        Set the infos fields.
        """
        self._infos_fields = value

    def _is_valid_field(self, field: list) -> bool:
        """
        Check if a field is valid.
        """
        is_in_format: bool = False
        is_in_info: bool = False

        if self._format_fields:
            is_in_format = field in self._format_fields

        if self._infos_fields:
            is_in_info = field in self._infos_fields

        return is_in_format or is_in_info

    def _is_valid_metadata(
        self, field, index: int = None, header: str = None, unit: str = None
    ) -> bool:
        """
        Check if a metadata is valid.
        """
        is_valid_index: bool = False if index else True
        is_valid_header: bool = False if header else True
        is_valid_unit: bool = False if unit else True

        if index:

            is_valid_index = True

        if header:

            if header == "format":
                is_valid_header = (
                    field in self._format_fields if self._format_fields else False
                )
            else:
                is_valid_header = (
                    field in self._infos_fields if self._infos_fields else False
                )

        if unit:

            is_valid_unit = True

        return is_valid_index and is_valid_header and is_valid_unit

    def _format_term(self, term: dict) -> str:
        """Format a single term into a template string.

        Args:
            term (dict): Term dictionary from ExpressionVisitor

        Returns:
            str: Formatted term string
        """
        if self._format_fields or self._infos_fields:

            if not self._is_valid_field(term["field"]):
                raise UnexpectedInput(
                    f"Key {term['field']} not in FORMAT or INFOS fields"
                )

        result = term["field"]

        # Add metadata if present
        if term["metadata"]:
            metadata = []
            if term["metadata"]["index"] is not None:
                if not self._is_valid_metadata(
                    field=term["field"], index=term["metadata"]["index"]
                ):
                    raise UnexpectedInput(
                        f"Index metadata in term {term} is not valid."
                    )
                metadata.append(str(term["metadata"]["index"]))
            if term["metadata"]["header"]:
                if not self._is_valid_metadata(
                    field=term["field"], header=term["metadata"]["header"]
                ):
                    raise UnexpectedInput(
                        f"{term['field']} is not in {term['metadata']['header'].upper()}."
                    )
                metadata.append(term["metadata"]["header"])
            if term["metadata"]["unit"]:
                if not self._is_valid_metadata(
                    field=term["field"], unit=term["metadata"]["unit"]
                ):
                    raise UnexpectedInput(f"Unit metadata in term {term} is not valid.")
                metadata.append(term["metadata"]["unit"])

            if metadata:
                result += f"[{','.join(metadata)}]"

        return result

    def _format_expression(self, expr: dict) -> str:
        """
        Format an expression into a template string.

        Args:
            expr (dict): Expression dictionary from ExpressionVisitor

        Returns:
            str: Formatted expression string
        """
        if expr["type"] == "term":
            return self._format_term(expr)

        # Format each term and join with operator
        terms = [self._format_expression(term) for term in expr["terms"]]
        return f"({terms[0]} {expr['operator']} {terms[1]})"

    def to_template(self, expression: dict) -> str:
        """
        Convert an expression structure to a template string.

        Args:
            expression (dict): Expression dictionary from ExpressionVisitor

        Returns:
            str: Template string
        """
        return self._format_expression(expression)

    def to_python(self, expression: dict) -> str:
        """
        Convert an expression structure to Python code.

        Args:
            expression (dict): Expression dictionary from ExpressionVisitor

        Returns:
            str: Python code string
        """
        if expression["type"] == "term":
            term = self._format_term(expression)
            # Convert VCF format notation to Python dictionary access
            return term.replace("[", "[").replace("]", "]")

        terms = [self.to_python(term) for term in expression["terms"]]
        return f"({terms[0]} {expression['operator']} {terms[1]})"


class ConfigParser:
    """
    A config file (YAML) parser.
    """

    # Schema for the config file (YAML)
    SCHEMA = {
        "caller": {  # Describe a caller
            "type": "dict",
            "schema": {
                "name": {  # Name of the caller
                    "type": "string",
                    "empty": False,
                    "required": True,
                    "forbidden": [  # Built-in variant callers
                        "BCFTools",
                        "Varscan",
                        "Vardict",
                        "Pindel",
                        "Haplotypecaller",
                        "Filt3r",
                        "DeepVariant",
                    ],
                },
                "info": {  # Info fields of the VCF file produced by the caller
                    "type": "string",
                    "empty": False,
                    "required": False,
                    "nullable": True,
                    "default": None,
                    "coerce": lambda x: x.replace(";", ","),
                },
                "format": {  # Format fields of the VCF file produced by the caller
                    "type": "string",
                    "empty": False,
                    "required": True,
                    # "regex": "^[A-Z]{1,}(:[A-Z]{1,})*$",
                    "coerce": lambda x: x.replace(":", ","),
                },
                "genotype": {  # Genotype field
                    "type": "dict",
                    "schema": {
                        "extract": {
                            "type": "string",
                            "empty": False,
                            "nullable": False,
                            "required": True,
                        }  # How to extract the genotype from the VCF file
                    },
                    "required": True,
                    "nullable": False,
                    "empty": False,
                    "dependencies": "format",
                },
                "vaf": {  # Variant allele frequency field
                    "type": "dict",
                    "schema": {
                        "extract": {
                            "type": "string",
                            "empty": False,
                            "nullable": False,
                            "required": True,
                        }  # How to extract the variant allele frequency from the VCF file
                    },
                    "required": True,
                    "empty": False,
                    "dependencies": "format",
                },
                "depth": {  # Depth field
                    "type": "dict",
                    "schema": {
                        "extract": {
                            "type": "string",
                            "empty": False,
                            "nullable": False,
                            "required": True,
                        }  # How to extract the depth from the VCF file
                    },
                    "required": True,
                    "empty": False,
                    "dependencies": "format",
                },
                "rrc": {  # Reference read count field
                    "type": "dict",
                    "schema": {
                        "forward": {
                            "type": "dict",
                            "schema": {
                                "extract": {  # How to extract the forward allele counts from the VCF file
                                    "type": "string",
                                    "required": True,
                                    "empty": False,
                                    "nullable": True,
                                    "default": None,
                                }
                            },
                            "empty": False,
                            "nullable": False,
                            "required": False,
                        },
                        "reverse": {
                            "type": "dict",
                            "schema": {
                                "extract": {  # How to extract the reverse allele counts from the VCF file
                                    "type": "string",
                                    "required": True,
                                    "empty": False,
                                    "nullable": True,
                                    "default": None,
                                }
                            },
                            "empty": False,
                            "nullable": False,
                            "required": False,
                        },
                        "total": {
                            "type": "dict",
                            "schema": {
                                "extract": {  # How to extract the total allele counts from the VCF file
                                    "type": "string",
                                    "required": True,
                                    "empty": False,
                                    "nullable": True,
                                    "default": None,
                                }
                            },
                            "empty": False,
                            "nullable": False,
                            "required": False,
                        },
                    },
                    "required": True,
                    "empty": False,
                    "nullable": True,
                    "dependencies": "format",
                },
                "arc": {  # Alternate read count field
                    "type": "dict",
                    "schema": {
                        "forward": {
                            "type": "dict",
                            "schema": {
                                "extract": {  # How to extract the forward allele counts from the VCF file
                                    "type": "string",
                                    "required": True,
                                    "empty": False,
                                    "nullable": True,
                                    "default": None,
                                }
                            },
                            "empty": False,
                            "nullable": False,
                            "required": False,
                        },
                        "reverse": {
                            "type": "dict",
                            "schema": {
                                "extract": {  # How to extract the reverse allele counts from the VCF file
                                    "type": "string",
                                    "required": True,
                                    "empty": False,
                                    "nullable": True,
                                    "default": None,
                                }
                            },
                            "empty": False,
                            "nullable": False,
                            "required": False,
                        },
                        "total": {
                            "type": "dict",
                            "schema": {
                                "extract": {  # How to extract the total allele counts from the VCF file
                                    "type": "string",
                                    "required": True,
                                    "empty": False,
                                    "nullable": True,
                                    "default": None,
                                }
                            },
                            "empty": False,
                            "nullable": False,
                            "required": False,
                        },
                    },
                    "required": True,
                    "empty": False,
                    "nullable": False,
                    "dependencies": "format",
                },
            },
        }
    }

    # Grammar for the config file (YAML)
    GRAMMAR = r"""
            expression: term (OPERATOR term)*
            term: FIELD indexing? | "(" expression ")"
            indexing: "[" metadata ("," metadata)* "]"
            metadata: INDEX | HEADER | UNIT
            INDEX: INT
            HEADER: "format" | "info"
            UNIT: "%"
            FIELD: STRING            
            OPERATOR: "+"
                        | "-"
                        | "/"
                        | "*"

            %import common.WORD -> STRING
            %import common.INT
            %import common.WS_INLINE

            %ignore WS_INLINE
            """

    # Parser for the DSL (Domain Specific Language)
    DSL_PARSER = Lark(GRAMMAR, start="expression", parser="lalr")

    def __init__(self, path: str):
        """
        Initialize the ConfigParser.

        Args:
            path (str): Path to the config file (YAML).
        """

        self.path: str = path

        self.validator = cerberus.Validator(
            self.SCHEMA
        )  # Validator for the config file schema.

    def valid_schema(self, document) -> bool:
        """
        Validate the config file schema.

        Args:
            document (dict): Config file document.

        Returns:
            bool: True if the config file schema is valid, False otherwise.
        """
        return self.validator.validate(document)

    def parse(self, document) -> dict:
        """
        Parse the config file and return the parsed parameters.
        This method is used to parse the config file after it has been validated.

        Args:
            document (dict): Config file document.

        Returns:
            dict: Parsed parameters.
        """

        parameters = copy(
            document
        )  # Copy the document to avoid modifying the original.

        transformer = TreeToExpression()  # Transformer for the DSL.

        visitor = ExpressionVisitor()  # Visitor for the DSL.

        formatter = ExpressionTemplate()  # Formatter for the DSL.

        # Iterate over the config file fields
        for field in parameters["caller"]:
            # Each following conditional statement is to validate the value of the current field.
            if field == "name":

                # Check that the name only contains alphanumeric letters (a-z) and (0-9), or underscores (_). A valid identifier cannot start with a number, or contain any spaces.
                if not parameters["caller"][field].isidentifier():

                    raise UnexpectedInput(
                        "Name value is not a valid identifier. It must only contain alphanumeric letters (a-z) and (0-9), or underscores (_). A valid identifier cannot start with a number, or contain any spaces."
                    )

            elif field == "info":

                # Check that the info value is a valid VCF info field.
                # A valid VCF info field is a string that only contains uppercase letters (A-Z) separated by commas.
                if not re.match(
                    r"^[A-Z]{1,}(,[A-Z]{1,})*$", parameters["caller"][field]
                ):

                    raise UnexpectedInput(
                        "Info value is not consistent with requested format. It must only contain uppercase letters (A-Z) separated by semicolons."
                    )

                formatter.infos = parameters["caller"][field].split(",")

            elif field == "format":

                if not re.match(
                    r"^[A-Z]{1,}(,[A-Z]{1,})*$", parameters["caller"][field]
                ):

                    raise UnexpectedInput(
                        "Format value is not consistent with VCF format. It must only contain uppercase letters (A-Z) separated by colons."
                    )

                formatter.format = parameters["caller"][field].split(",")

            # For the genotype, depth, and vaf fields, we need to generate the abstract syntax tree (AST) from the DSL.
            # These fields are required, and must be present in the config file.
            elif field in ["genotype", "depth", "vaf"]:

                # Generate the abstract syntax tree (AST) from the DSL.
                ast = self.DSL_PARSER.parse(parameters["caller"][field]["extract"])

                # Transform the AST to an linear expression.
                parameters["caller"][field]["extract"] = transformer.transform(ast)

                # Format the linear expression to a template string.
                formatter.to_template(
                    visitor.visit_EXPRESSION(parameters["caller"][field]["extract"])
                )
            # The following conditional statement handle optional fields (RRC and ARC).
            else:
                # Is the first key of the field present ?
                if parameters["caller"][field]:

                    # Iterate over the possible subfields of the field.
                    for subfield in ["forward", "reverse", "total"]:

                        # Is the subfield present ?
                        if subfield in parameters["caller"][field]:

                            # Is the extract key present ?
                            if parameters["caller"][field][subfield]["extract"]:

                                # Generate the abstract syntax tree (AST) from the DSL.
                                ast = self.DSL_PARSER.parse(
                                    parameters["caller"][field][subfield]["extract"]
                                )

                                # Transform the AST to an linear expression.
                                parameters["caller"][field][subfield]["extract"] = (
                                    transformer.transform(ast)
                                )

                                # Format the linear expression to a template string.
                                formatter.to_template(
                                    visitor.visit_EXPRESSION(
                                        parameters["caller"][field][subfield]["extract"]
                                    )
                                )

                        else:

                            # If the subfield is not present, set the extract key to None.
                            parameters["caller"][field][subfield] = {"extract": None}

                else:

                    # If the field is not present, set the forward, reverse, and total subfields to None.
                    parameters["caller"][field] = {
                        "forward": {"extract": None},
                        "reverse": {"extract": None},
                        "total": {"extract": None},
                    }

        return parameters

    @staticmethod
    def pretty_print_errors(errors: dict) -> None:
        """
        Pretty print schema validation errors in a user-friendly format.

        Args:
            errors (dict): Dictionary of validation errors from Cerberus
        """

        error_messages = []

        def collect_field_errors(field_path, field_errors):
            if isinstance(field_errors, dict):
                for key, value in field_errors.items():
                    new_path = f"{field_path}.{key}" if field_path else key
                    collect_field_errors(new_path, value)
            else:
                field_display = Text()
                field_display.append("Field", style="bold red")
                field_display.append(f" {field_path}", style="yellow")
                if isinstance(field_errors, list):
                    for error in field_errors:
                        field_display.append(f"\n • {error}", style="red")
                    error_messages.append(field_display)
                else:
                    error_messages.append(
                        field_display.append(f"\n • {error}", style="red")
                    )

        collect_field_errors("", errors)

        if error_messages:
            error_text = Text()
            for message in error_messages:
                error_text.append(message)
            panel = Panel(
                error_text,
                title="Configuration Validation Errors",
                border_style="red",
                padding=(1, 2),
            )
            print_stderr(panel)
        else:
            print_stderr("No specific error details available")

    def load(self):
        """
        Load the config file and return the parsed parameters.
        """

        try:

            with open(self.path, mode="r") as f:

                configs = yaml.safe_load_all(
                    f
                )  # Using safe_load_all for untrusted input.

                for config in configs:

                    # Validate the config file schema.
                    if self.valid_schema(document=config):

                        return self.parse(document=self.validator.document)

                    else:

                        logger.error(self.validator.errors)

                        self.pretty_print_errors(self.validator.errors)

                        raise ConfigError("Config file schema is not valid.")

        except FileNotFoundError as e:

            logger.error(f"YAML config file {self.path} not found on filesystem.")

            raise ConfigError(
                f"YAML config file {self.path} not found on filesystem."
            ) from e

        except yaml.YAMLError as e:

            logger.error(f"Error when parsing YAML config file {self.path}")

            raise ConfigError(f"Error when parsing YAML config file {self.path}") from e

        except cerberus.DocumentError as e:

            logger.error(
                f"Error with document when validating YAML config file schema {e}"
            )

            raise ConfigError(
                "Error with document when validating YAML config file schema"
            ) from e

        except UnexpectedInput as e:

            raise ConfigError("Value is not consistent with config DSL.") from e

        except Exception as e:

            if isinstance(e, ConfigError):

                raise

            logger.error(f"An unexpected error has occurred with YAML config file: {e}")

            raise ConfigError(
                "An unexpected error has occurred with YAML config file"
            ) from e
