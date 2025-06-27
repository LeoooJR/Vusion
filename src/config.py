import re
import cerberus
from copy import copy
from exceptions import ConfigError
from lark import Lark, Transformer, UnexpectedInput
from loguru import logger
import yaml
from dataclasses import dataclass
from typing import List, Union, Optional

try:
    from icecream import ic
except ImportError:  # Graceful fallback if IceCream isn't installed.
    ic = lambda *a: None if not a else (a[0] if len(a) == 1 else a)  # noqa

@dataclass
class Metadatas:
    header: Optional[str] = None
    index: Optional[int] = None
    unit: Optional[str] = None

@dataclass
class Term:
    field: str
    metadata: Optional[Metadatas] = None

@dataclass
class Expression:
    terms: List[Union[Term, 'Expression']]
    operator: Optional[str] = None

class TreeToExpression(Transformer):
    def STRING(self, token):
        return str(token.value)
    
    def INDEX(self, token):
        return int(token.value)
    
    def HEADER(self, token):
        return str(token.value)
    
    def UNIT(self, token):
        return str(token.value)
    
    def OPERATOR(self, token):
        return str(token.value)
    
    def FIELD(self, token):
        return str(token.value)
    
    def metadata(self, item):
        return item[0] 
    
    def indexing(self, items: list):
        header = None
        index = None
        unit = None
        
        for item in items:
            if isinstance(item, int):
                index = item
            elif isinstance(item, str) and item == '%':
                unit = item
            else:
                header = item
        
        return Metadatas(header=header, index=index, unit=unit)
    
    def term(self, items):
        if len(items) == 1:
            return Term(field=items[0])
        return Term(field=items[0], metadata=items[1])
    
    def expression(self, items):
        if len(items) == 1:
            return items[0]
            
        terms = []
        current_operator = None
        
        for item in items:
            if isinstance(item, str) and item in ['+', '-', '*', '/']:
                current_operator = item
            else:
                if current_operator and terms:
                    terms.append(Expression(terms=[terms[-1], item], operator=current_operator))
                    current_operator = None
                else:
                    terms.append(item)
                    
        return terms[-1] if len(terms) > 1 else terms[0]
    
class ExpressionVisitor:
    """Class to handle visiting of expression."""

    def visit_EXPRESSION(self, expression) -> dict:
        if isinstance(expression, Term):
            if isinstance(expression.field, Expression):
                return {
                    "type": "expression",
                    "operator": expression.field.operator,
                    "terms": [self.visit_EXPRESSION(term) for term in expression.field.terms]
                }
            else:
                return {
                    "type": "term",
                    "field": expression.field,
                    "metadata": {
                        "header": expression.metadata.header if expression.metadata else None,
                        "index": expression.metadata.index if expression.metadata else None,
                        "unit": expression.metadata.unit if expression.metadata else None
                    } if expression.metadata else None
                }
        
        return {
            "type": "expression",
            "operator": expression.operator,
            "terms": [self.visit_EXPRESSION(term) for term in expression.terms]
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

        return getattr(self, "_format_fields", None)

    @format.setter
    def format(self, value):

        self._format_fields = value

    @property
    def infos(self):

        return getattr(self, '_infos_field', None)
    
    @infos.setter
    def infos(self, value):

        self._infos_fields = value

    def _is_valid_field(self, field: list) -> bool:

        is_in_format: bool = False
        is_in_info: bool = False

        if self._format_fields:
            is_in_format = field in self._format_fields
        
        if self._infos_fields:
            is_in_info = field in self._infos_fields

        return is_in_format or is_in_info
    
    def _is_valid_metadata(self, field, index: int = None, header: str = None, unit: str = None) -> bool:

        is_valid_index: bool = False if index else True
        is_valid_header: bool = False if header else True
        is_valid_unit: bool = False if unit else True

        if index:

            is_valid_index = True

        if header:

            if header == "format":
                is_valid_header = field in self._format_fields if self._format_fields else False 
            else:
                is_valid_header = field in self._infos_fields if self._infos_fields else False

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
                raise UnexpectedInput(f"Key {term['field']} not in FORMAT or INFOS fields")
            
        result = term["field"]
        
        # Add metadata if present
        if term["metadata"]:
            metadata = []
            if term["metadata"]["index"] is not None:
                if not self._is_valid_metadata(field=term["field"], index=term["metadata"]["index"]):
                    raise UnexpectedInput(f"Index metadata in term {term} is not valid.")
                metadata.append(str(term["metadata"]["index"]))
            if term["metadata"]["header"]:
                if not self._is_valid_metadata(field=term["field"], header=term["metadata"]["header"]):
                    raise UnexpectedInput(f"{term["field"]} is not in {term["metadata"]["header"].upper()}.")
                metadata.append(term["metadata"]["header"])
            if term["metadata"]["unit"]:
                if not self._is_valid_metadata(field=term["field"], unit=term["metadata"]["unit"]):
                    raise UnexpectedInput(f"Unit metadata in term {term} is not valid.")
                metadata.append(term["metadata"]["unit"])
                
            if metadata:
                result += f"[{','.join(metadata)}]"
                
        return result
    
    def _format_expression(self, expr: dict) -> str:
        """Format an expression into a template string.
        
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
        """Convert an expression structure to a template string.
        
        Args:
            expression (dict): Expression dictionary from ExpressionVisitor
            
        Returns:
            str: Template string
        """
        return self._format_expression(expression)
    
    def to_python(self, expression: dict) -> str:
        """Convert an expression structure to Python code.
        
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

    SCHEMA = {
                "caller": {
                    "type": "dict",
                    "schema": {
                        "name": {
                            "type": "string",
                            "empty": False,
                            "required": True,
                            "forbidden": [
                                "BCFTools",
                                "Varscan",
                                "Vardict",
                                "Pindel",
                                "Haplotypecaller",
                                "Filt3r",
                                "DeepVariant",
                            ],
                        },
                        "info": {
                            "type": "string",
                            "empty": False,
                            "required": False,
                            "nullable": True,
                            "default": None,
                            "coerce": lambda x: x.replace(";", ",")
                        },
                        "format": {
                            "type": "string",
                            "empty": False,
                            "required": True,
                            # "regex": "^[A-Z]{1,}(:[A-Z]{1,})*$",
                            "coerce": lambda x: x.replace(":", ","),
                        },
                        "genotype": {
                            "type": "dict",
                            "schema": {
                                "extract": {"type": "string", "empty": False, "nullable": False}
                            },
                            "required": True,
                            "nullable": False,
                            "empty": False,
                            "dependencies": "format"
                        },
                        "vaf": {
                            "type": "dict",
                            "schema": {
                                "extract": {"type": "string", "empty": False, "nullable": False}
                            },
                            "required": True,
                            "empty": False,
                            "dependencies": "format",
                        },
                        "depth": {
                            "type": "dict",
                            "schema": {
                                "extract": {"type": "string", "empty": False, "nullable": False}
                            },
                            "required": True,
                            "empty": False,
                            "dependencies": "format",
                        },
                        "rrc": {
                            "type": "dict",
                            "schema": {
                                "forward": {
                                    "type": "dict",
                                    "schema": {
                                        "extract": {
                                            "type": "string",
                                            "required": True,
                                            "empty": False,
                                            "nullable": True,
                                            "default": None
                                        }
                                    },
                                    "empty": False,
                                    "nullable": False,
                                    "required": False
                                },
                                "reverse": {
                                    "type": "dict",
                                    "schema": {
                                        "extract": {
                                            "type": "string",
                                            "required": True,
                                            "empty": False,
                                            "nullable": True,
                                            "default": None
                                        }
                                    },
                                    "empty": False,
                                    "nullable": False,
                                    "required": False
                                },
                                "total": {
                                    "type": "dict",
                                    "schema": {
                                        "extract": {
                                            "type": "string",
                                            "required": True,
                                            "empty": False,
                                            "nullable": True,
                                            "default": None
                                        }
                                    },
                                    "empty": False,
                                    "nullable": False,
                                    "required": False
                                },
                            },
                            "required": True,
                            "empty": False,
                            "nullable": True,
                            "dependencies": "format",
                        },
                        "arc": {
                            "type": "dict",
                            "schema": {
                                "forward": {
                                    "type": "dict",
                                    "schema": {
                                        "extract": {
                                            "type": "string",
                                            "required": True,
                                            "empty": False,
                                            "nullable": True,
                                            "default": None
                                        }
                                    },
                                    "empty": False,
                                    "nullable": False,
                                    "required": False
                                },
                                "reverse": {
                                    "type": "dict",
                                    "schema": {
                                        "extract": {
                                            "type": "string",
                                            "required": True,
                                            "empty": False,
                                            "nullable": True,
                                            "default": None
                                        }
                                    },
                                    "empty": False,
                                    "nullable": False,
                                    "required": False
                                },
                                "total": {
                                    "type": "dict",
                                    "schema": {
                                        "extract": {
                                            "type": "string",
                                            "required": True,
                                            "empty": False,
                                            "nullable": True,
                                            "default": None
                                        }
                                    },
                                    "empty": False,
                                    "nullable": False,
                                    "required": False
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
        
    DSL_PARSER = Lark(GRAMMAR, start="expression", parser="lalr")

    def __init__(self, path: str):

        self.path: str = path

        self.validator = cerberus.Validator(
            self.SCHEMA
        )

    def valid_schema(self, document):

        return self.validator.validate(document)
        
    def valid_dsl(self, document):

        fdocument = copy(document)

        transformer = TreeToExpression()

        visitor = ExpressionVisitor()

        formatter = ExpressionTemplate()

        for field in fdocument["caller"]:

            if field == "name":

                if not fdocument["caller"][field].isidentifier():

                    raise UnexpectedInput("Name value is not a valid identifier.")

            elif field == "info":

                if not re.match(r"^[A-Z]{1,}(,[A-Z]{1,})*$", fdocument["caller"][field]):

                    raise UnexpectedInput("Info value is not consistent with requested format.")
                
                formatter.infos = fdocument["caller"][field].split(',')

            elif field == "format":

                if not re.match(r"^[A-Z]{1,}(,[A-Z]{1,})*$", fdocument["caller"][field]):

                    raise UnexpectedInput("Format value is not consistent with VCF format.")
                
                formatter.format = fdocument["caller"][field].split(',')

            elif field in ["genotype", "depth", "vaf"]:

                ast = self.DSL_PARSER.parse(fdocument["caller"][field]["extract"])

                fdocument["caller"][field]["extract"] = transformer.transform(ast)

                formatter.to_template(visitor.visit_EXPRESSION(fdocument["caller"][field]["extract"]))

            else:

                if fdocument["caller"][field]:

                    for subfield in ["forward", "reverse", "total"]:

                        if subfield in fdocument["caller"][field]:

                            if fdocument["caller"][field][subfield]["extract"]:

                                ast = self.DSL_PARSER.parse(fdocument["caller"][field][subfield]["extract"])

                                fdocument["caller"][field][subfield]["extract"] = transformer.transform(ast)

                                formatter.to_template(visitor.visit_EXPRESSION(fdocument["caller"][field][subfield]["extract"]))

                        else:

                            fdocument["caller"][field][subfield] = {"extract": None}

                else:

                    fdocument["caller"][field] = {"forward": {"extract": None}, "reverse": {"extract": None}, "total": {"extract": None}}

        return fdocument
    
    @staticmethod
    def pretty_print_errors(errors):
        """Pretty print schema validation errors in a user-friendly format.
        
        Args:
            errors (dict): Dictionary of validation errors from Cerberus
        """
        print("\nConfiguration Validation Errors:")
        print("=" * 40)
        
        def print_field_errors(field_path, field_errors):
            if isinstance(field_errors, dict):
                for key, value in field_errors.items():
                    new_path = f"{field_path}.{key}" if field_path else key
                    print_field_errors(new_path, value)
            else:
                print(f"\nField: {field_path}")
                if isinstance(field_errors, list):
                    for error in field_errors:
                        print(f"  - {error}\n")
                else:
                    print(f"  - {field_errors}")
        
        print_field_errors("", errors)
        print("\n" + "=" * 40)

    def load(self):

        try:

            with open(self.path, mode="r") as f:

                configs = yaml.safe_load_all(f)

                for config in configs:

                    if self.valid_schema(document=config):

                        return self.valid_dsl(document=self.validator.document)

                    else:

                        logger.error(self.validator.errors)

                        self.pretty_print_errors(self.validator.errors)

                        raise ConfigError("Config file schema is not valid.")
                
        except FileNotFoundError:

            logger.error(f"YAML config file {self.path} not found on filesystem.")

            raise ConfigError(f"YAML config file {self.path} not found on filesystem.")
                
        except yaml.YAMLError as e:

            logger.error(f"Error when parsing YAML config file {self.path}")

            raise ConfigError(f"Error when parsing YAML config file {self.path}")
                
        except cerberus.DocumentError as e:

            logger.error(f"Error with document when validating YAML config file schema {ic.format(e)}")

            raise ConfigError(f"Error with document when validating YAML config file schema {ic.format(e)}")
                
        except UnexpectedInput as e:

            raise ConfigError(f"Value is not consistent with config DSL. {ic.format(e)}")

        except Exception as e:

            if isinstance(e, ConfigError):

                raise

            logger.error(f"An unexpected error has occurred with YAML config file: {ic.format(e)}")

            raise ConfigError(f"An unexpected error has occurred with YAML config file: {ic.format(e)}")