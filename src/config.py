import cerberus
from copy import copy
from exceptions import ConfigError
from lark import Lark, Transformer, UnexpectedInput
from loguru import logger
import yaml
from dataclasses import dataclass
from typing import List, Union, Optional

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
                        },
                        "format": {
                            "type": "string",
                            "empty": False,
                            "required": True,
                            # "regex": r"^([A-Z]{1,}:)+[A-Z]{1,}$",
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

        for field in fdocument["caller"]:

            if field == "name":

                pass

            elif field == "info":

                pass

            elif field == "format":

                pass

            elif field in ["genotype", "depth", "vaf"]:

                ast = self.DSL_PARSER.parse(fdocument["caller"][field]["extract"])

                fdocument["caller"][field]["extract"] = TreeToExpression().transform(ast)

            else:

                for subfield in fdocument["caller"][field]:

                    if fdocument["caller"][field][subfield]["extract"]:

                        ast = self.DSL_PARSER.parse(fdocument["caller"][field][subfield]["extract"])

                        fdocument["caller"][field][subfield]["extract"] = TreeToExpression().transform(ast)

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

            logger.error(f"Error with document when validating YAML config file schema {e}")

            raise ConfigError(f"Error with document when validating YAML config file schema {e}")
                
        except UnexpectedInput as e:

            raise ConfigError(f"Value is not consistent with config DSL. {e}")

        except Exception as e:

            if isinstance(e, ConfigError):

                raise

            logger.error(f"An unexpected error has occurred with YAML config file: {e}")

            raise ConfigError(f"An unexpected error has occurred with YAML config file: {e}")