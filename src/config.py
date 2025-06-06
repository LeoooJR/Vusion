import cerberus
from copy import copy
from exceptions import ConfigError
from lark import Lark, Transformer, UnexpectedInput
from loguru import logger
import yaml

class TreeToExpression(Transformer):
    
    # def INT(self, token):
    #     return token.update(value=int(token))
    
    def STRING(self, token):

        return str(token.value)
    
    def METADATA(self, token):

        if token.value.isnumeric():

            return int(token.value)
        
        return str(token.value)
    
    OPERATOR = str
    FIELD = str
    indexing = list
    term = list
    expression = list

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
            FIELD: STRING
            indexing: "[" METADATA ("," METADATA)* "]"
            METADATA: INT | STRING
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

                print(TreeToExpression().transform(ast))

                fdocument["caller"][field]["extract"] = TreeToExpression().transform(ast)

            else:

                for subfield in fdocument["caller"][field]:

                    if fdocument["caller"][field][subfield]["extract"]:

                        ast = self.DSL_PARSER.parse(fdocument["caller"][field][subfield]["extract"])

                        fdocument["caller"][field][subfield]["extract"] = TreeToExpression().transform(ast)

        return fdocument
    
    def pretty_print_errors(errors):

        pass

    def load(self):

        with open(self.path, mode="r") as f:

            try:

                configs = yaml.safe_load_all(f)

                for config in configs:

                    if self.valid_schema(document=config):

                        return self.valid_dsl(document=self.validator.document)

                    else:

                        logger.error(self.validator.errors)

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