import cerberus
import jinja2
from loguru import logger
import os
import yaml

class ConfigParser():

    validator = cerberus.Validator({"caller": {"type": "dict",
                                                "schema": {"name": {"type": "string",
                                                                    "empty": False,
                                                                    "required": True,
                                                                    "forbidden": ["BCFTools", "Varscan", "Vardict", "Pindel", "Haplotypecaller", "Filt3r", "DeepVariant"]},
                                                           "info": {"type": "string",
                                                                    "empty": False,
                                                                    "required": False,
                                                                    "nullable": True,
                                                                    "default": None}, 
                                                           "format": {"type": "string",
                                                                      "empty": False,
                                                                      "required": True,
                                                                      "regex": "^([A-Za-z]+:[A-Za-z]+)+$",
                                                                      "coerce": lambda x: x.replace(':',',')}, 
                                                           "vaf": {"type": "dict",
                                                                    "schema": {
                                                                        "extract": {"type": "string",
                                                                                    "empty": False}
                                                                    },
                                                                    "required": True,
                                                                    "empty": False,
                                                                    "dependencies": "format"}, 
                                                            "depth": {"type": "dict",
                                                                      "schema": {
                                                                        "extract": {"type": "string",
                                                                                    "empty": False}
                                                                    },
                                                                    "required": True,
                                                                    "empty": False,
                                                                    "dependencies": "format"}, 
                                                            "rrc": {"type": "dict",
                                                                    "schema": {
                                                                        "plus": {"type": "string",
                                                                                 "required": False,
                                                                                 "empty": False,
                                                                                 "nullable": True,
                                                                                 "default": None},
                                                                        "minus": {"type": "string",
                                                                                  "required": False,
                                                                                  "empty": False,
                                                                                  "nullable": True,
                                                                                  "default": None},
                                                                        "total": {"type": "string",
                                                                                  "required": False,
                                                                                  "empty": False,
                                                                                  "nullable": True,
                                                                                  "default": None}
                                                                    }, 
                                                                    "required": True,
                                                                    "empty": False,
                                                                    "dependencies": "format"}, 
                                                            "arc": {"required": True, 
                                                                    "empty": False, 
                                                                    "dependencies": "format"}
                                                            }}})

    def __init__(self, path: str):

        self.path: str = path

    def load(self):

        with open(self.path, mode='r') as f:

            config = yaml.safe_load_all(f)

        try:

            if self.validator.validate(config):

                # Path to the template directory
                ressources = os.path.join(
                    os.path.dirname(os.path.abspath(__file__)), "templates"
                )

                env = jinja2.Environment(loader=jinja2.FileSystemLoader(ressources))

                # Load the header template
                template = env.get_template("caller")

                with open(f"{config['name']}.py", mode='w') as plugin:

                    plugin.writelines(template.render(caller=config['name'], format=config["format"]))

            else:

                logger.error(self.validator.errors)

                raise SystemExit("Error when parsing the config file.")
        
        except cerberus.DocumentError as e:

            logger.error(e)

            raise SystemExit(e)

