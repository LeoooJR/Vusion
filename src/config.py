import cerberus
import jinja2
import os
import yaml

class ConfigParser():

    validator = cerberus.Validator({"caller": {"type": "dict",
                                                "schema": {"name": {"type": "string"}, 
                                                           "format": {"type": "list", 
                                                                      "schema": {"type": "string"}}, 
                                                            "VAF": {}, 
                                                            "depth": {}, 
                                                            "RRC": {}, 
                                                            "ARC": {}
                                                            }}})

    def __init__(self, path: str):

        self.path: str = path

    def load(self):

        with open(self.path, mode='r') as f:

            config = yaml.safe_load_all(f)

        # Path to the template directory
        ressources = os.path.join(
            os.path.dirname(os.path.abspath(__file__)), "templates"
        )

        env = jinja2.Environment(loader=jinja2.FileSystemLoader(ressources))

        # Load the header template
        template = env.get_template("caller")

        with open(f"{config['name']}.py", mode='w') as plugin:

            plugin.writelines(template.render(caller=config['name'], format=config["format"]))

