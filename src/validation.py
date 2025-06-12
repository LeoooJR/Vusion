from argparse import Action

class ValidateThresholdsAction(Action):

    def __call__(self, parser, namespace, values, option_string=None):
        """Validate the thresholds in --thresholds option.
        
        Args:
            parser: The argument parser.
            namespace: The namespace to store the parsed value.
            values: The value to validate.
            option_string: The option string that triggered this action.
            
        Raises:
            ArgumentTypeError: If the value is not consistent with required format.
        """

        # Maximum and minimum threshold values
        MAX_THRESHOLD: float = 100.0
        MIN_THRESHOLD: float = 0.0
        
        thresholds: list[str] = values.split(',')

        # Check that we have 10 values in thresholds
        if len(thresholds) != 10:
            parser.error(f"Invalid number of values in --threshold option.")
            
        #Check that all values can be converted to floats
        try:
            thresholds: list[float] = list(map(float, thresholds)) 
        except ValueError:
            parser.error(f"Invalid values in --thresholds option.")
            
        if not any(list(map(lambda value: value >= MIN_THRESHOLD and value <= MAX_THRESHOLD, thresholds))):
            parser.error(f"Option --thresholds values cannot be equal or higher than 100, or lower than 0.")

        # Check that 6 first given threshold are unique
        if len(thresholds[0:6]) != len(set(thresholds[0:6])):
            parser.error(f"Option --thresholds six first values must be unique.")

        # Check that second group of threshold values are unique
        if len(thresholds[6:9]) != len(set(thresholds[6:9])):
            parser.error(f"Option --thresholds values 7, 8 and 9 must be unique.")
        
        # Sort the first 6 values and the last 3 values
        # This is done to make sure that the values are in the right order
        thresholds[0:6] = sorted(thresholds[0:6])
        thresholds[6:9] = sorted(thresholds[6:9])

        setattr(namespace, self.dest, thresholds)

class ValidateVCFSAction(Action):

    def __call__(self, parser, namespace, values, option_string = None):
        """Validate the VCF input in --vcf option.
        
        Args:
            parser: The argument parser.
            namespace: The namespace to store the parsed value.
            values: The value to validate.
            option_string: The option string that triggered this action.
            
        Raises:
            ArgumentTypeError: If the value is not consistent with required format.
        """

        infs: list[str] = values[0].split(',')

        if not len(infs) in [2, 3]:
            parser.error(f"Wrong number of argument in --vcf option {values}.")
            
        if not infs[0].isalpha():
            parser.error(f"The Variant Caller identifier in --vcf option must consist of letters only {values}.")

        vcfs = getattr(namespace, self.dest)

        if vcfs:

            vcfs.append(infs)

        else:

            setattr(namespace, self.dest, [infs])

class ValidateIndelsLengthAction(Action):

    def __call__(self, parser, namespace, values, option_string = None):
        """Validate the minimum INDELS length value.
        
        Args:
            parser: The argument parser.
            namespace: The namespace to store the parsed value.
            values: The value to validate.
            option_string: The option string that triggered this action.
            
        Raises:
            ArgumentTypeError: If the value is not a positive unsigned integer.
        """

        try:

            l = float(values)

            if l <= 0.0:

                parser.error("Minimum INDELS length value must be a positive unsigned integer.")
            
            setattr(namespace, self.dest, l)
            
        except ValueError:

            parser.error("Minimum INDELS length value must be a positive unsigned integer.")

class ValidateSBMAction(Action):

    def __call__(self, parser, namespace, values, option_string = None):
        """Validate the strand bias metric limit value.
        
        Args:
            parser: The argument parser.
            namespace: The namespace to store the parsed value.
            values: The value to validate.
            option_string: The option string that triggered this action.
            
        Raises:
            ArgumentTypeError: If the value is not a positive unsigned value.
        """

        try:

            sbm = float(values)

            if sbm <= 0.0:

                parser.error("Strand biais metric limit must be a positive unsigned value.")

            setattr(namespace, self.dest, sbm)

        except ValueError:

            parser.error("Strand biais metric limit must be a positive unsigned value.")