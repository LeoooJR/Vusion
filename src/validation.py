from argparse import Action

class ValidateThresholdsAction(Action):
    """
    Action to validate the thresholds in --thresholds option.
    """

    def __call__(self, parser, namespace, values, option_string=None):
        """
        Validate the thresholds in --thresholds option.
        
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
        
        # Retrieve the thresholds from the command line as a list of strings
        thresholds: list[str] = values.split(',')

        # Check that we have 10 values in the thresholds
        size: int = len(thresholds)
        if size != 10:
            # Raise an error if the number of values is not 10
            parser.error(f"Invalid number of values in --threshold option: {size}. Expected 10 values.")
            
        # Check that all values can be converted to floats
        try:
            # Convert the values to floats
            thresholds: list[float] = list(map(float, thresholds)) 
        except ValueError as e:
            # Raise an error if the values cannot be converted to floats
            parser.error(f"Invalid values in --thresholds option: {e}. Expected floats.")
            
        # Check that all values are between 0 and 100
        if not any(list(map(lambda value: value >= MIN_THRESHOLD and value <= MAX_THRESHOLD, thresholds))):
            # Raise an error if the values are not between 0 and 100
            parser.error(f"Option --thresholds values cannot be lower than {MIN_THRESHOLD} or higher than {MAX_THRESHOLD}.")

        # Check that 6 first given threshold are unique
        if len(thresholds[0:6]) != len(set(thresholds[0:6])):
            # Raise an error if the values are not unique
            parser.error(f"Option --thresholds six first values must be unique.")

        # Check that second group of threshold values are unique
        if len(thresholds[6:9]) != len(set(thresholds[6:9])):
            # Raise an error if the values are not unique
            parser.error(f"Option --thresholds values 7, 8 and 9 must be unique.")
        
        # Sort the first 6 values and the last 3 values
        # This is done to make sure that the values are in the right order
        thresholds[0:6] = sorted(thresholds[0:6])
        thresholds[6:9] = sorted(thresholds[6:9])

        # Store the thresholds in the namespace object
        setattr(namespace, self.dest, thresholds)

class ValidateVCFSAction(Action):
    """
    Action to validate the VCF input in --vcf option.
    """

    def __call__(self, parser, namespace, values, option_string = None):
        """
        Validate the VCF input in --vcf option.
        
        Args:
            parser: The argument parser.
            namespace: The namespace to store the parsed value.
            values: The value to validate.
            option_string: The option string that triggered this action.
            
        Raises:
            ArgumentTypeError: If the value is not consistent with required format.
        """

        # Retrieve the VCF input from the command line as a list of strings
        infs: list[str] = values[0].split(',')

        size: int = len(infs)
        # Check that the number of arguments is 2 or 3
        if not size in [2, 3]:
            # Raise an error if the number of arguments is not 2 or 3
            parser.error(f"Wrong number of argument in --vcf option {values}: {size}. Expected 2 or 3 values.")
            
        # Check that the first argument is alpha only
        if not infs[0].isalpha():
            # Raise an error if the first argument is not alpha only
            parser.error(f"The Variant Caller identifier in --vcf option must consist of letters only: {values}.")

        # Retrieve the VCFs from the namespace object as a list
        vcfs = getattr(namespace, self.dest)

        # If previous VCFs are already set
        if vcfs:
            # Append the new VCFs to the list of VCFs
            vcfs.append(infs)

        # If no previous VCFs are set
        else:
            # Set the VCFs to the list of VCFs
            setattr(namespace, self.dest, [infs])

class ValidateIndelsLengthAction(Action):
    """
    Action to validate the minimum INDELS length value in --length-indels option.
    """

    def __call__(self, parser, namespace, values, option_string = None):
        """
        Validate the minimum INDELS length value in --length-indels option.
        
        Args:
            parser: The argument parser.
            namespace: The namespace to store the parsed value.
            values: The value to validate.
            option_string: The option string that triggered this action.
            
        Raises:
            ArgumentTypeError: If the value is not a positive unsigned integer.
        """

        try:
            # Retrieve the minimum INDELS length value from the command line as a float
            l: float = float(values)

            # Check that the minimum INDELS length value is a positive unsigned integer
            if l <= 0.0:

                # Raise an error if the minimum INDELS length value is not a positive unsigned integer
                parser.error(f"Minimum INDELS length value must be a positive unsigned integer: {l}")
                
            # Store the minimum INDELS length value in the namespace object
            setattr(namespace, self.dest, l)
        
        # Catch an error if the value cannot be converted to a float
        except ValueError as e:

            # Raise an error if the minimum INDELS length value is not a positive unsigned integer
            parser.error(f"Minimum INDELS length value must be a positive unsigned integer: {e}")

class ValidateSBMAction(Action):
    """
    Action to validate the strand bias metric limit value in --sbm-homozygous option.
    """

    def __call__(self, parser, namespace, values, option_string = None):
        """
        Validate the strand bias metric limit value in --sbm-homozygous option.
        
        Args:
            parser: The argument parser.
            namespace: The namespace to store the parsed value.
            values: The value to validate.
            option_string: The option string that triggered this action.
            
        Raises:
            ArgumentTypeError: If the value is not a positive unsigned value.
        """

        try:
            # Retrieve the strand bias metric limit value from the command line as a float
            sbm: float = float(values)

            # Check that the strand bias metric limit value is a positive unsigned value
            if sbm <= 0.0:

                # Raise an error if the strand bias metric limit value is not a positive unsigned value
                parser.error(f"Strand biais metric limit must be a positive unsigned value: {sbm}")

            # Store the strand bias metric limit value in the namespace object
            setattr(namespace, self.dest, sbm)

        # Catch an error if the value cannot be converted to a float
        except ValueError as e:

            # Raise an error if the strand bias metric limit value is not a positive unsigned value
            parser.error(f"Strand biais metric limit must be a positive unsigned value: {e}")