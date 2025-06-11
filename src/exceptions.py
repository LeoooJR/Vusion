class VCFError(ValueError):
    """
    Exception raised for errors in the VCF file.
    """
    pass

class FastaIndexError(ValueError):
    """
    Exception raised for errors in the FASTA index file.
    """
    pass

class PileupError(ValueError):
    """
    Exception raised for errors in the pileup file.
    """
    pass

class ConfigError(ValueError):
    """
    Exception raised for errors in the config file.
    """
    pass
   

class VariantCallerError(ValueError):
    """
    Exception raised for errors related to the variant caller.
    """
    pass

class CheckSumFileError(ValueError):
    """
    Exception raised for errors related to the checksum file.
    """
    pass

class VariantCallerPluginError(ValueError):
    """
    Exception raised for errors related to the variant caller plugin file.
    """
    pass