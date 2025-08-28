class FileError(Exception):
    """
    Exception raised for errors in the file.
    """
    pass

class VCFError(FileError):
    """
    Exception raised for errors in the VCF file.
    """
    pass

class FastaIndexError(FileError):
    """
    Exception raised for errors in the FASTA index file.
    """
    pass

class PileupError(FileError):
    """
    Exception raised for errors in the Pileup file.
    """
    pass

class ConfigError(FileError):
    """
    Exception raised for errors in the config file.
    """
    pass
   

class VariantCallerError(Exception):
    """
    Exception raised for errors related to the variant caller.
    """
    pass

class CheckSumFileError(FileError):
    """
    Exception raised for errors related to the checksum file.
    """
    pass

class VariantCallerPluginError(FileError):
    """
    Exception raised for errors related to the variant caller plugin file.
    """
    pass