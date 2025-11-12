"""
A module to define the exceptions.
"""


class FileError(Exception):
    """
    Exception raised for errors in the file.
    """


class VCFError(FileError):
    """
    Exception raised for errors in the VCF file.
    """


class FastaIndexError(FileError):
    """
    Exception raised for errors in the FASTA index file.
    """


class PileupError(FileError):
    """
    Exception raised for errors in the Pileup file.
    """


class ConfigError(FileError):
    """
    Exception raised for errors in the config file.
    """


class VariantCallerError(Exception):
    """
    Exception raised for errors related to the variant caller.
    """


class CheckSumFileError(FileError):
    """
    Exception raised for errors related to the checksum file.
    """


class VariantCallerPluginError(FileError):
    """
    Exception raised for errors related to the variant caller plugin file.
    """
