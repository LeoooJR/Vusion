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

class VariantCallerError(ValueError):
    """
    Exception raised for errors related to the variant caller.
    """
    pass