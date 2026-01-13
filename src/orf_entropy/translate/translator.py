"""Translation of nucleotide sequences to amino acids."""

from dataclasses import dataclass
from typing import List

try:
    from genetic_codes import GeneticCode
except ImportError:
    # If pygenetic-code is not installed, provide a fallback or raise a helpful error
    GeneticCode = None

from ..config import DEFAULT_GENETIC_CODE_TABLE
from ..errors import TranslationError
from ..orf.types import OrfRecord


@dataclass
class ProteinRecord:
    """Represents a translated protein from an ORF.
    
    Attributes:
        orf: The OrfRecord that was translated
        aa_sequence: The amino acid sequence
        aa_length: Length of the amino acid sequence
    """

    orf: OrfRecord
    aa_sequence: str
    aa_length: int
    
    def __post_init__(self) -> None:
        """Validate protein attributes."""
        if self.aa_length != len(self.aa_sequence):
            raise ValueError(
                f"aa_length {self.aa_length} doesn't match sequence length "
                f"{len(self.aa_sequence)}"
            )


def translate_orf(orf: OrfRecord, table_id: int = DEFAULT_GENETIC_CODE_TABLE) -> ProteinRecord:
    """Translate an ORF to a protein sequence.
    
    Uses the pygenetic-code library for translation with NCBI genetic codes.
    Ambiguous codons (containing N or other IUPAC codes) are translated to 'X'.
    
    Args:
        orf: OrfRecord to translate
        table_id: NCBI genetic code table ID (default: from config)
        
    Returns:
        ProteinRecord with translated sequence
        
    Raises:
        TranslationError: If translation fails
    """
    if GeneticCode is None:
        raise TranslationError(
            "pygenetic-code package not installed. "
            "Install with: pip install pygenetic-code"
        )
    
    try:
        # Load genetic code
        genetic_code = GeneticCode.from_ncbi_id(table_id)
        
        # Translate sequence
        # Note: genetic_code.translate() should handle ambiguous codons
        aa_sequence = genetic_code.translate(orf.nt_sequence)

        if aa_sequence != orf.aa_sequence:
            raise ValueError(f"Translated amino acid sequence is not the same as read from the file:\nProvided:\n{orf.aa_sequence}\nTranslated:\n{aa_sequence}")
        
        # Remove stop codon (*) if present at the end
        if aa_sequence.endswith("*"):
            aa_sequence = aa_sequence[:-1]
        
        return ProteinRecord(
            orf=orf,
            aa_sequence=aa_sequence,
            aa_length=len(aa_sequence),
        )
        
    except Exception as e:
        raise TranslationError(f"Failed to translate ORF {orf.orf_id}: {e}")


def translate_orfs(
    orfs: List[OrfRecord], table_id: int = DEFAULT_GENETIC_CODE_TABLE
) -> List[ProteinRecord]:
    """Translate multiple ORFs to protein sequences.
    
    Args:
        orfs: List of OrfRecord objects to translate
        table_id: NCBI genetic code table ID
        
    Returns:
        List of ProteinRecord objects
    """
    return [translate_orf(orf, table_id=table_id) for orf in orfs]
