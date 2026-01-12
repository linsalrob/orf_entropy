"""ORF finder wrapper using get_orfs binary."""

import subprocess
import tempfile
from pathlib import Path
from typing import Dict, List, Union

from ..config import DEFAULT_GENETIC_CODE_TABLE, GET_ORFS_BINARY
from ..errors import OrfFinderError
from .types import OrfRecord


def find_orfs(
    sequences: Dict[str, str],
    table_id: int = DEFAULT_GENETIC_CODE_TABLE,
    min_nt_length: int = 90,
    binary_path: str = GET_ORFS_BINARY,
) -> List[OrfRecord]:
    """Find ORFs in DNA sequences using get_orfs binary.
    
    This function wraps the external get_orfs binary (https://github.com/linsalrob/get_orfs).
    The binary must be installed and available in PATH or specified via binary_path.
    
    Args:
        sequences: Dictionary mapping sequence IDs to DNA sequences
        table_id: NCBI genetic code table ID (default: 11, bacterial)
        min_nt_length: Minimum ORF length in nucleotides (default: 90)
        binary_path: Path to get_orfs binary (default: from config/environment)
        
    Returns:
        List of OrfRecord objects
        
    Raises:
        OrfFinderError: If get_orfs binary is not found or fails
    """
    # Check if binary exists
    try:
        result = subprocess.run(
            [binary_path, "-h"],
            capture_output=True,
            text=True,
            timeout=5,
        )
        if result.returncode != 0:
            raise OrfFinderError(f"get_orfs binary check failed: {result.stderr}")
    except FileNotFoundError:
        raise OrfFinderError(
            f"get_orfs binary not found at: {binary_path}\n"
            "Please install from https://github.com/linsalrob/get_orfs\n"
            "or set GET_ORFS_PATH environment variable."
        )
    except subprocess.TimeoutExpired:
        raise OrfFinderError("get_orfs binary check timed out")
    
    # Write sequences to temporary FASTA file
    with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False) as tmp_fasta:
        tmp_fasta_path = tmp_fasta.name
        for seq_id, sequence in sequences.items():
            tmp_fasta.write(f">{seq_id}\n{sequence}\n")
    
    try:
        # Run get_orfs
        # Expected command: get_orfs -f fasta_file -t table_id -m min_length
        cmd = [
            binary_path,
            "-f", tmp_fasta_path,
            "-t", str(table_id),
            "-m", str(min_nt_length),
        ]
        
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=300,  # 5 minute timeout
        )
        
        if result.returncode != 0:
            raise OrfFinderError(f"get_orfs failed: {result.stderr}")
        
        # Parse output
        orfs = _parse_get_orfs_output(result.stdout, sequences, table_id)
        return orfs
        
    finally:
        # Clean up temporary file
        Path(tmp_fasta_path).unlink(missing_ok=True)


def _parse_get_orfs_output(
    output: str, sequences: Dict[str, str], table_id: int
) -> List[OrfRecord]:
    """Parse get_orfs output into OrfRecord objects.
    
    Expected output format (tab-separated):
    seq_id  start  end  strand  frame  has_start  has_stop  sequence
    
    Note: This is a placeholder implementation. The actual parsing depends on
    the output format of get_orfs. This may need to be adjusted.
    
    Args:
        output: stdout from get_orfs
        sequences: Original DNA sequences (for validation)
        table_id: Genetic code table used
        
    Returns:
        List of OrfRecord objects
    """
    orfs = []
    
    for line in output.strip().split("\n"):
        if not line or line.startswith("#"):
            continue
        
        parts = line.split("\t")
        if len(parts) < 7:
            continue
        
        try:
            parent_id = parts[0]
            start = int(parts[1])
            end = int(parts[2])
            strand = parts[3]
            frame = int(parts[4])
            has_start = parts[5].lower() in ("true", "1", "yes")
            has_stop = parts[6].lower() in ("true", "1", "yes")
            
            # Extract sequence
            if parent_id not in sequences:
                continue
            
            parent_seq = sequences[parent_id]
            nt_sequence = parent_seq[start:end]
            
            # Generate unique ORF ID
            orf_id = f"{parent_id}_orf_{start}_{end}_{strand}_f{frame}"
            
            orf = OrfRecord(
                parent_id=parent_id,
                orf_id=orf_id,
                start=start,
                end=end,
                strand=strand,
                frame=frame,
                nt_sequence=nt_sequence,
                table_id=table_id,
                has_start_codon=has_start,
                has_stop_codon=has_stop,
            )
            orfs.append(orf)
            
        except (ValueError, IndexError) as e:
            # Skip malformed lines
            continue
    
    return orfs
