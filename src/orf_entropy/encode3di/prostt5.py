"""ProstT5 encoder for amino acid to 3Di structural token conversion."""

from dataclasses import dataclass
from typing import List, Literal, Optional

try:
    import torch
    from transformers import AutoModel, AutoTokenizer
except ImportError:
    torch = None
    AutoModel = None
    AutoTokenizer = None

from ..config import (
    AUTO_DEVICE,
    CPU_DEVICE,
    CUDA_DEVICE,
    DEFAULT_BATCH_SIZE,
    DEFAULT_PROSTT5_MODEL,
    MPS_DEVICE,
)
from ..errors import DeviceError, EncodingError, ModelError
from ..translate.translator import ProteinRecord


@dataclass
class ThreeDiRecord:
    """Represents a 3Di structural encoding of a protein.
    
    Attributes:
        protein: The ProteinRecord that was encoded
        three_di: The 3Di token sequence
        method: Method used for encoding (always "prostt5_aa2fold")
        model_name: Name of the ProstT5 model used
        inference_device: Device used for inference ("cuda", "mps", or "cpu")
    """

    protein: ProteinRecord
    three_di: str
    method: Literal["prostt5_aa2fold"]
    model_name: str
    inference_device: str


class ProstT5ThreeDiEncoder:
    """Encoder for converting amino acid sequences to 3Di structural tokens.
    
    Uses the ProstT5 model from HuggingFace to predict 3Di tokens directly
    from protein sequences without requiring 3D structures.
    """

    def __init__(
        self,
        model_name: str = DEFAULT_PROSTT5_MODEL,
        device: Optional[str] = None,
    ):
        """Initialize the ProstT5 encoder.
        
        Args:
            model_name: HuggingFace model identifier
            device: Device to use ("cuda", "mps", "cpu", or None for auto-detect)
            
        Raises:
            ModelError: If PyTorch or Transformers are not installed
            DeviceError: If specified device is not available
        """
        if torch is None or AutoModel is None:
            raise ModelError(
                "PyTorch and Transformers are required for 3Di encoding. "
                "Install with: pip install torch transformers"
            )
        
        self.model_name = model_name
        self.device = self._select_device(device)
        self.model = None
        self.tokenizer = None
        
    def _select_device(self, device_hint: Optional[str]) -> str:
        """Select the best available device for inference.
        
        Args:
            device_hint: User-specified device or None for auto-detection
            
        Returns:
            Device string ("cuda", "mps", or "cpu")
            
        Raises:
            DeviceError: If specified device is not available
        """
        # If user specified a device other than "auto", use it
        if device_hint and device_hint != AUTO_DEVICE:
            if device_hint == CUDA_DEVICE and not torch.cuda.is_available():
                raise DeviceError("CUDA requested but not available")
            elif device_hint == MPS_DEVICE and not (
                hasattr(torch.backends, "mps") and torch.backends.mps.is_available()
            ):
                raise DeviceError("MPS requested but not available")
            return device_hint
        
        # Auto-detect best device
        if torch.cuda.is_available():
            return CUDA_DEVICE
        elif hasattr(torch.backends, "mps") and torch.backends.mps.is_available():
            return MPS_DEVICE
        return CPU_DEVICE
    
    def _load_model(self) -> None:
        """Load the ProstT5 model and tokenizer.
        
        Raises:
            ModelError: If model loading fails
        """
        if self.model is not None:
            return  # Already loaded
        
        try:
            self.tokenizer = AutoTokenizer.from_pretrained(self.model_name)
            self.model = AutoModel.from_pretrained(self.model_name)
            self.model.to(self.device)
            self.model.eval()  # Set to evaluation mode
        except Exception as e:
            raise ModelError(f"Failed to load ProstT5 model {self.model_name}: {e}")
    
    def encode(
        self,
        aa_sequences: List[str],
        batch_size: int = DEFAULT_BATCH_SIZE,
    ) -> List[str]:
        """Encode amino acid sequences to 3Di tokens.
        
        Args:
            aa_sequences: List of amino acid sequences
            batch_size: Batch size for processing
            
        Returns:
            List of 3Di token sequences (one per input sequence)
            
        Raises:
            EncodingError: If encoding fails
        """
        self._load_model()
        
        three_di_sequences = []
        
        try:
            # Process in batches
            for i in range(0, len(aa_sequences), batch_size):
                batch = aa_sequences[i:i+batch_size]
                batch_results = self._encode_batch(batch)
                three_di_sequences.extend(batch_results)
        
        except Exception as e:
            raise EncodingError(f"Failed to encode sequences: {e}")
        
        return three_di_sequences
    
    def _encode_batch(self, batch: List[str]) -> List[str]:
        """Encode a batch of sequences.
        
        This is a simplified implementation. The actual ProstT5 model output
        processing may need to be adjusted based on the specific model version.
        
        Args:
            batch: List of amino acid sequences
            
        Returns:
            List of 3Di sequences
        """
        # Tokenize input
        inputs = self.tokenizer(
            batch,
            return_tensors="pt",
            padding=True,
            truncation=True,
            max_length=1024,
        )
        inputs = {k: v.to(self.device) for k, v in inputs.items()}
        
        # Run inference
        with torch.no_grad():
            outputs = self.model(**inputs)
        
        # Extract 3Di predictions
        # Note: This is a placeholder. The actual extraction depends on
        # how ProstT5 outputs 3Di tokens. May need adjustment.
        three_di_results = self._extract_3di_from_outputs(outputs, batch)
        
        return three_di_results
    
    def _extract_3di_from_outputs(
        self, outputs: any, original_sequences: List[str]
    ) -> List[str]:
        """Extract 3Di sequences from model outputs.
        
        Placeholder implementation that returns dummy 3Di sequences.
        This needs to be replaced with actual ProstT5 output processing.
        
        Args:
            outputs: Model outputs
            original_sequences: Original amino acid sequences
            
        Returns:
            List of 3Di sequences
        """
        # TODO: Implement actual 3Di extraction from ProstT5 outputs
        # For now, return placeholder (same length as input, all 'A')
        return ["A" * len(seq) for seq in original_sequences]
    
    def encode_proteins(
        self,
        proteins: List[ProteinRecord],
        batch_size: int = DEFAULT_BATCH_SIZE,
    ) -> List[ThreeDiRecord]:
        """Encode protein records to 3Di records.
        
        Args:
            proteins: List of ProteinRecord objects
            batch_size: Batch size for processing
            
        Returns:
            List of ThreeDiRecord objects
        """
        # Extract sequences
        aa_sequences = [p.aa_sequence for p in proteins]
        
        # Encode
        three_di_sequences = self.encode(aa_sequences, batch_size=batch_size)
        
        # Create records
        records = []
        for protein, three_di in zip(proteins, three_di_sequences):
            record = ThreeDiRecord(
                protein=protein,
                three_di=three_di,
                method="prostt5_aa2fold",
                model_name=self.model_name,
                inference_device=self.device,
            )
            records.append(record)
        
        return records
