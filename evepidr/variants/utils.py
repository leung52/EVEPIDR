def mutate_sequence(protein_sequence: str, start: int, end: int, mutation: str) -> str:
    """
    Mutates a given sequence by replacing the part from start to end with mutation.
  
    Parameters:
    protein_sequence (str): The protein sequence to be mutated.
    start (int): The start position of the mutation (1-based index).
    end (int): The end position of the mutation (inclusive, 1-based index).
    mutation (str): The mutation sequence to insert at start position.
  
    Returns:
    str: The mutated sequence.
    """
    # Adjusting for 1-based index
    start_index = start - 1
    end_index = end
  
    # Ensure the start and end positions are within the original sequence
    if start_index < 0 or end_index > len(protein_sequence) or start_index > end_index:
        raise ValueError("Invalid start or end position for the sequence mutation.")
  
    # Replace the sequence between start and end with new_sequence
    mutated_sequence = protein_sequence[:start_index] + mutation + protein_sequence[end_index:]
  
    return mutated_sequence

def three_to_one_aa_code(code: str) -> str:
    """
    """
    three_to_single = {
        "Ala": "A",
        "Arg": "R",
        "Asn": "N",
        "Asp": "D",
        "Cys": "C",
        "Glu": "E",
        "Gln": "Q",
        "Gly": "G",
        "His": "H",
        "Ile": "I",
        "Leu": "L",
        "Lys": "K",
        "Met": "M",
        "Phe": "F",
        "Pro": "P",
        "Ser": "S",
        "Thr": "T",
        "Trp": "W",
        "Tyr": "Y",
        "Val": "V"
    }

    one_code = three_to_single.get(code)
    
    if not one_code:
        raise ValueError(f"{code} is not an amino acid code")
    else:
        return one_code

def adjust_clinvar_classification(clinvar_classification: str) -> str:
    """
    """
    groups = {
        ["Benign", "Likely benign", "protective"]: 'Bengin',
        ["Likely pathogenic", "Pathogenic", "Likely pathogenic, low penetrance", "Pathogenic, low penetrance", "Likely risk allele", "Established risk allele", "association"]: 'Pathogenic'
    }
    for key, value in groups.items():
        if clinvar_classification in key:
            return value
    return "Other"
