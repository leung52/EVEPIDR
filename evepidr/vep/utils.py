def batchify(sequences: list, batch_size: int) -> list:
  """
  Split sequences into batches of a specified size.

  Parameters:
  sequences (list): List of sequences.
  batch_size (int): Size of each batch.

  Returns:
  List of batches, where each batch is a list of sequences.
  """
  batches = [sequences[i:i + batch_size] for i in range(0, len(sequences), batch_size)]
  return batches
