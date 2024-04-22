import pandas as pd
import requests
import json

from evepidr.variants.utils import mutate_sequence

def get_canonical_sequence_from_uniprot(accession
