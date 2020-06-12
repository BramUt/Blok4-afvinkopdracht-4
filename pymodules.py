import re
from Bio.Blast import NCBIWWW, NCBIXML


def is_dna(seq):
    """Checks for characters that aren't valid Nucleotides for DNA"""
    if re.search("[^ATCG]", seq):
        return False
    else:
        return True


def is_rna(seq):
    """Checks for characters that aren't valid Nucleotides for RNA"""
    if re.search("[^AUCG]", seq):
        return False
    else:
        return True


def is_prot(seq):
    """Checks for any character that aren't valid amino acids"""
    if re.search("[^ACDEFGHIKLMNPQRSTVWY]", seq):
        return False
    else:
        return True


def get_gene(seq):
    """Zoekt met behulp van BLASTp naar het gen van een gegeven eiwitsequentie.

    Input:  seq - str, eiwitsequentie.

    Output: accessie_code - str, accessie code van het eiwit.
    """
    # BLAST met blastp
    output_raw = NCBIWWW.qblast("blastp", "nr", sequence=seq, hitlist_size=1)
    output = list(NCBIXML.parse(output_raw))
    # Eerste record wordt gepakt.
    record = output[0]
    if record.alignments:
        # Zoekt accessiecode in de header.
        accessie_code = re.search(r"(?<=\|).+?(?=\|)",
                                  record.alignments[0].title)[0]
        return accessie_code
    else:
        return "Not found"


# mgarasvlsggeldrwekirlrpggkkkyklkhivwasrelerfavnpglletsegcrqilgqlqpslqtgseelrslyntvatlycvhqrieikdtkealdkieeeqnkskkkaqqaaadtghsnqvsqnypivqniqgqmvhqaisprtlnawvkvveekafspevipmfsalsegatpqdlntmlntvgghqaamqmlketineeaaewdrvhpvhagpiapgqmreprgsdiagttstlqeqigwmtnnppipvgeiykrwiilglnkivrmysptsildirqgpkepfrdyvdrfyktlraeqasqevknwmtetllvqnanpdcktilkalgpaatleemmtacqgvggpghkarvlaeamsqvtnsatimmqrgnfrnqrkivkcfncgkeghtarncraprkkgcwkcgkeghqmkdcterqanflgkiwpsykgrpgnflqsrpeptappeesfrsgvetttppqkqepidkelypltslrslfgndpssq
