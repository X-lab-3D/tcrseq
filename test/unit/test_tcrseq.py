from tcrseq import get_full_tcr_sequence


def test_full_tcr_sequence():

    full_sequence = get_full_tcr_sequence("TRBV1", "TRBJ1", "CTSSQAEFAFANTEA")

    assert full_sequence.replace(".", "") == "DTGITQTPKYLVTAMGSKRTMKREHLGHDSMYWYRQKAKKSLEFMFYYNCKEFIENKTVPNHFTPECPDSSRLYLHVVALQQEDSAAYLCTSSQAEFAFANTEAFFGQGTRLTVV"
