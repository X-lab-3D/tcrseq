from tcrseq import get_full_tcr_sequence


def test_full_tcr_sequence():

    full_sequence = get_full_tcr_sequence("TRBV1", "TRBJ1", "CTSSQAEFAFANTEA")
    assert full_sequence.replace(".", "") == "DTGITQTPKYLVTAMGSKRTMKREHLGHDSMYWYRQKAKKSLEFMFYYNCKEFIENKTVPNHFTPECPDSSRLYLHVVALQQEDSAAYLCTSSQAEFAFANTEAFFGQGTRLTVV"

    full_sequence = get_full_tcr_sequence("TRBV28*01", "TRBJ2-1*01", "CASSLTGTGFKQFF")
    assert full_sequence.replace(".", "") == "DVKVTQSSRYLVKRTGEKVFLECVQDMDHENMFWYRQDPGLGLRLIYFSYDVKMKEKGDIPEGYSVSREKKERFSLILESASTNQTSMYLCASSLTGTGFKQFFGPGTRLTVL"

    full_sequence = get_full_tcr_sequence("TRBV26*01", "TRBJ2-3*01", "CASSPGDSAETLYF")
    assert full_sequence.replace(".", "") == "DAVVTQFPRHRIIGTGKEFILQCSQNMNHVTMYWYRQDPGLGLKLVYYSPGTGSTEKGDISEGYHVS*NTIASFPLTLKSASTNQTSVYLCASSPGDSAETLYFGPGTRLTVL"
