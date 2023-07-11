from numpy import array

def local_emf():
    # EMF vector at AE building (WMM-2020 model)
    b_emf_wmm2020 = array([0.7216, 19.1629, -45.4592])

    # EMF vector at AE building (IGRF2020 model)
    b_emf_igrf2020 = array([0.7313, 19.1870, -45.4557])

    # Take average of the two as baseline B_EMF:
    b_emf_g = 0.5 * (b_emf_wmm2020 + b_emf_igrf2020)

    return b_emf_g
