import tools.molten_salts as ms


def W(k_d, D, thick, K_S, P_H2):
    """
    Calculates the value of W using the given parameters.

    Parameters:
    k_d (float): The value of dissociation.
    D (float): The value of Diffusivity in the metal.
    thick (float): The value of thickness of the metal.
    K_S (float): The value of the sievert's constant.
    P_H2 (float): The value of P_H2.

    Returns:
    float: The calculated value of W.
    For W<<1, the regime is surface limited.
    """
    return 2 * k_d * thick * P_H2**0.5 / (D * K_S)


def H(k_t, k_H, k_d):
    """
    Calculate the value of H dimensionless number based on the given parameters.

    Parameters:
    k_t (float): The value of mass transport coeff
    k_H (float): The value of Henri's constant.
    k_d (float): The value of dissociation.

    Returns:
    float: The calculated value of H.
    If H>>1, and H/W>>1 the regime is mass transport limited.
    """
    return k_d / (k_t * k_H)


def get_regime(k_d, D, thick, K_S, P_H2, k_t, k_H):
    """
    Get the regime based on the value of H and W.

    Returns:
    str: The regime based on the value of H and W.
    """
    W = ms.W(k_d, D, thick, K_S, P_H2)
    H = ms.H(k_t, k_H, k_d)
    # print("H is equal to", H)
    # print("W is equal to", W)
    # print("H/W is equal to", H / W)
    # if H > 10 and H / W > 10:
    #     # print("Mass transport limited")
    #     return
    # elif H < 0.1 and W < 0.1:
    #     # print("Surface limited")
    #     return
    # elif W > 10 and H / W < 0.1:
    #     # print("Diffusion Limited")
    #     return
    # else:
    #     # print("Mixed regime")
    #     return
