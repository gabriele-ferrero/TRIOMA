import TRIOMA.tools.molten_salts as ms


def W(k_d, D, thick, K_S, c0, k_H):
    """
    Calculates the value of W using the given parameters.

    Parameters:
    k_d (float): The value of dissociation.
    D (float): The value of Diffusivity in the metal.
    thick (float): The value of thickness of the metal.
    K_S (float): The value of the sievert's constant.
    c0 (float): The value of the initial concentration of hydrogen in the molten salt.
    k_H (float): The value of Henri's constant.

    Returns:
    float: The calculated value of W.
    For W<<1, the regime is surface limited.
    """
    return 2 * k_d * thick * (c0 / k_H) ** 0.5 / (D * K_S)


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


def get_regime(k_d, D, thick, K_S, c0, k_t, k_H, print_var: bool = False):
    """
    Get the regime based on the value of H and W.
    k_d (float): The value of dissociation constant.
    D (float): The value of Diffusivity in the metal.
    thick (float): The value of thickness of the metal.
    K_S (float): The value of the sievert's constant in the metal.
    c0 (float): The value of the initial concentration of hydrogen in the molten salt.
    k_t (float): The value of mass transport coefficient.
    k_H (float): The value of Henri's constant.
    print_var (bool): If True, print the value of H, W and H/W.
    Returns:
    updates the values of Compontent.H and Component.W.
    if print_var is True, it will print the regime based on the value of H and W.
    str: The regime based on the value of H and W.
    """
    W = ms.W(k_d, D, thick, K_S, c0=c0, k_H=k_H)
    H = ms.H(k_t, k_H, k_d)
    result = None  ## initialize the result variable
    match (H, W):
        case (H, W) if H > 10 and H / W > 10:
            result = "Mass transport limited"
        case (H, W) if H < 0.1 and W < 0.1:
            result = "Surface limited"
        case (H, W) if W > 10 and H / W < 0.1:
            result = "Diffusion Limited"
        case _:
            result = "Mixed regime"
    if print_var:
        print("H is equal to", H)
        print("W is equal to", W)
        print("H/W is equal to", H / W)
        print(result)
    return result
