import tools.liquid_metals as lm


def W(k_r, D, thick, K_S, P_H2):
    """
    Calculates the value of W using the given parameters.

    Parameters:
    k_r (float): The value of recombination coefficient.
    D (float): The value of Diffusivity in the metal.
    thick (float): The value of thickness of the metal.
    K_S (float): The value of the sievert's constant.
    P_H2 (float): The value of P_H2.

    Returns:
    float: The calculated value of W.
    For W<<1, the regime is surface limited.
    flux is 0.5*k_d*P_H2^0.5
    """
    return k_r / D * K_S * thick * P_H2**0.5


def partition_param(D, k_t, K_S_S, K_S_L, t):
    """
    Calculate the partition parameter.

    Parameters:
    D (float): Diffusion coefficient.
    k_t (float): mass transfer coeff.
    K_S_S (float): Solid Sievert constant.
    K_S_L (float): Liquid Sievert constant.
    t (float): Thickness.

    Returns:
    float: The partition parameter.
    for partition<<1 membrane effects are slower than mass transfer
    when W<<1 surface limited
    # when W>>1 Diffusion limited limited
    if partition >>1 and W>1 mass transfer limited
    if partition >>1 and W<<1 surface limited and mass transfer limited
    """
    return D / k_t * K_S_S / (K_S_L * t)


def get_regime(D, k_t, K_S_S, K_S_L, k_r, thick, P_H2):
    """
    Get the regime based on the value of H and W.

    Returns:
    str: The regime based on the value of H and W.
    """
    W = lm.W(k_r, D, thick, K_S_S, P_H2)
    partition_param = lm.partition_param(D, k_t, K_S_S, K_S_L, thick)
    print("H is equal to", partition_param, "and W is equal to", W)
    if partition_param > 10 and W > 10:
        print("Mass transport limited")
        return
    elif partition_param < 0.1 and W < 0.1:
        print("Surface limited")
        return
    elif W > 10 and partition_param < 0.1:
        print("Diffusion Limited")
        return
    elif partition_param > 10 and W < 0.1:
        print("Transport and surface limited regime")
        return
    else:
        print("Mixed regime")
        return
