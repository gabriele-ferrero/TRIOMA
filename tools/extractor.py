# import numpy
# import tools.correlations as cor


# def extractor_lm(Z, R, G_l, G_gas, pl_in, pl_out, T, p_t, K_S):
#     """_summary_

#     Args:
#         Z (_type_): Height
#         R (_type_): Radius
#         G_l (_type_): Liquid flowrate
#         G_gas (_type_): Gas flowrate
#         pl_in (_type_): T pressure inlet
#         pl_out (_type_): T pressure outlet
#         T (_type_): Temperature
#         p_t pressure Pa of the column
#         K_S Sievert's constant

#     Returns:
#         B_l liquid load
#         k_la mass transfer coefficient in packed column
#     """
#     Area = numpy.pi * R**2
#     B_l = (G_l) / Area * 3600 / 0.8314
#     u_l = G_l / Area / 0.8314  # Liquid velocity
#     R_const = 8.31  # Gas constant
#     cg_in = 0
#     R_g = G_gas / G_l
#     pg_in = cg_in * R_const * T
#     A_p = (K_S * p_t * 0.0224 / R_g) ** 0.5 * (
#         pl_in**0.5 - pl_out**0.5 + pg_in * R_g / (K_S * p_t * 0.0224)
#     ) ** 0.5
#     kl_ap = u_l / Z * numpy.log((pl_in**0.5 - A_p) / (pl_out**0.5 - A_p))
#     return [B_l, kl_ap]


# def extractor_ms(Z, R, G_l, G_gas, pl_in, pl_out, T, p_t, K_H):
#     """_summary_

#     Args:
#         Z (_type_): Height
#         R (_type_): Radius
#         G_l (_type_): Liquid flowrate
#         G_gas (_type_): Gas flowrate
#         pl_in (_type_): T pressure inlet
#         pl_out (_type_): T pressure outlet
#         T (_type_): Temperature
#         p_t pressure Pa of the column
#         K_H Henry's constant
#     To be adapted for MS still

#     Returns:
#         B_l liquid load
#         k_la mass transfer coefficient in packed column
#     """
#     Area = numpy.pi * R**2
#     B_l = (G_l) / Area * 3600 / 0.8314
#     u_l = G_l / Area / 0.8314  # Liquid velocity
#     R_const = 8.31  # Gas constant
#     cg_in = 0
#     R_g = G_gas / G_l
#     pg_in = cg_in * R_const * T
#     A_p = (K_H * p_t * 0.0224 / R_g) ** 0.5 * (
#         pl_in**0.5 - pl_out**0.5 + pg_in * R_g / (K_H * p_t * 0.0224)
#     ) ** 0.5
#     kl_ap = u_l / Z * numpy.log((pl_in**0.5 - A_p) / (pl_out**0.5 - A_p))
#     return [B_l, kl_ap]


# def pack_corr(a, d, D, eta, v):
#     k_l = (
#         0.0051
#         * (v / eta / a) ** (2 / 3)
#         * (D / eta) ** 0.5
#         * (a * d) ** 0.4
#         * (1 / eta / 9.81) ** (-1 / 3)
#     )
#     return k_l


# def corr_packed(Re, Sc, d, rho_L, mu_L, L, D):
#     """_summary_
#     Args
#         Re (_type_): Reynolds
#         Sc (_type_): Schmidt
#         d (_type_): ring diameter
#         rho_L (_type_): Liquid density
#         mu_L (_type_): viscosity
#         D diffusion coeff
#         L= characteristic length
#     Returns:
#         _type_: _description_
#         Warning: verification of this must be done
#     """
#     beta = 0.32  # Raschig rings 0.25
#     g = 9.81
#     Sh = beta * Re**0.59 * Sc**0.5 * (d**3 * g * rho_L**2 / mu_L**2) ** 0.17
#     return cor.get_k_from_Sh(Sh, L, D)
