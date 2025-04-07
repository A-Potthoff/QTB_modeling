from modelbase.ode import Model

from .rate_laws import *

# # QSSA from matuszynska script
# def ps2states(PQ, PQred, ps2cs, Q, PSIItot, k2, kF, _kH, Keq_PQred, kPQred, pfd, kH0):
#     L = ps2cs * pfd               #? mass_action_2s
#     kH = kH0 + _kH * Q            #???? alm ?
#     k3p = kPQred * PQ                             # tick
#     k3m = kPQred * PQred / Keq_PQred              # tick

#     Bs = []

#     if isinstance(kH, float) and isinstance(PQ, np.ndarray):
#         kH = np.repeat(kH, len(PQ))

#     for L, kH, k3p, k3m in zip(L, kH, k3p, k3m):
#         M = np.array(
#             [
#                 [-L - k3m, kH + kF, k3p, 0],
#                 [L, -(kH + kF + k2), 0, 0],
#                 [0, 0, L, -(kH + kF)],
#                 [1, 1, 1, 1],
#             ]
#         )
#         A = np.array([0, 0, 0, PSIItot])
#         B0, B1, B2, B3 = np.linalg.solve(M, A)
#         Bs.append([B0, B1, B2, B3])
#     return np.array(Bs).T

# def vPS2(B1, k2):
#     """reaction rate constant for photochemistry"""
#     return 0.5 * k2 * B1


# m.add_reaction(
#     rate_name="vPS2",
#     function=vPS2,
#     stoichiometry={"PQ": -1, "H": 2 / m.get_parameter("bH")},
#     modifiers=["B1"],
#     dynamic_variables=["B1"],  # doesn't depend on PQ
#     parameters=["k2"],
# )

# unchanged from matuszynska script
def fluorescence(Q, B0, B2, ps2cs, k2, kF, kH, kH0):
    return (ps2cs * kF * B0) / (kF + k2 + kH * Q) + (ps2cs * kF * B2) / (kF + kH * Q)

def add_PSII(m: Model) -> Model:

    m.add_compounds(["B0", "B1", "B2"]) # in matuszynska script!

    m.add_algebraic_module(
        module_name="B3_alm",
        function=moiety_3,
        compounds=["B0", "B1", "B2"],
        derived_compounds=["B3"],
        parameters=["PSIItot"]
    )

    m.add_algebraic_module(     # unchanged from matuszynska script
        module_name="fluorescence_alm",
        function=fluorescence,
        compounds=["Q", "B0", "B2", "ps2cs"],
        derived_compounds=["Fluo"],
        parameters=["k2", "kF", "kH_factor", "kH0"],
    )


    m.add_reaction_from_args(
        rate_name="vB01",
        function = mass_action_2s,
        stoichiometry={"B0": -1, "B1": +1},
        args = ["B0", "ps2cs", "pfd"] # kLII = ps2cs * pfd  and  vB01 = B0 * kLII 
    )

    m.add_reaction_from_args(
        rate_name="vB10Q",
        function = kquencher,
        stoichiometry={"B1": -1, "B0": +1},
        args=["B1", "Q", "kH_factor", "kH0"]
    )

    m.add_reaction_from_args(
        rate_name="vB10F",
        function=mass_action_1s,
        stoichiometry={"B1": -1, "B0": 1},
        args=["B1", "kF"]
    )

    m.add_reaction_from_args(
        rate_name="vB12",
        function=mass_action_1s,
        stoichiometry={"B1": -1, "B2": 1},
        args=["B1", "k2"] # k2 is "kPchem" or "kP"
    )

    m.add_reaction_from_args(
        rate_name="vB20",
        function=mass_action_22_rev,
        stoichiometry={"B2": -1, "PQ": -0.5, "B0": 1}, # PQred is alm
        args=["B2", "PQ", "PQred", "B0", "kPQred", "Keq_PQred"]
    )

    m.add_reaction_from_args(
        rate_name="vB23",
        function = mass_action_2s,
        stoichiometry={"B2": -1}, # "B3" is alm
        args = ["B2", "ps2cs", "pfd"] # kLII = ps2cs * pfd  and  vB23 = B2 * kLII 
    )

    m.add_reaction_from_args(
        rate_name="vB32F",
        function=mass_action_1s,
        stoichiometry={"B2": 1}, # "B3": -1 (alm)
        args=["B3", "kF"]
    )

    m.add_reaction_from_args(
        rate_name="vB32Q",
        function=kquencher,
        stoichiometry={"B2": 1}, #"B3": -1 (alm)
        args=["B3", "Q", "kH_factor", "kH0"]
    )

    #FIXME: INSERT 0.5 stoichiometries where necessary!

    return m
