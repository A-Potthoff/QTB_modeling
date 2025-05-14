from .rate_laws import *
from modelbase.ode import Model

# This is the model of the PSI from Saadat et al. 2021, but this time modeled using ODE instead of QSSA
# Note that this is not the 4-state model of PSI that is used in the other models, here we still model 3 states!


p = {
    "PSItot": 2.5,
    "PCtot": 4.0,  # Bohme1987 but other sources give different values - seems to depend greatly on organism and conditions
    "Fdtot": 5.0,  # Bohme1987
    "kPCox": 2500.0,  # a rough estimate: half life of PC->P700 should be ~0.2ms
    "kFdred": 2.5e5,  # a rough estimate: half life of PC->P700 should be ~2micro-s
    "O2ext": 8.0,  # corresponds to 250 microM cor to 20%

    "E0_PC": 0.380,
    "E0_P700": 0.480,
    "E0_FA": -0.550,
    "E0_Fd": -0.430,
    # physical constants
    "F": 96.485,  # Faraday constant
    "R": 8.3e-3,  # universal gas constant
    "T": 298.0,  # Temperature in K - for now assumed to be constant at 25 C
    # light
    "pfd": 100.0,
    "Ton": 0.0,
    "Toff": 1800,
    "dT": 120,
    "ox": True,  # 1. means True, switched on

    "k": 10.0 ** 8.0 * 8,
}

compounds = ["P700FA", "P700pFAm", "PC", "Fd", "ps2cs"]

def Keq_FAFd(E0_FA, F, E0_Fd, RT):
    DG1 = -E0_FA * F
    DG2 = -E0_Fd * F
    DG = -DG1 + DG2
    K = np.exp(-DG / RT)
    return K


def Keq_PCP700(E0_PC, F, E0_P700, RT):
    DG1 = -E0_PC * F
    DG2 = -E0_P700 * F
    DG = -DG1 + DG2
    K = np.exp(-DG / RT)
    return K

def vPS1(P700FA, ps2cs, pfd):
    """reaction rate constant for open PSI"""
    return (1 - ps2cs) * pfd * P700FA


def get_model():
    m = Model(parameters=p, compounds=compounds)

    m.add_derived_parameter(
        parameter_name="RT",
        function=proportional,
        parameters=["R", "T"],
    )

    # m.add_derived_parameter(
    #     parameter_name="dG_pH",
    #     function=dg_ph,
    #     parameters=["R", "T"],
    # )

    m.add_derived_parameter(
        parameter_name="Keq_FAFd",
        function=Keq_FAFd,
        parameters=["E0_FA", "F", "E0_Fd", "RT"],
    )

    m.add_derived_parameter(
        parameter_name="Keq_PCP700",
        function=Keq_PCP700,
        parameters=["E0_PC", "F", "E0_P700", "RT"],
    )

    m.add_algebraic_module(
        module_name="pc_alm",
        function=moiety_1,
        compounds=["PC"],
        derived_compounds=["PCred"],
        parameters=["PCtot"],
    )

    m.add_algebraic_module(
        module_name="fd_alm",
        function=moiety_1,
        compounds=["Fd"],
        derived_compounds=["Fdred"],
        parameters=["Fdtot"],
    )

    # ! The model

    m.add_algebraic_module(
        module_name="P700pFA_alm",
        function=moiety_2,
        compounds=["P700FA", "P700pFAm"],
        derived_compounds=["P700pFA"],
        parameters=["PSItot"]
    )

    m.add_reaction(
        rate_name="vPS1",
        function=vPS1,
        stoichiometry={"P700FA": -1, "P700pFAm": 1},
        modifiers=["P700FA", "ps2cs"],  # redundant line, does not change the model, tested with and without and simulation results were identical
        dynamic_variables=["P700FA", "ps2cs"],
        parameters=["pfd"],
    )

    m.add_reaction_from_args(
        rate_name = "v4_to_P700pFA",
        function = mass_action_22_rev,
        stoichiometry = {"P700pFAm": -1, "Fd": -1},
        # modifiers = ["Fdred", "P700pFA"],
        # parameters = ["kFdred", "Keq_FAFd"],
        args = ["P700pFAm", "Fd", "P700pFA", "Fdred", "kFdred", "Keq_FAFd"]
    )

    m.add_reaction_from_args(
        rate_name = "v5_to_P700FA",
        function = mass_action_22_rev,
        stoichiometry = {"P700FA": +1, "PC": +1},
        # modifiers = ["PCred", "P700pFA"],
        # parameters=["kPCox", "Keq_PCP700"],
        args = ["P700pFA", "PCred", "P700FA", "PC", "kPCox", "Keq_PCP700"]
    )

    # useful normalizations of concentrations ! 
    
    m.add_algebraic_module(
        module_name="fd_redoxstate",
        function=normalize_concentration,
        compounds=["Fdred"],
        derived_compounds=["Fd_redoxstate"],
        parameters=["Fdtot"],
    )
    
    m.add_algebraic_module(
        module_name="pc_redoxstate",
        function=normalize_concentration,
        compounds=["PCred"],
        derived_compounds=["PC_redoxstate"],
        parameters=["PCtot"],
    )
    
    # PSI

    m.add_algebraic_module(
        module_name="rel_P700pFA_alm",
        function=normalize_concentration,
        compounds=["P700pFA"],
        derived_compounds=["rel_P700pFA"],
        parameters=["PSItot"],
        args = ["P700pFA", "PSItot"]
    )

    m.add_algebraic_module(
        module_name="rel_P700FA_alm",
        function=normalize_concentration,
        compounds=["P700FA"],
        derived_compounds=["rel_P700FA"],
        parameters=["PSItot"],
        args = ["P700FA", "PSItot"]
    )

    m.add_algebraic_module(
        module_name="rel_P700pFAm_alm",
        function=normalize_concentration,
        compounds=["P700pFAm"],
        derived_compounds=["rel_P700pFAm"],
        parameters=["PSItot"],
        args = ["P700pFAm", "PSItot"]
    )

    return m