from .rate_laws import *
from modelbase.ode import Model


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

    "k": 10.0 ** 8.0 * 8,
    "kMehler": 0,
}


compounds = ["PC", "Fd", "ps2cs"]

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

def ps1analytic_mehler(PC, PCred, Fd, Fdred, ps2cs, PSItot, kFdred, KeqF, KeqC, kPCox, pfd, k0, O2):
    """
    QSSA calculates open state of PSI
    depends on reduction states of plastocyanin and ferredoxin
    C = [PC], F = [Fd] (ox. forms)
    y0 = P700FA
    y1 = P700+FA-
    y2 = P700+FA
    """
    kLI = (1 - ps2cs) * pfd

    y0 = (
        KeqC
        * KeqF
        * PCred
        * PSItot
        * kPCox
        * (Fd * kFdred + O2 * k0)
        / (
            Fd * KeqC * KeqF * PCred * kFdred * kPCox
            + Fd * KeqF * kFdred * (KeqC * kLI + PC * kPCox)
            + Fdred * kFdred * (KeqC * kLI + PC * kPCox)
            + KeqC * KeqF * O2 * PCred * k0 * kPCox
            + KeqC * KeqF * PCred * kLI * kPCox
            + KeqF * O2 * k0 * (KeqC * kLI + PC * kPCox)
        )
    )

    y1 = (
        PSItot
        * (Fdred * kFdred * (KeqC * kLI + PC * kPCox) + KeqC * KeqF * PCred * kLI * kPCox)
        / (
            Fd * KeqC * KeqF * PCred * kFdred * kPCox
            + Fd * KeqF * kFdred * (KeqC * kLI + PC * kPCox)
            + Fdred * kFdred * (KeqC * kLI + PC * kPCox)
            + KeqC * KeqF * O2 * PCred * k0 * kPCox
            + KeqC * KeqF * PCred * kLI * kPCox
            + KeqF * O2 * k0 * (KeqC * kLI + PC * kPCox)
        )
    )
    y2 = PSItot - y0 - y1

    return y0, y1, y2

def vFd_red(Fd, Fdred, P700pFAm, P700pFA, kFdred, Keq_FAFd):
    """rate of the redcution of Fd by the activity of PSI
    used to be equall to the rate of PSI but now
    alternative electron pathway from Fd allows for the production of ROS
    hence this rate has to be separate
    """
    return kFdred * Fd * P700pFAm - kFdred / Keq_FAFd * Fdred * P700pFA

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
        module_name="ps1states",
        function=ps1analytic_mehler,
        compounds=["PC", "PCred", "Fd", "Fdred", "ps2cs"], #removed "LHC" (see below)
        derived_compounds=["P700FA", "P700pFAm", "P700pFA"],
        parameters=[
            "PSItot",
            "kFdred",
            "Keq_FAFd",
            "Keq_PCP700",
            "kPCox",
            "pfd",
            "kMehler",
            "O2ext",
        ],
        args=[
            "PC",
            "PCred",
            "Fd",
            "Fdred", 
            #"LHC", #removed cause it is not used in the function
            "ps2cs",
            "PSItot",
            "kFdred",
            "Keq_FAFd",
            "Keq_PCP700",
            "kPCox",
            "pfd",
            "kMehler",
            "O2ext",
        ],
    )

    m.add_reaction(
        rate_name="vPS1",
        function=vPS1,
        stoichiometry={"PC": 1},
        modifiers=["P700FA", "ps2cs"],
        dynamic_variables=["P700FA", "ps2cs"],
        parameters=["pfd"],
    )

    m.add_reaction(
        rate_name="vFdred",
        function=vFd_red,
        stoichiometry={"Fd": -1},
        modifiers=["Fdred", "P700pFAm", "P700pFA"],
        parameters=["kFdred", "Keq_FAFd"],
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
