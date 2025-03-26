from modelbase.ode import Model

from .rate_laws import vPS1


def ps1analytic_mehler(PC, PCred, Fd, Fdred, LHC, ps2cs, PSItot, kFdred, KeqF, KeqC, kPCox, pfd, k0, O2):
    """
    QSSA calculates open state of PSI
    depends on reduction states of plastocyanin and ferredoxin
    C = [PC], F = [Fd] (ox. forms)
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


def vFd_red(Fd, Fdred, A1, A2, kFdred, Keq_FAFd):
    """rate of the redcution of Fd by the activity of PSI
    used to be equall to the rate of PSI but now
    alternative electron pathway from Fd allows for the production of ROS
    hence this rate has to be separate
    """
    return kFdred * Fd * A1 - kFdred / Keq_FAFd * Fdred * A2


def vAscorbate(A, H, kf1, kr1, kf2, kr2, kf3, kf4, kr4, kf5, XT):
    """lumped reaction of ascorbate peroxidase
    the cycle stretched to a linear chain with
    two steps producing the MDA
    two steps releasing ASC
    and one step producing hydrogen peroxide"""
    nom = A * H * XT
    denom = (
        A * H * (1 / kf3 + 1 / kf5)
        + A / kf1
        + H / kf4
        + H * kr4 / (kf4 * kf5)
        + H / kf2
        + H * kr2 / (kf2 * kf3)
        + kr1 / (kf1 * kf2)
        + kr1 * kr2 / (kf1 * kf2 * kf3)
    )
    return nom / denom


def vMDAreduct(NADPH, MDA, kcatMDAR, KmMDAR_NADPH, KmMDAR_MDA, MDAR0):
    """Compare Valero et al. 2016"""
    nom = kcatMDAR * MDAR0 * NADPH * MDA
    denom = KmMDAR_NADPH * MDA + KmMDAR_MDA * NADPH + NADPH * MDA + KmMDAR_NADPH * KmMDAR_MDA
    return nom / denom


def vMehler(A, O2ext, kMehler):
    """Draft Mehler reaction inspired from PSI reaction.
    This reaction is lumping the reduction of O2 instead of Fd
    resulting in Superoxide, as well as the Formation of H2O2 in one reaction.
    The entire reaction is scaled by the arbitrary parameter kMehler
    """
    return A * kMehler * O2ext


def vGR(NADPH, GSSG, kcat_GR, GR0, KmNADPH, KmGSSG):
    nom = kcat_GR * GR0 * NADPH * GSSG
    denom = KmNADPH * GSSG + KmGSSG * NADPH + NADPH * GSSG + KmNADPH * KmGSSG
    return nom / denom


def vDHAR(DHA, GSH, kcat_DHAR, DHAR0, KmDHA, K, KmGSH):
    nom = kcat_DHAR * DHAR0 * DHA * GSH
    denom = K + KmDHA * GSH + KmGSH * DHA + DHA * GSH
    return nom / denom


def v3ASC(MDA, k3):
    return k3 * MDA ** 2


def ascorbate_moiety(MDA, DHA, ASCtotal):
    return ASCtotal - MDA - DHA


def glutathion_moiety(GSSG, GStotal):
    return GStotal - 2 * GSSG


def add_mehler(m) -> Model:
    m.add_parameters(
        {
            "kf1": 10000.0,
            "kr1": 220.0,
            "kf2": 10000.0,
            "kr2": 4000.0,
            "kf3": 2510.0,  #
            "kf4": 10000.0,
            "kr4": 4000.0,
            "kf5": 2510.0,  #
            "XT": 0.07,  # according to Valero
            "kMehler": 1.0,
            # V09 µM -> mM
            "kcat_GR": 595,
            "kcat_DHAR": 142,
            "k1APX": 12 / 1e-3,
            "k2APX": 50 / 1e-3,
            "k3APX": 2.1 / 1e-3,
            "k4APX": 0.7 / 1e-3,
            "k5APX": 0.01,
            "k3": 0.5 / 1e-3,
            "k4": 0.1 / 1e-3,
            "k5": 0.2 / 1e-3,
            "k6": 0.2 / 1e-3,
            "k7": 0.7 / 1e-3,
            "k8": 2e-6 / 1e-3,
            "KmNADPH": 3e-3,
            "KmGSSG": 2e2 * 1e-3,
            "KmDHA": 70e-3,
            "KmGSH": 2.5e3 * 1e-3,
            "K": 5e5 * (1e-3) ** 2,  # ?
            "GR0": 1.4e-3,
            "DHAR0": 1.7e-3,
            "APX0": 70e-3,
            "Glutathion_total": 10,
            "Ascorbate_total": 10,
            # V16 µM->mM and h->s
            "kcatMDAR": 1080000 / (60 * 60),
            "KmMDAR_NADPH": 23e-3,
            "KmMDAR_MDA": 1.4e-3,
            "MDAR0": 2e-3,
        }
    )
    m.add_compounds(
        [
            "MDA",
            "H2O2",
            "DHA",
            "GSSG",
        ]
    )

    m.add_algebraic_module(
        module_name="ascorbate_alm",
        function=ascorbate_moiety,
        compounds=["MDA", "DHA"],
        derived_compounds=["ASC"],
        parameters=["Ascorbate_total"],
    )

    m.add_algebraic_module(
        module_name="glutathion_alm",
        function=glutathion_moiety,
        compounds=["GSSG"],
        derived_compounds=["GSH"],
        parameters=["Glutathion_total"],
    )

    m.update_algebraic_module(
        module_name="ps1states",
        function=ps1analytic_mehler,
        compounds=["PC", "PCred", "Fd", "Fdred", "LHC", "ps2cs"],
        derived_compounds=["A0", "A1", "A2"],
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
            "LHC",
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

    m.update_reaction(
        rate_name="vPS1",
        function=vPS1,
        stoichiometry={"PC": 1},
        modifiers=["A0", "ps2cs"],
        dynamic_variables=["A0", "ps2cs"],
        parameters=["pfd"],
    )

    m.add_reaction(
        rate_name="vFdred",
        function=vFd_red,
        stoichiometry={"Fd": -1},
        modifiers=["Fdred", "A1", "A2"],
        parameters=["kFdred", "Keq_FAFd"],
    )

    m.add_reaction(
        rate_name="vAscorbate",
        function=vAscorbate,
        stoichiometry={"H2O2": -1, "MDA": 2},
        modifiers=["ASC"],
        dynamic_variables=["ASC", "H2O2"],
        parameters=["kf1", "kr1", "kf2", "kr2", "kf3", "kf4", "kr4", "kf5", "XT"],
    )

    m.add_reaction(
        rate_name="vMDAreduct",
        function=vMDAreduct,
        stoichiometry={"NADPH": -1, "MDA": -2},
        parameters=["kcatMDAR", "KmMDAR_NADPH", "KmMDAR_MDA", "MDAR0"],
    )

    m.add_reaction(
        rate_name="vMehler",
        function=vMehler,
        stoichiometry={
            "H2O2": 1 * m.get_parameter("convf")
        },  # required to convert as rates of PSI are expressed in mmol/mol Chl
        modifiers=["A1"],
        parameters=["O2ext", "kMehler"],
    )

    m.add_reaction(
        rate_name="vGR",
        function=vGR,
        stoichiometry={"NADPH": -1, "GSSG": -1},
        dynamic_variables=["NADPH", "GSSG"],
        parameters=["kcat_GR", "GR0", "KmNADPH", "KmGSSG"],
    )

    m.add_reaction(
        rate_name="vDHAR",
        function=vDHAR,
        stoichiometry={"DHA": -1, "GSSG": 1},
        modifiers=["GSH"],
        dynamic_variables=["DHA", "GSH"],
        parameters=["kcat_DHAR", "DHAR0", "KmDHA", "K", "KmGSH"],
    )

    m.add_reaction(
        rate_name="v3ASC",
        function=v3ASC,
        stoichiometry={"MDA": -2, "DHA": 1},
        dynamic_variables=["MDA"],
        parameters=["k3"],
    )
    return m
