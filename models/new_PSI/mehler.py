from modelbase.ode import Model

# Mehler reaction
# The Mehler reaction is a side reaction of the photosynthetic electron transport chain
# that reduces oxygen, instead of Fd
# this leads to the formation of superoxide and hydrogen peroxide. (ROS)
# This serves as an alternative electron sink to the PSI


def vAscorbate(A, H, kf1, kr1, kf2, kr2, kf3, kf4, kr4, kf5, XT):
    """lumped reaction of ascorbate peroxidase
    the cycle stretched to a linear chain with
    two steps producing the MDA
    two steps releasing ASC
    and one step producing hydrogen peroxide.
    A = Ascorbate in this case!"""
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


def vMehler(PSI_red_acceptor, O2ext, kMehler):
    """
    acceptor_side of the PSI is the side of the Mehler reaction
    This reaction is lumping the reduction of O2 instead of Fd
    resulting in Superoxide, as well as the Formation of H2O2 in one reaction.
    The entire reaction is scaled by the arbitrary parameter kMehler
    """
    #! Assumes linear dependence on both reduced PSI and external O2 concentration. This could be debated as we have a stoichiometry of 2 for the PSI_red_acceptor (giving 2 electrons to O2). Because in vivo, this is sequential and not acchieved through dimerization, we assume that the rate is linear with respect to the reduced PSI acceptor.
    return kMehler * O2ext * PSI_red_acceptor

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

    # MEHLER REACTIONS

    m.add_reaction(
        rate_name="v3_Mehler",
        function=vMehler,
        stoichiometry={
            "H2O2": +1 * m.get_parameter("convf"), # required to convert as rates of PSI are expressed in mmol/mol Chl
            "P700FA": +2,
            "P700FA-": -2
        },
        dynamic_variables=["P700FA-"],
        parameters=["O2ext", "kMehler"]
    )

    m.add_reaction(
        rate_name="v4_Mehler",
        function=vMehler,
        stoichiometry={
            "H2O2": +1 * m.get_parameter("convf"), # required to convert as rates of PSI are expressed in mmol/mol Chl
            # "P700+FA": +2, # is a derived compound!
            "P700+FA-": -2
        },  
        dynamic_variables=["P700+FA-"],
        parameters=["O2ext", "kMehler"]
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
