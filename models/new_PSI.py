from modelbase.ode import Model

from .rate_laws import mass_action_22_rev
from .rate_laws import moiety_3
from .rate_laws import vPS1
from .rate_laws import normalize_concentration

# as written in matuszynska.py:
# m.add_algebraic_module(
#     module_name="ps1states",
#     function=ps1states,
#     compounds=["PC", "PCred", "Fd", "Fdred", "LHC", "ps2cs"],
#     derived_compounds=["A1"],
#     parameters=["PSItot", "kFdred", "Keq_FAFd", "Keq_PCP700", "kPCox", "pfd"],
# )

# commented out in matuszynska.py:
    # m.add_reaction( # ! attention
#     rate_name="vPS1",
#     function=vPS1,
#     stoichiometry={"Fd": -1, "PC": 1},
#     modifiers=["A1", "ps2cs"],
#     dynamic_variables=["A1", "ps2cs"],  # doesn't depend on Fd
#     parameters=["pfd"],
# )

# def vFd_red(Fd, Fdred, A1, A2, kFdred, Keq_FAFd): 
#     """rate of the redcution of Fd by the activity of PSI
#     A1 -> A2 (old model) or
#     A1 -> A3 / A2 -> A0 (new model)
#     used to be equall to the rate of PSI but now
#     alternative electron pathway from Fd allows for the production of ROS
#     hence this rate has to be separate
#     """
#     return kFdred * Fd * A1 - kFdred / Keq_FAFd * Fdred * A2


    # m.add_reaction(
    #     rate_name="vFdred",
    #     function=vFd_red,
    #     stoichiometry={"Fd": -1},
    #     modifiers=["Fdred", "A1", "A2"],
    #     parameters=["kFdred", "Keq_FAFd"],
    # )
    
def add_PSI(m) -> Model:
    m.add_compounds(["P700FA", "P700+FA-", "P700FA-"])

    m.add_algebraic_module(
        module_name="P700+FA_alm",
        function=moiety_3,
        compounds=["P700FA-", "P700FA", "P700+FA-"],
        derived_compounds=["P700+FA"],
        parameters=["PSItot"]
    )

    m.update_reaction(
        rate_name="vPS1",
        function=vPS1,
        stoichiometry={"P700FA": -1, "P700+FA-": 1},
        #modifiers=["P700FA", "ps2cs"],                      # ! try if commenting out changes anything
        dynamic_variables=["P700FA", "ps2cs"],              # ?
        parameters=["pfd"],
    )

    # m.add_reaction_from_args(
    #     rate_name="v2_to_P700FA-",
    #     function=mass_action_22_rev,
    #     stoichiometry={"P700+FA-": -1, "P700FA-": +1, "PC": +1}, # "PCred": -1 not included because computed by moiety
    #     # modifiers=["PCred"], # has to be included here then!
    #     # parameters=["kPCox", "Keq_PCP700"],
    #     args=["P700+FA-", "PCred", "PC", "P700FA-", "kPCox", "Keq_PCP700"]
    # )

    m.add_reaction(
        rate_name="v2_to_P700FA-",
        function=mass_action_22_rev,
        stoichiometry={"P700+FA-": -1, "PC": +1, "P700FA-": +1}, # "PCred": -1 not included because computed by moiety
        dynamic_variables=["P700+FA-", "PCred", "PC", "P700FA-"],
        parameters=["kPCox", "Keq_PCP700"],
        # modifiers=["PCred"], # has to be included here then!
        # args=["P700+FA-", "PCred", "PC", "P700FA-", "kPCox", "Keq_PCP700"]
    )
# UserWarning: Supplied dynamic variables {'P700FA-', 'PC'} for rate v2_to_P700FA- that aren't in substrates or modifiers

    m.add_reaction_from_args(
        rate_name="v3_to_P700FA",
        function=mass_action_22_rev,
        stoichiometry={"P700FA-": -1, "Fd": -1, "P700FA": +1},
        # modifiers=["Fdred"],  
        # parameters=["kFdred", "Keq_FAFd"],
        args=["P700FA-", "Fd", "P700FA", "Fdred", "kFdred", "Keq_FAFd"]
    )

    m.add_reaction_from_args(
        rate_name = "v4_to_P700+FA",
        function = mass_action_22_rev,
        stoichiometry = {"P700+FA-": -1, "Fd": -1},
        # modifiers = ["Fdred", "P700+FA"],
        # parameters = ["kFdred", "Keq_FAFd"],
        args = ["P700+FA-", "Fd", "P700+FA", "Fdred", "kFdred", "Keq_FAFd"]
    )

    m.add_reaction_from_args(
        rate_name = "v5_to_P700FA",
        function = mass_action_22_rev,
        stoichiometry = {"P700FA": +1, "PC": +1},
        # modifiers = ["PCred", "P700+FA"],
        # parameters=["kPCox", "Keq_PCP700"],
        args = ["P700+FA", "PCred", "P700FA", "PC", "kPCox", "Keq_PCP700"]
    )

    # normalization of PSI compounds

    m.add_algebraic_module(
        module_name="rel_P700+FA_alm",
        function=normalize_concentration,
        compounds=["P700+FA"],
        derived_compounds=["rel_P700+FA"],
        parameters=["PSItot"],
        args = ["P700+FA", "PSItot"]
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
        module_name="rel_P700FA-_alm",
        function=normalize_concentration,
        compounds=["P700FA-"],
        derived_compounds=["rel_P700FA-"],
        parameters=["PSItot"],
        args = ["P700FA-", "PSItot"]
    )

    m.add_algebraic_module(
        module_name="rel_P700+FA-_alm",
        function=normalize_concentration,
        compounds=["P700+FA-"],
        derived_compounds=["rel_P700+FA-"],
        parameters=["PSItot"],
        args = ["P700+FA-", "PSItot"]
    )

    return m




# ? how do I deal with A3 not being computed explicitly? Do I remove it from the stoichiometry?
# ? can I reuse kPCox and kFdred?
# ? can I reuse Keq_PCP700 and Keq_FAFd?
# ? How do I define vPS1?
# ? what is ps2cs?


