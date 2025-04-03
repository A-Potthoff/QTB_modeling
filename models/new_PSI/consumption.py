from modelbase.ode import Model

from .rate_laws import mass_action_1s


def add_consumption(model: Model) -> Model:
    model.add_parameters(
        {
            "k_ex_atp": 0.2,
            "k_ex_nadph": 0.2,
        }
    )
    model.add_reaction(
        rate_name="vEX_ATP",
        function=mass_action_1s,
        stoichiometry={"ATP": -1},
        parameters=["k_ex_atp"],
    )
    model.add_reaction(
        rate_name="vEX_NADPH",
        function=mass_action_1s,
        stoichiometry={"NADPH": -1},
        parameters=["k_ex_nadph"],
    )
    return model
