from modelbase.ode import Model

from .rate_laws import mass_action_1s, mass_action_2s, moiety_1, proportional


def add_thioredoxin(model: Model) -> Model:
    for i in ["fCBB", "V1", "V6", "V9", "V13", "Vst"]:
        model.remove_derived_parameter(i)
    model.remove_parameters(["Km_fcbb", "Vmax_fcbb"])

    model.add_compounds(["TR_ox", "E_inactive"])
    model.add_parameters(
        {
            "thioredoxin_tot": 1,
            "e_cbb_tot": 6,
            "k_fd_tr_reductase": 1,
            "k_e_cbb_activation": 1,
            "k_e_cbb_relaxation": 0.1,
        }
    )
    model.add_algebraic_module(
        module_name="thioredoxin_alm",
        function=moiety_1,
        compounds=["TR_ox"],
        derived_compounds=["TR_red"],
        parameters=["thioredoxin_tot"],
    )
    model.add_algebraic_module(
        module_name="e_cbb_alm",
        function=moiety_1,
        compounds=["E_inactive"],
        derived_compounds=["E_active"],
        parameters=["e_cbb_tot"],
    )

    for rate_name, par in {
        "vRuBisCO": "V1",
        "vFBPase": "V6",
        "v9": "V9",
        "v13": "V13",
        "vStarch": "Vst",
    }.items():
        model.add_algebraic_module(
            module_name=par,
            function=proportional,
            derived_compounds=[par],
            compounds=["E_active"],
            parameters=[f"{par}_base"],
        )
        rate = model.rates[rate_name].copy()
        modifiers = rate["modifiers"] + [par]
        dynamic_variables = rate["dynamic_variables"] + [par]
        parameters = rate["parameters"]
        parameters.remove(par)
        model.update_rate(
            rate_name=rate_name,
            substrates=rate["substrates"],
            products=rate["products"],
            modifiers=modifiers,
            parameters=parameters,
            dynamic_variables=dynamic_variables,
            args=rate["args"],
            reversible=rate["reversible"],
        )

    model.add_reaction(
        rate_name="vFdTrReductase",
        function=mass_action_2s,
        stoichiometry={"TR_ox": -1, "Fd": 1},
        modifiers=["Fdred"],
        parameters=["k_fd_tr_reductase"],
    )
    model.add_reaction(
        rate_name=f"vE_activation",
        function=mass_action_2s,
        stoichiometry={"E_inactive": -5, "TR_ox": 5},
        modifiers=["TR_red"],
        parameters=["k_e_cbb_activation"],
    )
    model.add_reaction(
        rate_name="vE_inactivation",
        function=mass_action_1s,
        stoichiometry={"E_inactive": 5},
        modifiers=["E_active"],
        parameters=["k_e_cbb_relaxation"],
    )
    return model
