from .consumption import add_consumption
from .matuszynska import get_matusznyska
from .mehler import add_mehler
from .thioredoxin import add_thioredoxin
from .new_PSII import add_PSII

from .rate_laws import normalize_concentration
from .rate_laws import normalize_2_concentrations

def get_model():
    m = get_matusznyska()
    m = add_PSII(m)
    m = add_mehler(m)
    m = add_consumption(m)
    m = add_thioredoxin(m)

    # useful normalizations of concentrations ! 

    m.add_algebraic_module(
        module_name="pq_redoxstate",
        function=normalize_concentration,
        compounds=["PQred"],
        derived_compounds=["PQ_redoxstate"],
        parameters=["PQtot"],
    )
    
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
    
    m.add_algebraic_module(
        module_name="nadp_redoxstate",
        function=normalize_concentration,
        compounds=["NADPH"],
        derived_compounds=["NADP_redoxstate"],
        parameters=["NADPtot"],
    )

    m.add_algebraic_module(
        module_name="energystate",
        function=normalize_concentration,
        compounds=["ATP"],
        derived_compounds=["ATP_norm"],
        parameters=["APtot"],
    )

    # PSI

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
    m.add_algebraic_module(
        module_name="rel_P700_alm",
        function=normalize_2_concentrations,
        compounds=["P700FA-", "P700FA"],
        derived_compounds=["rel_P700"],
        parameters=["PSItot"]
    )
    m.add_algebraic_module(
        module_name="rel_P700+_alm",
        function=normalize_2_concentrations,
        compounds=["P700+FA-", "P700+FA"],
        derived_compounds=["rel_P700+"],
        parameters=["PSItot"]
    )

    # PSII

    m.add_algebraic_module(
        module_name="rel_B0_alm",
        function=normalize_concentration,
        compounds=["B0"],
        derived_compounds=["rel_B0"],
        parameters=["PSIItot"],
        args = ["B0", "PSIItot"]
    )
    m.add_algebraic_module(
        module_name="rel_B1_alm",
        function=normalize_concentration,
        compounds=["B1"],
        derived_compounds=["rel_B1"],
        parameters=["PSIItot"],
        args = ["B1", "PSIItot"]
    )
    m.add_algebraic_module(
        module_name="rel_B2_alm",
        function=normalize_concentration,
        compounds=["B2"],
        derived_compounds=["rel_B2"],
        parameters=["PSIItot"],
        args = ["B2", "PSIItot"]
    )
    m.add_algebraic_module(
        module_name="rel_B3_alm",
        function=normalize_concentration,
        compounds=["B3"],
        derived_compounds=["rel_B3"],
        parameters=["PSIItot"],
        args = ["B3", "PSIItot"]
    )

    return m
