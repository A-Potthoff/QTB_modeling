from .consumption import add_consumption
from .matuszynska import get_matusznyska
from .mehler import add_mehler
from .thioredoxin import add_thioredoxin

from .rate_laws import normalize_concentration

def get_model():
    m = get_matusznyska()
    m = add_mehler(m)
    m = add_consumption(m)
    m = add_thioredoxin(m)

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

    return m
