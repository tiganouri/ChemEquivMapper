from .step1_refmet import Step1RefMet
from .step2_chebi import Step2ChebiLookup
from .step3_ro_reactome import Step3RO
from .step4_isa_specialization import Step4IsASpecialization
from .step5_isa_generalisation import Step5IsAGeneralisation

__all__ = [
    "Step1RefMet",
    "Step2ChebiLookup",
    "Step3RO",
    "Step4IsASpecialization",
    "Step5IsAGeneralisation",
]

# Optional: handy registry for pipelines
STEP_REGISTRY = {
    "step1_refmet": Step1RefMet,
    "step2_chebi": Step2ChebiLookup,
    "step3_ro": Step3RO,
    "step4_isa_specialization": Step4IsASpecialization,
    "step5_isa_generalisation": Step5IsAGeneralisation,
}
