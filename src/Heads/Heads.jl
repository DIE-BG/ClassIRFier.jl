"""
    module Heads

Tools for classifying simulated impulse response functions (IRFs).

# Description
The `Heads` module provides methods to classify IRFs based on their characteristics. These classifications are intended for preprocessing and labeling simulated IRFs, which are later used to train models in `ClassIRFier.jl`.

# Exports
- `head_sign`: Classifies an IRF according to the dominance of its positive or negative area.

# Usage
Import this module and use the provided methods to assign labels to simulated IRFs before model training.
"""
module Heads

export head_sign
include(joinpath(@__DIR__, "head_sign.jl"))

end
