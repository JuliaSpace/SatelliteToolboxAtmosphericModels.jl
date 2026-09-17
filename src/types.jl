## Description #############################################################################
#
# Abstract types shared by the atmospheric models.
#
############################################################################################

export AbstractAtmosphericModelOutput

"""
    abstract type AbstractAtmosphericModelOutput end

Supertype of the structures returned by the atmospheric models that provide the number
density of the species in addition to the total density.

# Implementation

Every concrete subtype must define the private functions `_model_name`, returning the
short name of the model (e.g. `"JR1971"`), `_model_description`, returning the full name
of the model (e.g. `"Jacchia-Roberts 1971"`), and `_show_fields`, returning a tuple with the
description of each field shown by `show`. Those definitions provide the shared `show`
methods.
"""
abstract type AbstractAtmosphericModelOutput end
