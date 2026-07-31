## Description #############################################################################
#
# This file integrates the atmospheric models with Zygote.jl. The NRLMSISE-00 model
# contains mutations that Zygote.jl cannot differentiate. Hence, we define a rule that
# computes the model Jacobian using ForwardDiff.jl, which is efficient here given the small
# number of inputs. Notice that the rule for `AtmosphericModels._get_doy` is defined in the
# ChainRulesCore.jl extension, which is always loaded together with this one.
#
############################################################################################

module SatelliteToolboxAtmosphericModelsZygoteExt

using SatelliteToolboxAtmosphericModels

using ForwardDiff

using Zygote.ChainRulesCore

function ChainRulesCore.rrule(
    ::typeof(AtmosphericModels.nrlmsise00),
    jd::Number,
    h::Number,
    ϕ_gd::Number,
    λ::Number,
    F10ₐ::Number,
    F10::Number,
    ap::Number;
    flags::AtmosphericModels.Nrlmsise00Flags = AtmosphericModels.Nrlmsise00Flags(),
    include_anomalous_oxygen::Bool = true,
    P::Union{Nothing, AbstractMatrix} = nothing
)
    y = AtmosphericModels.nrlmsise00(
        jd,
        h,
        ϕ_gd,
        λ,
        F10ₐ,
        F10,
        ap;
        flags = flags,
        include_anomalous_oxygen = include_anomalous_oxygen,
        P = P
    )

    fields = fieldnames(AtmosphericModels.Nrlmsise00Output)

    function nrlmsise00_pullback(Δ)
        # Compute the Jacobian of all output fields with respect to the inputs in a single
        # forward sweep. Notice that we must not forward the user-provided matrix `P` here
        # since it cannot store the dual numbers used by ForwardDiff.jl.
        jac = ForwardDiff.jacobian(
            x -> begin
                out = AtmosphericModels.nrlmsise00(
                    x...;
                    flags = flags,
                    include_anomalous_oxygen = include_anomalous_oxygen,
                    P = nothing
                )
                collect(getfield(out, f) for f in fields)
            end,
            [jd, h, ϕ_gd, λ, F10ₐ, F10, ap]
        )

        Δvec = [Δ[i] for i in 1:length(fields)]
        vjp  = Δvec' * jac

        # The pullback must return one tangent per primal argument (the function and the
        # seven positional arguments). The keyword arguments must not receive tangents.
        return (NoTangent(), vjp...)
    end

    return y, nrlmsise00_pullback
end

end
