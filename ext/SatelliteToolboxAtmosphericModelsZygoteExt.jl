module SatelliteToolboxAtmosphericModelsZygoteExt

using SatelliteToolboxAtmosphericModels

using ForwardDiff

using Zygote.ChainRulesCore

function ChainRulesCore.rrule(::typeof(AtmosphericModels._get_doy), jd::Number)

    y = AtmosphericModels._get_doy(jd)

    function _get_doy_pullback(Δ::Number)
        return (NoTangent(), Δ)
    end

    return y, _get_doy_pullback

end

function ChainRulesCore.rrule(
    ::typeof(AtmosphericModels.nrlmsise00),
    jd::JT,
    h::HT,
    ϕ_gd::PT,
    λ::LT,
    F10ₐ::FT,
    F10::FT2,
    ap::T_AP;
    flags::AtmosphericModels.Nrlmsise00Flags = AtmosphericModels.Nrlmsise00Flags(),
    include_anomalous_oxygen::Bool = true,
    P::Union{Nothing, Matrix} = nothing
) where {JT<:Number, HT<:Number, PT<:Number, LT<:Number, FT<:Number, FT2<:Number, T_AP<:Union{Number, AbstractVector}}

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

    function nrlmsise00_pullback(Δ)

        fields = fieldnames(AtmosphericModels.Nrlmsise00Output)

        jac = reduce(hcat, ForwardDiff.gradient(
            (x) -> getfield(
                AtmosphericModels.nrlmsise00(
                    x...;
                    flags = flags,
                    include_anomalous_oxygen = include_anomalous_oxygen,
                    P = P
                ),
                field
            ),
            [jd, h, ϕ_gd, λ, F10ₐ, F10, ap]
        ) for field in fields)'

        Δvec = [Δ[i] for i in 1:length(fields)]

        vjp = Δvec' * jac

        return (NoTangent(), vjp..., NoTangent(), NoTangent(), NoTangent())
    end

    return y, nrlmsise00_pullback

end

end
