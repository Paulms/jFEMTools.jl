struct TrialFunction{T<:AbstractElement}
    element::T
    components::Int
end

TrialFunction(el::AbstractElement) = TrialFunction(el,1)

function getnlocaldofs(u::TrialFunction, cell::Cell)
    getnlocaldofs(u.element, cell)*u.components
end

getncomponents(u::TrialFunction) = u.components
