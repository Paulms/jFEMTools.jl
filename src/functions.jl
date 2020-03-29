struct TrialFunction{dim,T}
    fs::AbstractDiscreteFunctionSpace{dim,T}
end

@inline getfunctionspace(u::TrialFunction) = u.fs

# Test Functions

struct TestFunction{dim,T}
    fs::AbstractDiscreteFunctionSpace{dim,T}
end

@inline getfunctionspace(u::TestFunction) = u.fs