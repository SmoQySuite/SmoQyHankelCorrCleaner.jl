module SmoQyHankelCorrCleaner

using LinearAlgebra
using Statistics
import JuMP
import HiGHS
import OnlineStats

include("init_hankel_matrix.jl")

include("projections.jl")

include("dykstra.jl")

include("hankel_correlation_cleaner.jl")
export hankel_correlation_cleaner, hankel_correlation_cleaner!

include("jackknife_hankel_correlation_cleaner.jl")
export jackknife_hankel_correlation_cleaner, jackknife_hankel_correlation_cleaner!

include("bootstrap_hankel_correlation_cleaner.jl")
export bootstrap_hankel_correlation_cleaner, bootstrap_hankel_correlation_cleaner!

end
