include("../processes/Urca.jl")
using .Urca
using Profile, ProfileView

rate_wrapper(960, 7e16, 0.1)
VSCodeServer.@profview rate_wrapper(960, 6e16, 0.1)