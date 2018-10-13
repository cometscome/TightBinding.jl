using Test,TightBinding

sumt = 52.438262341764066
@test round(sum(abs.(TightBinding.test_2Dsquare()))) == round(sumt)
