using Test,TightBinding
using Plots
using LinearAlgebra




@testset "Real to k" begin
    @testset "1D" begin
        @time energies = TightBinding.test_1D()
        sumt = 26.21913117088203
        @test round(sum(abs.(energies))*100) == round(sumt*100)
    end
    @testset "2Dsquare" begin
        sumt = 52.438262341764066
        @test round(sum(abs.(TightBinding.test_2Dsquare()))*100) == round(sumt*100)
    end
    @testset "Graphene" begin
        @time hamk = TightBinding.test_2DGraphene()
        sumt = 5.975026051875501
        @test round(sum(abs.(hamk))*100) == round(sumt*100)         
#        @time test_Graphene()
    end
    @testset "Iron pnictides" begin
        @time hamk =  TightBinding.test_pnictides()
        sumt = 2.680985878019015        
        @test round(sum(abs.(hamk))*100) == round(sumt*100)   

        @time hamk = TightBinding.test_pnictides_5orbitals()
        sumt = 3.0141792528265308
        @test round(sum(abs.(hamk))*100) == round(sumt*100)   

        #@time test_Fe2()
        #@time test_Fe5()
    end

    
end



@testset "k to real" begin
    @testset "Topological insulator" begin
        @time mat_e = TightBinding.test_surface()
        sumt = 17011.010972794524
        @test round(sum(abs.(mat_e))*100) == round(sumt*100)
    end
end
