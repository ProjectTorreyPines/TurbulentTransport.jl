@testset "InputTGLFs Collection" begin
    @testset "construction" begin
        input1 = load_sample_input()
        input2 = load_sample_input()
        input2.BETAE = 0.005

        collection = TurbulentTransport.InputTGLFs([input1, input2])

        @test collection isa TurbulentTransport.InputTGLFs
        @test length(collection.tglfs) == 2
    end

    @testset "getindex" begin
        input1 = load_sample_input()
        input2 = load_sample_input()
        input2.BETAE = 0.005

        collection = TurbulentTransport.InputTGLFs([input1, input2])

        @test collection[1] === input1
        @test collection[2] === input2
        @test collection[1].BETAE != collection[2].BETAE
    end

    @testset "getproperty returns vector" begin
        input1 = load_sample_input()
        input2 = load_sample_input()
        input2.BETAE = 0.005
        input2.Q_LOC = 3.0

        collection = TurbulentTransport.InputTGLFs([input1, input2])

        # Getting a property should return a vector of all values
        betae_values = collection.BETAE
        @test betae_values isa Vector
        @test length(betae_values) == 2
        @test betae_values[1] == input1.BETAE
        @test betae_values[2] == input2.BETAE

        q_values = collection.Q_LOC
        @test q_values[1] == input1.Q_LOC
        @test q_values[2] == 3.0
    end

    @testset "setproperty! with vector" begin
        input1 = load_sample_input()
        input2 = load_sample_input()

        collection = TurbulentTransport.InputTGLFs([input1, input2])

        # Set different values for each element
        collection.BETAE = [0.001, 0.002]

        @test collection[1].BETAE == 0.001
        @test collection[2].BETAE == 0.002
    end

    @testset "setproperty! with scalar (broadcast)" begin
        input1 = load_sample_input()
        input2 = load_sample_input()

        collection = TurbulentTransport.InputTGLFs([input1, input2])

        # Set same value for all elements
        collection.BETAE = 0.01

        @test collection[1].BETAE == 0.01
        @test collection[2].BETAE == 0.01
    end

    @testset "setproperty! length mismatch error" begin
        input1 = load_sample_input()
        input2 = load_sample_input()

        collection = TurbulentTransport.InputTGLFs([input1, input2])

        # Setting with wrong length should error
        @test_throws AssertionError collection.BETAE = [0.001, 0.002, 0.003]
    end
end
