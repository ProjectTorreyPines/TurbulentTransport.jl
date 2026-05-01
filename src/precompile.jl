using PrecompileTools: @compile_workload

@compile_workload begin
    sample_input = joinpath(pkgdir(TurbulentTransport), "test", "data", "sample_input.tglf")
    input_tglf = load(InputTGLF(), sample_input)
    apply_presets!(input_tglf)
    run_tglfnn(input_tglf; model_filename="sat3_em_d3d_azf-1", warn_nn_train_bounds=false)
end
