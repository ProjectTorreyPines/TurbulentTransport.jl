"""
    save(input::Union{InputTGLF, InputCGYRO}, filename::AbstractString)

Write input_tglf/input_cgyro to file in input.tglf/input.cgyro/input.qlgyro format to be read by TGLF/CGYRO
"""
function save(input::Union{InputTGLF,InputCGYRO,InputQLGYRO}, filename::AbstractString)
    open(filename, "w") do io
        for key in fieldnames(typeof(input))
            if startswith(String(key), "_")
                continue
            end
            try
                value = getfield(input, key)
                if ismissing(value)
                    continue
                elseif isa(value, Int)
                    println(io, "$(key)=$(convert(Int, value))")
                elseif isa(value, String)
                    println(io, "$(key)=$value")
                elseif isa(value, Bool)
                    println(io, "$(key)=.$value.")
                else
                    println(io, "$(key)=$(convert(Float64, value))")
                end
            catch e
                println("Error writing $key to input file")
                rethrow(e)
            end
        end
    end
end

"""
    compare_two_input_tglfs(itp_1::InputTGLF, itp_2::InputTGLF)

Compares two input_tglfs, prints the difference and stores the difference in a new InputTGLF
"""
function compare_two_input_tglfs(itp_1::InputTGLF, itp_2::InputTGLF)
    itp_diff = InputTGLF()
    for field in fieldnames(typeof(itp_diff))
        if typeof(getproperty(itp_1, field)) <: String
            setproperty!(itp_diff, field, getproperty(itp_1, field) * "  " * getproperty(itp_2, field))
        else
            setproperty!(itp_diff, field, getproperty(itp_1, field) - getproperty(itp_2, field))
        end
    end

    for key in fieldnames(typeof(itp_diff))
        itp_1_value = getproperty(itp_1, key)
        itp_diff_value = getproperty(itp_diff, key)
        if typeof(itp_diff_value) <: String || typeof(itp_diff_value) <: Missing
            continue
        end
        println("Difference for $key = $(round(itp_diff_value/itp_1_value*100,digits=2)) % itp_2 = $(itp_diff_value+itp_1_value), itp_1 = $itp_1_value ")
    end

    return itp_diff
end

function diff(A::InputTGLF, B::InputTGLF)
    differences = Symbol[]
    for field in fieldnames(typeof(A))
        if getfield(A, field) === getfield(B, field)
            # pass
        else
            push!(differences, field)
        end
    end
    return differences
end

function scan(input_tglf::InputTGLF; kw...)
    # Convert keyword arguments to a dictionary for easier handling
    kw_dict = Dict(kw)

    # Base case: if no keywords are left, return the current input_tglf
    if isempty(kw_dict)
        return [deepcopy(input_tglf)]
    end

    # Extract one keyword and its values
    key, values = pop!(kw_dict)
    intglfs = InputTGLF[]

    # Iterate over the values for the current keyword
    for v in values
        tmp = deepcopy(input_tglf)
        setproperty!(tmp, key, v)

        # Recursively call scan for the rest of the keywords
        intglfs = vcat(intglfs, scan(tmp; kw_dict...))
    end

    return intglfs
end

export compare_two_input_tglfs

"""
    parse_out_tglf_gbflux(lines::String; outnames=["Gam/Gam_GB", "Q/Q_GB", "Pi/Pi_GB", "S/S_GB"])

parse out.tglf.gbflux file into a dictionary with possibility of using custom names for outputs
"""
function parse_out_tglf_gbflux(lines::String; outnames::NTuple{4,String}=("Gam/Gam_GB", "Q/Q_GB", "Pi/Pi_GB", "S/S_GB"))
    vals = map(x -> parse(Float64, x), split.(lines))
    ns = Int(length(vals) / length(outnames))
    out = Dict()
    let k = 1
        for t in outnames
            out[t*"_elec"] = vals[k]
            k += 1
            out[t*"_all_ions"] = Float64[]
            for s in 1:ns-1
                out[t*"_ion$s"] = vals[k]
                push!(out[t*"_all_ions"], vals[k])
                k += 1
            end
            out[t*"_ions"] = sum(out[t*"_all_ions"])
        end
    end
    return out
end


"""
    convertInputFromCGYRO(InputCGYRO::InputTJLF)

convert input file from InputCGYRO to InputTJLF
#"""
function convert_input_CGYRO_to_TJLF(inputCGYRO::InputCGYRO, nky::Integer)
    #get number of species
    ns = inputCGYRO.N_SPECIES
    #generate TJLF structure
    inputTJLF = InputTJLF{Float64}(ns, nky)

    inputTJLF.NS=ns
    inputTJLF.RMIN_LOC = inputCGYRO.RMIN
    inputTJLF.RMAJ_LOC = inputCGYRO.RMAJ
    inputTJLF.Q_LOC = abs(inputCGYRO.Q)
    inputTJLF.Q_PRIME_LOC=inputCGYRO.S*(inputTJLF.Q_LOC/inputTJLF.RMIN_LOC)^2
    
    inputTJLF.VEXB_SHEAR=inputCGYRO.GAMMA_E
    inputTJLF.BETAE=inputCGYRO.BETAE_UNIT
    inputTJLF.XNUE=inputCGYRO.NU_EE
    
    pressure_sum=0.0
    pressure_grad_sum=0.0
    for i in 1:ns
        n=getproperty(inputCGYRO, Symbol("DENS_" * string(i)))
        T=getproperty(inputCGYRO, Symbol("TEMP_" * string(i)))
        dn=getproperty(inputCGYRO, Symbol("DLNNDR_" * string(i)))
        dT=getproperty(inputCGYRO, Symbol("DLNTDR_" * string(i)))
        pressure = n*T
        pressure_sum += pressure 
        pressure_grad_sum += dn*T+dT*n
    end
    #should calculate it properly, but BETAE_UNIT and beta_star are very close with BETA_STAR_SCALE=1
    beta_star=inputCGYRO.BETAE_UNIT# *pressure_grad_sum/pressure_sum
    
    inputTJLF.P_PRIME_LOC=(-beta_star/8/pi)*abs(inputTJLF.Q_LOC/inputTJLF.RMIN_LOC)
    #inputTJLF.DEBYE=inputCGYRO.LAMBDA_STAR
    inputTJLF.ZMAJ_LOC=inputCGYRO.ZMAG
    inputTJLF.DZMAJDX_LOC=inputCGYRO.DZMAG
    inputTJLF.DRMAJDX_LOC=inputCGYRO.SHIFT
    inputTJLF.KAPPA_LOC=inputCGYRO.KAPPA
    inputTJLF.S_KAPPA_LOC=inputCGYRO.S_KAPPA
    inputTJLF.DELTA_LOC=inputCGYRO.DELTA
    inputTJLF.S_DELTA_LOC=inputCGYRO.S_DELTA
    inputTJLF.ZETA_LOC=inputCGYRO.ZETA
    inputTJLF.S_ZETA_LOC=inputCGYRO.S_ZETA

    
    
    for i in 1:ns
       inputTJLF.VPAR[i]=abs(inputCGYRO.MACH)
      
    end

    zeff = 0.0

    # do not use electrons to calculate zeff
    for i in 1:ns-1
        zeff += (getproperty(inputCGYRO, Symbol("Z_" * string(i)))^2* getproperty(inputCGYRO, Symbol("DENS_" * string(i))))/getproperty(inputCGYRO, Symbol("DENS_" * string(ns)))      
      
    end
     
    inputTJLF.ZEFF=zeff    
    geom_array=["SHAPE_SIN", "SHAPE_S_SIN"]
    for i in 3:6
          for prefix in geom_array
           x = Symbol(prefix*string(i))
           if x in fieldnames(typeof(inputCGYRO))
            cgyro_geom = getproperty(inputCGYRO, x)
            setproperty!(inputTJLF, Symbol(prefix*string(i)), cgyro_geom)
           end
        end
    end

    
    geom_array=["SHAPE_COS", "SHAPE_S_COS"]
    for i in 0:6
          for prefix in geom_array
           x = Symbol(prefix*string(i))
           if x in fieldnames(typeof(inputCGYRO))
            cgyro_geom = getproperty(inputCGYRO, x)
            setproperty!(inputTJLF, Symbol(prefix*string(i)), cgyro_geom)
           end
        end
    end

    plasma_array=["Z_","MASS_","DENS_","TEMP_","DLNNDR_","DLNTDR_"]
    tglf_plasma_array=["ZS","MASS","AS","TAUS","RLNS","RLTS"]
    #tglf starts with electrons, while in CGYRO it is the last species
    for i in 2:ns
        for (k,prefix) in enumerate(plasma_array)
                
         cgyro_plasma = getproperty(inputCGYRO, Symbol(prefix*string(i-1)))
         cgyro_plasma_elec=getproperty(inputCGYRO, Symbol(prefix*string(ns)))
         tglf_dic=getproperty(inputTJLF, Symbol(tglf_plasma_array[k]))
         tglf_dic[i]=cgyro_plasma
         tglf_dic[1]=cgyro_plasma_elec
        end
      end
  
    
    inputTJLF.SIGN_BT=inputCGYRO.BTCCW
    inputTJLF.SIGN_IT=inputCGYRO.IPCCW


    return inputTJLF
   
end

"""
    convert InputTJLF into InputCGYRO file. Create the same structure as InputCGYRO function working with dd

convert input file from InputCGYRO to InputTJLF
#"""

function convert_input_TJLF_to_CGYRO(inputTJLF::InputTJLF)
    #get number of species
    ns=inputTJLF.NS
    #generate InputCGYRO structure
    inputCGYRO = InputCGYRO()

   
    inputCGYRO.RMIN = inputTJLF.RMIN_LOC
    inputCGYRO.RMAJ = inputTJLF.RMAJ_LOC 
    inputCGYRO.SHIFT = inputTJLF.DRMAJDX_LOC
    inputCGYRO.KAPPA =  inputTJLF.KAPPA_LOC
    inputCGYRO.S_KAPPA = inputTJLF.S_KAPPA_LOC
    inputCGYRO.DELTA =  inputTJLF.DELTA_LOC
    inputCGYRO.S_DELTA =  inputTJLF.S_DELTA_LOC
    inputCGYRO.ZETA = inputTJLF.ZETA_LOC
    inputCGYRO.S_ZETA = inputTJLF.S_ZETA_LOC
    inputCGYRO.ZMAG = inputTJLF.ZMAJ_LOC
    inputCGYRO.DZMAG = inputTJLF.DZMAJDX_LOC

    inputCGYRO.Q= inputTJLF.Q_LOC
    inputCGYRO.S = (inputTJLF.Q_LOC/inputTJLF.RMIN_LOC)^2/inputTJLF.Q_PRIME_LOC

    inputCGYRO.BTCCW = inputTJLF.SIGN_BT
    inputCGYRO.IPCCW = inputTJLF.SIGN_IT
    inputCGYRO.BETAE_UNIT =  inputTJLF.BETAE
    inputCGYRO.N_SPECIES = inputTJLF.NS

    inputCGYRO.GAMMA_E = inputTJLF.VEXB_SHEAR
    
    inputCGYRO.NU_EE = inputTJLF.XNUE
    inputCGYRO.MACH = inputTJLF.VPAR[1]
   
   
    geom_array=["SHAPE_SIN", "SHAPE_S_SIN"]
    for i in 3:6
          for prefix in geom_array
           x = Symbol(prefix*string(i))
           tglf_value = getproperty(inputTJLF, x)
           if (tglf_value !==missing)  &&   (x in fieldnames(typeof(inputCGYRO)))        
            setproperty!(inputCGYRO, Symbol(prefix*string(i)), tglf_value)
           end
        end
    end

    
    geom_array=["SHAPE_COS", "SHAPE_S_COS"]
    for i in 0:6
        for prefix in geom_array
            x = Symbol(prefix*string(i))
            tglf_value = getproperty(inputTJLF, x)
            if (tglf_value !==missing)   &&   (x in fieldnames(typeof(inputCGYRO)))             
             setproperty!(inputCGRYO, Symbol(prefix*string(i)), tglf_value)
            end
        end
    end

    plasma_array=["Z_","MASS_","DENS_","TEMP_","DLNNDR_","DLNTDR_"]
    tglf_plasma_array=["ZS","MASS","AS","TAUS","RLNS","RLTS"]
    #tglf starts with electrons, while in CGYRO it is the last species
    for i in 2:ns
        for (k,prefix) in enumerate(plasma_array)
         tglf_dic=getproperty(inputTJLF, Symbol(tglf_plasma_array[k])) 
           
         setproperty!(inputCGYRO, Symbol(prefix*string(i-1)), tglf_dic[i])
         setproperty!(inputCGYRO, Symbol(prefix*string(ns)), tglf_dic[1])         
         
        end
      end 
      

    return inputCGYRO
   
end
