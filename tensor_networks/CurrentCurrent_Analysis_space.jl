using ITensors
using Plots
using LaTeXStrings
using ITensorTDVP
using DelimitedFiles
using Printf
using QuadGK
using SpecialFunctions
using Base.MathConstants
using CSV
using DataFrames


function Schwinger_Hamiltonian(;N::Int64,m::Float64,g::Float64,a::Float64)
        #=
        Construct the Schwinger Hamiltonian as a sum of Pauli strings (via staggered fermions --> Jordan-Wigner)

        Args:
        N: number of lattice points
        m: fermion mass
        g: coupling constant
        a: lattice spacing

        Yield:
        hamil: Schwinger Hamiltonian as Pauli operator
        =#


        hamil = OpSum()
        for n in 1:N-1

            # mass term
            hamil .+= m*(-1)^n , "Sz", n ## note Sz convention: Sz = sigma_z / 2

            # kinetic term
            hamil .+= 1/(2*a), "S+",n,"S-",n+1
            hamil .+= 1/(2*a), "S-",n,"S+",n+1

            # gauge term
            gaugecoef= g^2*a/8
            for i=1:n, j=1:n
                hamil .+= gaugecoef*4, "Sz",i,"Sz",j
                hamil .+= gaugecoef*((-1)^i)*2, "Sz", j
                hamil .+= gaugecoef*((-1)^j)*2, "Sz", i
            end
        end

        hamil .+= m*(-1)^N , "Sz", N

        return hamil


end


function State_charge(psi)
    #=
    Calculate charge of the input state
    =#
    list = expect(psi, "Sz")
    totalcharge=0
    for n in 1:length(list)
        totalcharge+= list[n]/2. + (-1)^n/2.
    end
    return totalcharge
end



function initialiseState(N::Int64, sites)
    #=
    Construct the initial state for the lattice model
    \chi = up for even sites, down for odd sites

    Args:
    N: number of lattice points
    sites: array of iTensor index objects which have certain properties, e.g. spin=1/2 and conserve quantum numbers

    Yield:
    initialPsi: initial state
    =#
    state = [isodd(n) ? "Up" : "Dn" for n=1:N]
    initialPsi = MPS(sites,state)
    println("\nFlux of initial state: $(flux(initialPsi))")
    println("Charge of initial state: $(State_charge(initialPsi))\n")
    return initialPsi
end


function getGroundState(pauliOp, initPsi ; N::Int64, nsweeps::Int64, maxdim, tolerance,sites)
    #=
    Convert the given Pauli operator (operator sum) into a matrix product operator (MPO) tensor network
    Produce optimised ground state and energy

    Args:
    pauliOp: Hamiltonian as Pauli operator/ operator sum
    initPsi: initial guess for ground state
    N: number of lattice points
    nsweeps: number of DMRG sweeps
    maxdim: limits on maximum MPS bond dimension allowed for each sweep
    tolerance: truncation error threshold of each sweep
    sites: array of iTensor index objects which have certain properties, e.g. spin=1/2 and conserve quantum numbers

    Yield:
    Ham: Hamiltonian as MPO
    energy: ground state energy
    psi: ground state
    =#
    Ham = MPO(pauliOp, sites)
    energy, psi = dmrg(Ham, initPsi; nsweeps, maxdim, tolerance)

    return Ham, energy, psi
end

function timeEvolution(Ham, initialPsi ; numSweeps, maxDim, tolerance ,timeT, sites, backend=0)
    #=
    Perform time evolution of an initial state

    Args:
    Ham: Hamiltonian as MPO object
    initPsi: initial state
    numSweeps: number of tdvp sweeps
    maxDim: limits on maximum MPS bond dimension allowed for each sweep
    tolerance: truncation error threshold of each sweep
    timeT: evolution time
    sites: array of iTensor index objects which have certain properties, e.g. spin=1/2 and conserve quantum numbers

    Yield:
    psiT: time-evolved state
    =#

    println("Beginning time evolution of the input state.")
    #println(H)
    #println(psi0)
    if backend==0
        psiT = tdvp(
        Ham,
        (-im)*timeT/numSweeps,
        initialPsi;
        nsweeps=numSweeps,
        reverse_step=true,
        normalize=true,
        maxdim=maxDim,
        cutoff=tolerance,
        outputlevel=1)
    elseif backend==1
        psiT = tdvp(
        Ham,
        (-im)*timeT,
        initialPsi;
        nsweeps=numSweeps,
        reverse_step=true,
        normalize=true,
        maxdim=maxDim,
        cutoff=tolerance,
        outputlevel=0,
        solver_backend="exponentiate")
    end

    return psiT
end


function createCondCondOp(n)
    #=
    Construct the \bar{\psi}\psi operator
    To be used for calculating  <\bar{\psi}(x)\psi(x) \bar{\psi}(0)\psi(0)>

    Args:
    n: (staggered) lattice site
    =#

    operator = OpSum()
    operator .+= 1, "Sz",n
    operator .+= 1, "Sz",n-1

    return operator
end


function Initialisation(N,m,g,a)
    #=
    Initialise parameters and create the Schwinger Hamiltonian

    Args:
    N: lattice points
    m: fermion mass
    g: gauge coupling
    a: lattice spacing

    Yield:
    hamil: Schwinger Hamiltonian as OpSum object
    =#

    # 0. Initialise parameters
    println("Initialising parameters: ")
    println("\t$N lattice points, $m fermion mass, $g gauge coupling, $a staggered lattice spacing\n")

    # 1. Create Schwinger model
    println("Creating Schwinger Hamiltonian...")
    hamil = Schwinger_Hamiltonian(;N,m,g,a)

    return hamil
end


function createInitialState(N)
    # Create initial state with N sites
    println("Initialising state...")
    # global sites = siteinds("S=1/2", N; conserve_qns=true) # make an array of iTensor index objects which have spin=1/2 and conserve quantum numbers
    global sites = siteinds("S=1/2", N; conserve_qns=false)
    global initPsi = initialiseState(N, sites)

    return initPsi, sites
end


function GroundState(hamil, initPsi, N, nsweeps, maxdim, cutoff, sites)
    #=
    Find ground state

    Yield:
    MPO_Ham: Hamiltonian as MPO object
    groundEnergy: energy of ground state
    groundPsi: ground state of Hamiltonian
    =#
    println("Calculating ground state...")
    global MPO_Ham, groundEnergy, groundPsi = getGroundState(hamil, initPsi ; N, nsweeps, maxdim, tolerance=cutoff,sites)
    # println("Ground state of the system is calculated as: $(groundPsi)")
    return MPO_Ham, groundEnergy, groundPsi
end


function createJOperator(omega,n,a)
    #=
    Construct the J^\omega (n) operator
    To be used for calculating  <J(x) J(0)>

    Args:
    omega: tensor index
    n: (real, even) lattice site
    N: number of sites
    =#

    operator = OpSum()
    if omega==0
        if iseven(n)
            operator += 1/a, "Sz",n
            operator += 1/a, "Sz",n-1


            # operator += 1/a, "S+",n,"S-",n
            # operator += 1/a, "S+",n+1,"S-",n+1
            # operator += 1/a, "S+",n-1,"S-",n-1
            # operator += 1/a, "S+",n,"S-",n
        else
            println("Error odd indexation")
            exit()
        end
    elseif omega==1
        if iseven(n)
            im_unit = 1.0im
            operator -= im_unit/a, "S+",n,"S-",n-1
            operator += im_unit/a, "S+",n-1,"S-",n
            # operator -= im_unit/a, "S+",n-1,"S-",n
            # operator += im_unit/a, "S+",n,"S-",n-1
        else
            println("Error odd indexation")
            exit()
        end
    end

    return operator
end



function getEvolvedStates(ind,alpha,beta,a, n,m ,MPO_Ham, groundPsi, sites)
    #=
    index:
           0 for primary correlator of interest <\Omega|J^\mu (x) J^\alpha (x) J^\beta (0) J^\nu (0)|\Omega> #Not Implemented
           1 for <J(x)J(0)> with x as pure space
           2 for <J^2(x)J(0)> with x as pure space  #Not Implemented
           3 for <J(x)J(0)> with x- as light cone coord #Not Implemented
           4 for <J^2(x)J(0)> with x- as light cone coord #Not Implemented

           5 for <condensate(x) condensate(0)> with fixed t
    =#

    n = 2*n
    m = 2*m
    println("inside computation")
    println("n= ",n, " n0= ",m)
    fixedCutOff = 1E-8

    Psi_left_op = createJOperator(alpha,n,a)
    MPO_Psi_left = MPO(Psi_left_op, sites)
    Psi_left = apply(MPO_Psi_left, groundPsi; fixedCutOff)
    # Psi_left = normalize!(Psi_left)

    Psi_right_op = createJOperator(beta,m,a)
    MPO_Psi_right = MPO(Psi_right_op, sites)
    Psi_right = apply(MPO_Psi_right, groundPsi; fixedCutOff)
    # Psi_right = normalize!(Psi_right)

    return Psi_left, Psi_right
end

function Correlator_x(ind,N,MPO_Ham, groundPsi, asv; x0, maxdim, cutoff, sites)
    #=
    Construct the space-separated, i.e. fixed time, correlator <J^\alpha (x) J^\beta (0)>
    ind=1 for now
    =#
    x_array = collect(1:N ÷ 2)
    # print(x_array)
    correlator00Vector = zeros(Complex{Float64}, length(x_array)-x0)
    correlator01Vector = zeros(Complex{Float64}, length(x_array)-x0)
    correlator10Vector = zeros(Complex{Float64}, length(x_array)-x0)
    correlator11Vector = zeros(Complex{Float64}, length(x_array)-x0)
    points_array = zeros(Float64, length(x_array)-x0)



    x_Index=0
    for nx in x_array
        if (nx-x0>0)
            x_Index+=1
            println("Running nx: ", nx, " and nx0: ", x0, " and delta x: ", nx-x0)
            println("x index: ", x_Index)
            Psi_left00, Psi_right00= getEvolvedStates(ind,0,0,asv,nx,x0, MPO_Ham, groundPsi, sites)
            Psi_left01, Psi_right01= getEvolvedStates(ind,0,1,asv,nx,x0, MPO_Ham, groundPsi, sites)
            Psi_left10, Psi_right10= getEvolvedStates(ind,1,0,asv,nx,x0, MPO_Ham, groundPsi, sites)
            Psi_left11, Psi_right11= getEvolvedStates(ind,1,1,asv,nx,x0, MPO_Ham, groundPsi, sites)


            correlator00Vector[x_Index] = inner(Psi_left00, Psi_right00)
            correlator01Vector[x_Index] = inner(Psi_left01, Psi_right01)
            correlator10Vector[x_Index] = inner(Psi_left10, Psi_right10)
            correlator11Vector[x_Index] = inner(Psi_left11, Psi_right11)
            points_array[x_Index] = nx-x0
        end
    end

    # @show LRarray

    return points_array, correlator00Vector,correlator01Vector,correlator10Vector,correlator11Vector
end


function main(ind,N,m,g,a,x0=5)
    # 1. Create Hamiltonian
    hamil = Initialisation(N,m,g,a)

    # 2. Create initial state
    initPsi, sites = createInitialState(N)

    # 3. Find ground state
    nsweeps = 30 # number of DMRG sweeps
    maxdim = [10,20,100,150,200] # limits on maximum MPS bond dimension allowed for each sweep
    cutoff = 1E-10 # truncation threshold of each sweep
    MPO_Ham, groundEnergy, groundPsi = GroundState(hamil, initPsi, N, nsweeps, maxdim, cutoff, sites)


    ####################################################################################################################
    asv=a # for compatibility with Joao's code
    Nv = N # for compatibility with Joao's code
    nsweeps = 5
    maxdim = 250
    cutoff = 1E-8
    x_array= collect(x0:2:N-x0+1)
    correlator00Vector = Vector{Complex{Float64}}(undef, length(x_array))
    correlator01Vector = Vector{Complex{Float64}}(undef, length(x_array))
    correlator10Vector = Vector{Complex{Float64}}(undef, length(x_array))
    correlator11Vector = Vector{Complex{Float64}}(undef, length(x_array))

    lattice_vals, correlator00Vector,correlator01Vector,correlator10Vector,correlator11Vector = Correlator_x(ind,N,MPO_Ham, groundPsi, asv; x0, maxdim, cutoff, sites)
    println("Real separation values are:", lattice_vals)

    plot(lattice_vals, real(correlator00Vector),color="red",label=L"Re(C^{00}(x))")
    plot!(lattice_vals, real(correlator01Vector),color="blue",label=L"Re(C^{01}(x))")
    plot!(lattice_vals, real(correlator10Vector),color="green",label=L"Re(C^{10}(x))")
    plot!(lattice_vals, real(correlator11Vector),color="purple",label=L"Re(C^{11}(x))")

    # log scale on y-axis
    # plot!(yaxis=:log10)

    xlabel!("Lattice Site n")
    ylabel!(L"\mathcal{C}^{\alpha \beta}(x)")
    title_str = "CurrCurr_Real_CorrX_N_$(N)_m_$(m)_g_$(g)_a_$(asv)"
    filename_str = "Plots/Correlators/CurrentCurrent/ind$(ind)/Julia_plots/$title_str.pdf"
    title!(title_str)
    savefig(filename_str)


    plot(lattice_vals, imag(correlator00Vector),color="red",label=L"Im(C^{00}(x))")
    plot!(lattice_vals, imag(correlator01Vector),color="blue",label=L"Im(C^{01}(x))")
    plot!(lattice_vals, imag(correlator10Vector),color="green",label=L"Im(C^{10}(x))")
    plot!(lattice_vals, imag(correlator11Vector),color="purple",label=L"Im(C^{11}(x))")
    xlabel!("Lattice Site n")
    # ylabel!(L"\langle \bar\psi(x) \psi(x)\bar\psi(0)\psi(0)\rangle / \langle \bar \psi \psi \rangle^2 ")
    ylabel!(L"\mathcal{C}^{\alpha \beta}(x)")
    title_str = "CurrCurr_Imag_CorrX_N_$(N)_m_$(m)_g_$(g)_a_$(asv)"
    filename_str = "Plots/Correlators/CurrentCurrent/ind$(ind)/Julia_plots/$title_str.pdf"
    title!(title_str)
    savefig(filename_str)

    data_dict = Dict("Lattice_sites" => lattice_vals,
                 "correlator00Vector_imag" => imag(correlator00Vector),
                 "correlator01Vector_imag" => imag(correlator01Vector),
                 "correlator10Vector_imag" => imag(correlator10Vector),
                 "correlator11Vector_imag" => imag(correlator11Vector),
                 "correlator00Vector_real" => real(correlator00Vector),
                 "correlator01Vector_real" => real(correlator01Vector),
                 "correlator10Vector_real" => real(correlator10Vector),
                 "correlator11Vector_real" => real(correlator11Vector))
    df = DataFrame(data_dict)
    CSV.write("Plots/Correlators/CurrentCurrent/ind$(ind)/CSV/CurrCurr_CorrX_N_$(N)_m_$(m)_g_$(g)_a_$(asv).csv", df)

    return
end
##################################################################################################################################################################################################################################################

# N = 100
# x0 = 20


N = 20
x0=5

# g = 0.6
g = 0.5

a = 1.0


mPhys = 0.25
m=mPhys-g*g*a/8

ind=1
main(ind,N,m,g,a,x0)
