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
    println("Ground state of the system is calculated as: $(groundPsi)")
    return MPO_Ham, groundEnergy, groundPsi
end

function Corr_4pFo2Pf(n) ## this lattice n, not physical n
  operator = OpSum()


  if iseven(n)
    operator .+= 1 , "Sz", n
    operator .+= -1 , "Sz", n-1
  else
    println("Error odd indexation")
    exit()
  end


  return operator
end


function createJOperator(omega,n)
        #=
        Construct the J^\omega (n) operator
        To be used for calculating  <J(x) J(0)>

        Args:
        omega: tensor index
        n: (real, even) lattice site
        =#

        operator = OpSum()
        if omega==0
            if iseven(n)
                operator += 2, "Sz",n
                operator += 2, "Sz",n+1
                operator += 2, "Sz",n-1
                operator += 2, "Sz",n
                # operator += 1, "S+",n,"S-",n
                # operator += 1, "S+",n+1,"S-",n+1
                # operator += 1, "S+",n-1,"S-",n-1
                # operator += 1, "S+",n,"S-",n
            else
                println("Error odd indexation")
                exit()
            end
        elseif omega==1
            if iseven(n)
                im_unit = 1.0im
                operator -= im_unit, "S+",n,"S-",n+1
                operator += im_unit, "S+",n+1,"S-",n
                operator -= im_unit, "S+",n-1,"S-",n
                operator += im_unit, "S+",n,"S-",n-1
            else
                println("Error odd indexation")
                exit()
            end
        end

        return operator
end

function createState(n)
    return Corr_4pFo2Pf(n)
end



function getEvolvedStates(n,m ,MPO_Ham, groundPsi, sites)

    n=2*n
    m = 2*m
    println("inside computation")
    println("n= ",n, " n0= ",m)
    fixedCutOff = 1E-8


    # Psi_left_op = createCondCondOp(n)
    Psi_left_op = createState(n)
    MPO_Psi_left = MPO(Psi_left_op, sites)
    Psi_left = apply(MPO_Psi_left, groundPsi; fixedCutOff)
    # Psi_left = normalize!(Psi_left)

    # Psi_right_op = createCondCondOp(m)
    Psi_right_op = createState(m)
    MPO_Psi_right = MPO(Psi_right_op, sites)
    Psi_right = apply(MPO_Psi_right, groundPsi; fixedCutOff)
    # Psi_right = normalize!(Psi_right)

    return Psi_left, Psi_right
end

function Correlator_x(N,MPO_Ham, groundPsi, asv; x0, maxdim, cutoff, sites)
    #=
    index:
           5 for <condensate(x) condensate(0)> with fixed t
    =#

    # x_array= collect(x0:2:N-x0+1)
    x_array = collect(1:N ÷ 2)
    # print(x_array)
    LRarray = zeros(Float64, length(x_array)-x0)
    points_array = zeros(Float64, length(x_array)-x0)

    println("x array: $x_array")

    x_Index=0
    for nx in x_array
        if (nx-x0>0)
            x_Index+=1
            println("Running nx: ", nx, " and nx0: ", x0, " and delta x: ", nx-x0)
            println("x index: ", x_Index)
            Psi_left, Psi_right= getEvolvedStates(nx,x0, MPO_Ham, groundPsi, sites)
            LRarray[x_Index] = inner(Psi_left, Psi_right)/ (2*a)^2
            points_array[x_Index] = nx-x0
        end
    end

    # @show LRarray

    return points_array, LRarray
end


##################################################################################################################################################################################################################################################
function condensate_site!(list,psi)
  zs = expect(psi, "Sz")
  # push!(list, [(-1)^n * z for (n, z) in enumerate(zs)])
  push!(list,[(-1)^n for n in 1:length(zs)].*zs)

end

##################################################################################################################################################################################################################################################
function main(N,m,g,a,x0=5)
    # 1. Create Hamiltonian
    global hamil = Initialisation(N,m,g,a)


    # 2. Create initial state
    initPsi, sites = createInitialState(N)

    # 3. Find ground state
    nsweeps = 30 # number of DMRG sweeps
    maxdim = [10,20,100,150,200] # limits on maximum MPS bond dimension allowed for each sweep
    cutoff = 1E-10 # truncation threshold of each sweep
    global MPO_Ham, groundEnergy, groundPsi = GroundState(hamil, initPsi, N, nsweeps, maxdim, cutoff, sites)


    ####################################################################################################################
    asv=a
    Nv = N
    nsweeps = 5
    maxdim = 250
    cutoff = 1E-8
    x_array= collect(x0:2:N-x0+1)
    correlatorVector = zeros(Float64, length(x_array)-x0)


    l1=[]
    sites_physical = 1:Nv/2
    condensate_site!(l1,groundPsi)

    real_x_values, correlatorVector = Correlator_x(N,MPO_Ham, groundPsi, asv; x0, maxdim, cutoff, sites)
    # real_x_values = [((x-x0)/2)+1 for x in x_array]
    println("Real separation values are:", real_x_values)

    # Analytic function for S_F (x) / <psi bar psi>^2
    fourp_over_twop_cont(x) = real(cosh(2*pi*besselk(0, g/sqrt(pi)*x*a)  ) )### here one has to adjust for size of lattice space; if asv =1 is for the computational lattice then real lattice has a=2

    ### Euclidean signature means cosh-->cos?? !!
    fourp_over_twop_cont_values = [fourp_over_twop_cont(x*2) for x in real_x_values]
    # @show fourp_over_twop_cont_values
    # psiBarpsi = -m*(g^.577)/(2*pi)
    ####################################################################################################################
    # l1_sq = [x^2 for x in l1]

    # subset_real_x_values = real_x_values[1:end-x0]
    # subset_l2 = correlatorVector[1:end-x0]
    # subset_fourp_over_twop_cont = fourp_over_twop_cont_values[1:end-x0]
    # @show length(subset_real_x_values)
    # @show length(subset_l2)
    # @show length(subset_fourp_over_twop_cont)


    con_den =  l1[1][Int(Nv/2)]
    println("Condensate= ", con_den)
    # @show real(con_den^2 * subset_fourp_over_twop_cont)

    # @show correlatorVector
    plot(real_x_values, real(correlatorVector),color="red",label="CondCond")
    plot!(real_x_values, real(con_den^2 * fourp_over_twop_cont_values), color="blue",label="Analytic Result")

    # log scale on y-axis
    plot!(yaxis=:log10)

    xlabel!("Lattice Site n")
    # ylabel!(L"\langle \bar\psi(x) \psi(x)\bar\psi(0)\psi(0)\rangle / \langle \bar \psi \psi \rangle^2 ")
    ylabel!(L"\mathcal{C}(x)")
    title_str = "log_CorrX_N_$(N)_m_$(m)_g$(g)"
    filename_str = "Plots/Correlators/CondCond/ind5/Julia_plots/$title_str.pdf"
    title!(title_str)
    savefig(filename_str)

    plot(real_x_values, correlatorVector ./ (con_den)^2 ./ fourp_over_twop_cont_values,color="red",label="ratio")
    # factor of 4 bc we have (2a)^2?

    xlabel!("Lattice Site n")
    # ylabel!(L"\langle \bar\psi(x) \psi(x)\bar\psi(0)\psi(0)\rangle / \langle \bar \psi \psi \rangle^2 ")
    ylabel!(L"\mathcal{C}(x)")
    title_str = "CondCond_CorrX_N_$(N)_m_$(m)_g$(g)"
    filename_str = "Plots/Correlators/CondCond/ind5/Julia_plots/$title_str.pdf"
    title!(title_str)
    savefig(filename_str)

    storage_values = correlatorVector./(con_den)^2. / fourp_over_twop_cont_values
    data_dict = Dict("Separation Values" => real_x_values,
                 "Ratio Values" => storage_values)
    df = DataFrame(data_dict)
    CSV.write("Plots/Correlators/CondCond/ind6/CSV/CondCond_CorrX_N_$(N)_m_$(m)_g_$(g)_a_$(asv).csv", df)

    return
end
##################################################################################################################################################################################################################################################

N = 20
# N=50

g = 0.5
# g = 1.2

a = 1.0
# a = 0.6

# x0 = 15
x0=5

mPhys = 0
# mPhys = 0.25
m=mPhys-g*g*a/8

main(N,m,g,a,x0)
