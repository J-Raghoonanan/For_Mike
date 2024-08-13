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
        (-im)*timeT,
        initialPsi;
        nsweeps=numSweeps,
        reverse_step=true,
        normalize=true,
        maxdim=maxDim,
        cutoff=tolerance,
        outputlevel=1,
        solver_backend="exponentiate")
    elseif backend==1
        psiT = tdvp(
        Ham,
        (-im)*timeT,
        initialPsi;
        nsweeps=numSweeps,
        reverse_step=true,
        normalize=false,
        maxdim=maxDim,
        cutoff=tolerance,
        outputlevel=1,
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
    sites = siteinds("S=1/2", N; conserve_qns=false)
    initPsi = initialiseState(N, sites)

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

function getEvolvedStates(n,m ,MPO_Ham, groundPsi,initialLeftState,initialRightState, nsweeps, maxdim, cutoff, dt, sites)

    n=2*n
    m = 2*m
    println("inside computation")
    println("n= ",n, " n0= ",m)
    fixedCutOff = 1E-8

    #=
    We can incrementally evolve Psi_right since
    Psi_right = exp(-iHt) {J|omega>}
    so we only need to apply J to ground state once, then each time evolution is just
    an additional exp(-iH dt)

    We do the same thing for the ground state with Psi_left
    =#

    if nsweeps==0
        psiT = groundPsi #starting with t=0, so no time evolution
    else
        psiT = timeEvolution(MPO_Ham, initialLeftState ; numSweeps=1, maxDim=maxdim, tolerance=cutoff ,timeT=dt, sites=sites, backend=0)
        #backend = 0 ==> normalize = true
    end

    Psi_left_op = Corr_4pFo2Pf(n)
    MPO_Psi_left = MPO(Psi_left_op, sites)
    Psi_left = apply(MPO_Psi_left, psiT; fixedCutOff)
    # Psi_left = normalize!(Psi_left)

    if nsweeps==0
        Psi_right_op = Corr_4pFo2Pf(m)
        MPO_Psi_right = MPO(Psi_right_op, sites)
        Psi_right = apply(MPO_Psi_right, groundPsi; fixedCutOff)
        # Psi_right = normalize!(Psi_right)
        # start with t=0, so no time evolution here
    else
        Psi_right = timeEvolution(MPO_Ham, initialRightState; numSweeps=1, maxDim=maxdim, tolerance=cutoff ,timeT=dt, sites=sites, backend=1)
        #backend = 1 ==> normalize = false
    end


    return psiT, Psi_left, Psi_right
end

function Correlator_t(N,MPO_Ham, groundPsi,a,dt;nx, x0, cparameter=1, nsweeps, maxdim, cutoff, sites)

    time_array = zeros(Float64, nsweeps+1)
    Result_array = zeros(ComplexF64, length(time_array))
    println("Totaltime= ", dt*nsweeps)
    println("Total timesteps= ", nsweeps)
    println("dt value= ", dt)

    t_Index=0
    prevPsiLeft = groundPsi
    prevPsiRight = groundPsi
    for timeSubInterval in time_array
        t_Index+=1
        nsweep = t_Index-1
        time_step = dt*nsweep # want to start at time=0, but Julia is 1-indexed
        time_array[t_Index] = time_step
        println("Time now:",time_array[t_Index] )
        println("Number of Trotter steps:",t_Index-1 )
        println("Running nt: $(time_step) with nx: $nx and nx0: $x0, and delta x: $(nx-x0)")
        println()
        newPsiLeft, Psi_left, Psi_right= getEvolvedStates(nx,x0, MPO_Ham, groundPsi,prevPsiLeft,prevPsiRight, nsweep, maxdim, cutoff, dt, sites)

        Result_array[t_Index] = inner(Psi_left, Psi_right)/ (2*a)^2
        @show Result_array[t_Index]*4
        prevPsiLeft = newPsiLeft
        prevPsiRight = Psi_right
    end

    # @show LRarray

    return time_array, Result_array
end



##################################################################################################################################################################################################################################################
function condensate_site!(list,psi)
  zs = expect(psi, "Sz")
  push!(list,[(-1)^n for n in 1:length(zs)].*zs)

end


function Cond_Cond_lattice(psi0,psitime,psi_r,n_phys,n0_phys, as,sites,dt;  mmaxdim, tol)

  n = 2*n_phys
  n0 = 2*n0_phys


  Operator_n = MPO(Corr_4pFo2Pf(n), sites )
  Operator_n0 = MPO(Corr_4pFo2Pf(n0), sites )

  cutoff = 1E-8

  # psitime = Exp(-iH t)|GS>, with t the t=current_t- deltat


  psi_l = tdvp( MPO_Ham, (-im)*dt,
  psitime;
  nsweeps=1,
  reverse_step=true,
  normalize=true,
  maxdim=mmaxdim,
  cutoff=tol,
  outputlevel=1,
  solver_backend="exponentiate" )

  #=
  # here we use the fact that |GS> is eigenstate of H so U|Gs> is trivial

  psi_l= psi0

  =#

  psi_l_final = apply(Operator_n, psi_l; cutoff)
  #normalize!(psi_l )


  ########################################################

  #psi_r = apply(Operator_n0, psi0; cutoff)
  #normalize!(psi_r )

  # psi_r = O |GS>

  psi_r_final = tdvp( MPO_Ham, (-im)*dt,
  psi_r;
  nsweeps=1,
  reverse_step=true,
  normalize=false,
  maxdim=mmaxdim,
  cutoff=tol,
  outputlevel=1,
  solver_backend="exponentiate" )




  return psi_l, psi_r_final , inner(psi_l_final, psi_r_final)/4 ## psi_l is next time evolved GS; psi_r_final is  Op |GS> evolved up to t

end


function Correlator_Cond_Cond(psi0,x,x0, as,sites,totaltime,Ntimesteps;Nv,  mmaxdim, tol) # Cond_Cond_lattice_new(psi0,psitime,psi_r,n_phys,n0_phys, as,sites,dt;  mmaxdim, tol)

  dt= totaltime/Ntimesteps
  cutoff = 1E-8

  time_array = zeros(Float64, Ntimesteps+1)
  Result_array = zeros(ComplexF64, Ntimesteps+1)
  println("Totaltime= ", totaltime)
  println("Total timesteps= ", Ntimesteps)
  println("dt value= ", dt)

  n0 = 2*x0
  Operator_n0 = MPO(Corr_4pFo2Pf(n0), sites )
  psi_r = apply(Operator_n0, psi0; cutoff)

  n = 2*x
  Operator_n = MPO(Corr_4pFo2Pf(n), sites )
  psi_l_init = apply(Operator_n, psi0; cutoff)

  t_index=0

  prevPsiLeft=psi0
  prevPsiRight=psi_r


  for _ in time_array
        time_array[t_index+1] = dt*t_index
        println("Time now:",time_array[t_index+1] )
        println("Number of Trotter steps:",t_index )
        if t_index==0
          Result_array[t_index+1] = inner(psi_l_init, psi_r)/4
        else
            PsiLeft, PsiRight ,Result_array[t_index+1] = Cond_Cond_lattice(psi0,prevPsiLeft,prevPsiRight,x,x0, as,sites,dt;  mmaxdim, tol)
            prevPsiLeft = PsiLeft
            prevPsiRight = PsiRight
        end
        t_index+=1
  end

  # end

  return time_array, Result_array

end

##################################################################################################################################################################################################################################################
function main(N,m,g,a,x0)
    # 1. Create Hamiltonian
    global hamil = Initialisation(N,m,g,a)


    # 2. Create initial state
    initPsi, sites = createInitialState(N)

    # 3. Find ground state
    nsweeps = 20 # number of DMRG sweeps
    maxdim = [10, 20, 100, 150, 200] # limits on maximum MPS bond dimension allowed for each sweep
    cutoff = 1E-10 # truncation threshold of each sweep
    global MPO_Ham, groundEnergy, groundPsi = GroundState(hamil, initPsi, N, nsweeps, maxdim, cutoff, sites)
    ####################################################################################################################

    asv=a
    Nv = N
    nsweeps = 150
    maxdim = 350
    cutoff = 1E-10
    totaltime = 60
    dt = totaltime / nsweeps
    init_time_array = collect(0:1:nsweeps)
    correlatorVector = Vector{Complex{Float64}}(undef, length(init_time_array))
    nx = floor(Int, x0+2)
    cparameter = 1

    l1=[]
    sites_physical = 1:Nv/2
    condensate_site!(l1,groundPsi)
    con_den =  l1[1][Int(Nv/2)]
    println("Condensate= ", con_den)

    cutoff = 1E-10
    time_array, correlatorVector = Correlator_t(N,MPO_Ham, groundPsi, asv, dt; nx, x0, cparameter, nsweeps, maxdim, cutoff, sites)
    # @show correlatorVector
    # time_array, correlatorVector = Correlator_Cond_Cond(groundPsi,nx,x0, asv,sites,totaltime,nsweeps; Nv, mmaxdim=350, tol=1E-10)
    ####################################################################################################################

    plot(time_array, real(correlatorVector),color="red",label="Real CondCond")
    plot!(time_array, imag(correlatorVector),color="blue",label="Imaginary CondCond")

    # log scale on y-axis
    # plot!(yaxis=:log10)

    xlabel!("Time")
    # ylabel!(L"\langle \bar\psi(x) \psi(x)\bar\psi(0)\psi(0)\rangle / \langle \bar \psi \psi \rangle^2 ")
    ylabel!(L"\mathcal{C}(t)")
    title_str = "CondCond_CorrT_N_$(N)_m_$(m)_g$(g)"
    filename_str = "Plots/Correlators/CondCond/ind6/Julia_plots/$title_str.pdf"
    title!(title_str)
    savefig(filename_str)

    data_dict = Dict("time_values" => time_array,
                 "correlator_imag" => imag(correlatorVector)./(con_den)^2,
                 "correlator_real" => real(correlatorVector)./(con_den)^2)
    df = DataFrame(data_dict)
    CSV.write("Plots/Correlators/CondCond/ind6/CSV/CondCond_CorrT_N_$(N)_m_$(m)_g_$(g)_a_$(asv)_time_$(totaltime).csv", df)

    return
end

# Set 5
# N = 100
N=20
mPhys = 0.25
g = 0.5
a = 20/N
m=mPhys-g*g*a/8
# x0 = 15
# x0 = 20
x0=5


main(N,m,g,a,x0)
