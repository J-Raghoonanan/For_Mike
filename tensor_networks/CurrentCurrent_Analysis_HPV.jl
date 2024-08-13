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

function createJOperator(omega,n, a)
    #=
    Construct the J^\omega (n) operator
    To be used for calculating  <J(x) J(0)>

    Args:
    omega: tensor index
    n: staggered lattice site
    a: staggered lattice spacing
    =#

    operator = OpSum()
    if omega==0
        if iseven(n)
            operator .+= 1/a , "S+", n, "S-", n
            operator .+= 1/a, "S+", n-1, "S-", n-1
        else
            operator .+= 1/a , "S+", n+1, "S-", n+1
            operator .+= 1/a, "S+", n, "S-", n
            # println("Error odd indexation")
            # exit()
        end
    elseif omega==1
        im_unit = 1.0im
        if iseven(n)
            operator -= im_unit/a, "S+",n,"S-",n-1
            operator += im_unit/a, "S+",n-1,"S-",n
        else
            operator -= im_unit/a, "S+",n+1,"S-",n
            operator += im_unit/a, "S+",n,"S-",n+1
            # println("Error odd indexation")
            # exit()
        end
    elseif omega==5
        if iseven(n)
            operator .+= 1/a , "S+", n-1, "S-", n
            operator .-= 1/a, "S+", n, "S-", n-1
        else
            operator .+= 1/a , "S+", n, "S-", n+1
            operator .-= 1/a, "S+", n+1, "S-", n
        end
    end

    return operator
end


function getEvolvedStates(alpha,beta,a,n,m ,MPO_Ham, groundPsi,initialLeftState, initialRightState, nsweeps, maxdim, cutoff, dt, sites)
    #=
    =#

    # n=2*n
    # m = 2*m
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
        psiT = groundPsi
    else
        psiT = timeEvolution(MPO_Ham, initialLeftState ; numSweeps=1, maxDim=maxdim, tolerance=cutoff ,timeT=dt, sites=sites, backend=0)
    end
    #backend = 0 ==> normalize = true

    Psi_left_op = createJOperator(alpha,n,a)
    MPO_Psi_left = MPO(Psi_left_op, sites)
    Psi_left = apply(MPO_Psi_left, psiT; fixedCutOff)
    # Psi_left = normalize!(Psi_left)

    if nsweeps==0
        Psi_right_op = createJOperator(beta,m,a)
        MPO_Psi_right = MPO(Psi_right_op, sites)
        Psi_right = apply(MPO_Psi_right, groundPsi; fixedCutOff)
    else
        Psi_right = timeEvolution(MPO_Ham, initialRightState; numSweeps=1, maxDim=maxdim, tolerance=cutoff ,timeT=dt, sites=sites, backend=1)
    end
    #backend = 0 ==> normalize = false

    return psiT, Psi_left, Psi_right
end



# function Correlator_kplus(N,MPO_Ham, groundPsi, asv; x0, cparameter, kplus, nsweeps, maxdim, cutoff, sites)
#     #=
#     Compute the correlator C(k+)
#     =#
#     println()
#     println("Total Kplus: $kplus")
#     c=cparameter
#     totaltime= c*sqrt(2)/kplus
#     dt = totaltime/nsweeps

#     println("Total time: $totaltime ... Total time steps: $nsweeps ... dt value= $dt")

#     time_array = zeros(Float64, nsweeps+1)
#     x_array= collect(x0:2:N-x0+1)
#     correlator00Vector = zeros(Complex{Float64}, length(time_array), length(x_array))
#     correlator01Vector = zeros(Complex{Float64}, length(time_array), length(x_array))
#     correlator10Vector = zeros(Complex{Float64}, length(time_array), length(x_array))
#     correlator11Vector = zeros(Complex{Float64}, length(time_array), length(x_array))

#     # println("Time array: $time_array")
#     println("x array: $x_array")
#     println("Loop size= $(length(time_array)* length(x_array))")

#     t_Index=0
#     x_Index=0

#     for nx in x_array
#         x_Index+=1

#         prevPsiRight00Even = groundPsi
#         prevPsiRight00Odd = groundPsi
#         prevPsiRight01Even = groundPsi
#         prevPsiRight01Odd = groundPsi
#         prevPsiRight10Even = groundPsi
#         prevPsiRight10Odd = groundPsi
#         prevPsiRight11Even = groundPsi
#         prevPsiRight11Odd = groundPsi
#         prevPsiLeftEven = groundPsi
#         prevPsiLeftOdd = groundPsi

#         nx_odd = nx
#         nx_even = nx+1
#         for _ in time_array
#             t_Index+=1
#             nsweep = t_Index-1
#             time_step = dt*nsweep # want to start at time=0, but Julia is 1-indexed
#             println("t index: ",t_Index, " x_index: ", x_Index)

#             # ODD Site Calculation
#             println("Running nt: $time_step and nx odd: $nx_odd, and nx0: $x0, and delta x: $(nx_odd-x0)")

#             newPsiLeftOdd,Psi_left00,Psi_right00_odd = getEvolvedStates(0,0,a,nx_odd,x0,MPO_Ham,groundPsi,prevPsiLeftOdd,prevPsiRight00Odd,nsweep,maxdim,cutoff,dt,sites)
#             newPsiLeftOdd,Psi_left01,Psi_right01_odd = getEvolvedStates(0,1,a,nx_odd,x0,MPO_Ham,groundPsi,prevPsiLeftOdd,prevPsiRight01Odd,nsweep,maxdim,cutoff,dt,sites)
#             newPsiLeftOdd,Psi_left10,Psi_right10_odd = getEvolvedStates(1,0,a,nx_odd,x0,MPO_Ham, groundPsi,prevPsiLeftOdd,prevPsiRight10Odd,nsweep,maxdim,cutoff,dt,sites)
#             newPsiLeftOdd,Psi_left11,Psi_right11_odd = getEvolvedStates(1,1,a,nx_odd,x0,MPO_Ham,groundPsi,prevPsiLeftOdd,prevPsiRight11Odd,nsweep,maxdim,cutoff,dt,sites)

#             # Fourier factors
#             f1 = exp(im*kplus*nsweep*dt/sqrt(2)) #adding phase kplus-phase
#             f2 = exp(-im*(kplus)*(nx_odd-x0)*asv/sqrt(2))
#             correlator00Vector[t_Index, x_Index]=asv*dt*f1*f2* inner(Psi_left00, Psi_right00_odd)
#             correlator01Vector[t_Index, x_Index]=asv*dt*f1*f2* inner(Psi_left01, Psi_right01_odd)
#             correlator10Vector[t_Index, x_Index]=asv*dt*f1*f2* inner(Psi_left10, Psi_right10_odd)
#             correlator11Vector[t_Index, x_Index]=asv*dt*f1*f2* inner(Psi_left11, Psi_right11_odd)

#             # EVEN Site Calculation
#             println("Running nt: $time_step and nx even: $nx_even, and nx0: $x0, and delta x: $(nx_even-x0)")
#             println()

#             newPsiLeftEven,Psi_left00,Psi_right00_even = getEvolvedStates(0,0,a,nx_odd,x0,MPO_Ham,groundPsi,prevPsiLeftEven,prevPsiRight00Even,nsweep,maxdim,cutoff,dt,sites)
#             newPsiLeftEven,Psi_left01,Psi_right01_even = getEvolvedStates(0,1,a,nx_odd,x0,MPO_Ham,groundPsi,prevPsiLeftEven,prevPsiRight01Even,nsweep,maxdim,cutoff,dt,sites)
#             newPsiLeftEven,Psi_left10,Psi_right10_even = getEvolvedStates(1,0,a,nx_odd,x0,MPO_Ham, groundPsi,prevPsiLeftEven,prevPsiRight10Even,nsweep,maxdim,cutoff,dt,sites)
#             newPsiLeftEven,Psi_left11,Psi_right11_even = getEvolvedStates(1,1,a,nx_odd,x0,MPO_Ham,groundPsi,prevPsiLeftEven,prevPsiRight11Even,nsweep,maxdim,cutoff,dt,sites)
#             f2 = exp(-im*(kplus)*(nx_even-x0)*asv/sqrt(2))
#             correlator00Vector[t_Index, x_Index]+=asv*dt*f1*f2* inner(Psi_left00, Psi_right00_even)
#             correlator01Vector[t_Index, x_Index]+=asv*dt*f1*f2* inner(Psi_left01, Psi_right01_even)
#             correlator10Vector[t_Index, x_Index]+=asv*dt*f1*f2* inner(Psi_left10, Psi_right10_even)
#             correlator11Vector[t_Index, x_Index]+=asv*dt*f1*f2* inner(Psi_left11, Psi_right11_even)


#             prevPsiRight00Even = Psi_right00_even
#             prevPsiRight00Odd = Psi_right00_odd
#             prevPsiRight01Even = Psi_right01_even
#             prevPsiRight01Odd = Psi_right01_odd
#             prevPsiRight10Even = Psi_right10_even
#             prevPsiRight10Odd = Psi_right10_odd
#             prevPsiRight11Even = Psi_right11_even
#             prevPsiRight11Odd = Psi_right11_odd
#             prevPsiLeftEven = newPsiLeftEven
#             prevPsiLeftOdd = newPsiLeftOdd
#         end
#         t_Index=0
#     end

#     corr00 = 0.0 + 0.0* 1.0im # Initialize the total sum
#     corr01 = 0.0 + 0.0* 1.0im
#     corr10 = 0.0 + 0.0* 1.0im
#     corr11 = 0.0 + 0.0* 1.0im
#     for i in 1:length(time_array)
#         for j in 1:length(x_array)
#             corr00  += correlator00Vector[i, j]
#             corr01  += correlator01Vector[i, j]
#             corr10  += correlator10Vector[i, j]
#             corr11  += correlator11Vector[i, j]
#         end
#     end

#     return corr00,corr01,corr10,corr11
# end






function Correlator_HPV(N,MPO_Ham, groundPsi, asv; x0, cparameter,n0,n1, totaltime, nsweeps, maxdim, cutoff, sites)
    #=
    Compute the correlator
    =#
    println()
    c=cparameter
    # totaltime= c*sqrt(2)/kplus
    dt = totaltime/nsweeps

    println("Total time: $totaltime ... Total time steps: $nsweeps ... dt value= $dt")

    time_array = zeros(Float64, nsweeps+1)
    x_array= collect(x0:2:N-x0+1)
    correlator00Vector = zeros(Complex{Float64}, length(time_array), length(x_array))
    correlator01Vector = zeros(Complex{Float64}, length(time_array), length(x_array))
    correlator10Vector = zeros(Complex{Float64}, length(time_array), length(x_array))
    correlator11Vector = zeros(Complex{Float64}, length(time_array), length(x_array))

    # println("Time array: $time_array")
    println("x array: $x_array")
    println("Loop size= $(length(time_array)* length(x_array))")

    t_Index=0
    x_Index=0

    for nx in x_array
        x_Index+=1

        prevPsiRight00Even = groundPsi
        prevPsiRight00Odd = groundPsi
        prevPsiRight01Even = groundPsi
        prevPsiRight01Odd = groundPsi
        prevPsiRight10Even = groundPsi
        prevPsiRight10Odd = groundPsi
        prevPsiRight11Even = groundPsi
        prevPsiRight11Odd = groundPsi
        prevPsiLeftEven = groundPsi
        prevPsiLeftOdd = groundPsi

        nx_odd = nx
        nx_even = nx+1
        for _ in time_array
            t_Index+=1
            nsweep = t_Index-1
            time_step = dt*nsweep # want to start at time=0, but Julia is 1-indexed
            println("t index: ",t_Index, " x_index: ", x_Index)

            # ODD Site Calculation
            println("Running nt: $time_step and nx odd: $nx_odd, and nx0: $x0, and delta x: $(nx_odd-x0)")

            newPsiLeftOdd,Psi_left00,Psi_right00_odd = getEvolvedStates(0,0,a,nx_odd,x0,MPO_Ham,groundPsi,prevPsiLeftOdd,prevPsiRight00Odd,nsweep,maxdim,cutoff,dt,sites)
            newPsiLeftOdd,Psi_left01,Psi_right01_odd = getEvolvedStates(0,1,a,nx_odd,x0,MPO_Ham,groundPsi,prevPsiLeftOdd,prevPsiRight01Odd,nsweep,maxdim,cutoff,dt,sites)
            newPsiLeftOdd,Psi_left10,Psi_right10_odd = getEvolvedStates(1,0,a,nx_odd,x0,MPO_Ham, groundPsi,prevPsiLeftOdd,prevPsiRight10Odd,nsweep,maxdim,cutoff,dt,sites)
            newPsiLeftOdd,Psi_left11,Psi_right11_odd = getEvolvedStates(1,1,a,nx_odd,x0,MPO_Ham,groundPsi,prevPsiLeftOdd,prevPsiRight11Odd,nsweep,maxdim,cutoff,dt,sites)

            # Fourier factors
            # f1 = exp(im*kplus*nsweep*dt/sqrt(2)) #adding phase kplus-phase
            # f2 = exp(-im*(kplus)*(nx_odd-x0)*asv/sqrt(2))
            f1 = dt*exp(im*2*pi*n0*(t_Index/totaltime)) #adding phase kplus-phase
            f2 = (nx_odd-x0)*exp(-im*2*pi*n1*(nx_odd/N)*asv)
            correlator00Vector[t_Index, x_Index]=f1*f2* inner(Psi_left00, Psi_right00_odd)
            correlator01Vector[t_Index, x_Index]=f1*f2* inner(Psi_left01, Psi_right01_odd)
            correlator10Vector[t_Index, x_Index]=f1*f2* inner(Psi_left10, Psi_right10_odd)
            correlator11Vector[t_Index, x_Index]=f1*f2* inner(Psi_left11, Psi_right11_odd)

            # EVEN Site Calculation
            println("Running nt: $time_step and nx even: $nx_even, and nx0: $x0, and delta x: $(nx_even-x0)")
            println()

            newPsiLeftEven,Psi_left00,Psi_right00_even = getEvolvedStates(0,0,a,nx_odd,x0,MPO_Ham,groundPsi,prevPsiLeftEven,prevPsiRight00Even,nsweep,maxdim,cutoff,dt,sites)
            newPsiLeftEven,Psi_left01,Psi_right01_even = getEvolvedStates(0,1,a,nx_odd,x0,MPO_Ham,groundPsi,prevPsiLeftEven,prevPsiRight01Even,nsweep,maxdim,cutoff,dt,sites)
            newPsiLeftEven,Psi_left10,Psi_right10_even = getEvolvedStates(1,0,a,nx_odd,x0,MPO_Ham, groundPsi,prevPsiLeftEven,prevPsiRight10Even,nsweep,maxdim,cutoff,dt,sites)
            newPsiLeftEven,Psi_left11,Psi_right11_even = getEvolvedStates(1,1,a,nx_odd,x0,MPO_Ham,groundPsi,prevPsiLeftEven,prevPsiRight11Even,nsweep,maxdim,cutoff,dt,sites)
            f2 = (nx_even-x0)*exp(-im*2*pi*n1*(nx_even/N)*asv)
            correlator00Vector[t_Index, x_Index]+=f1*f2* inner(Psi_left00, Psi_right00_even)
            correlator01Vector[t_Index, x_Index]+=f1*f2* inner(Psi_left01, Psi_right01_even)
            correlator10Vector[t_Index, x_Index]+=f1*f2* inner(Psi_left10, Psi_right10_even)
            correlator11Vector[t_Index, x_Index]+=f1*f2* inner(Psi_left11, Psi_right11_even)


            prevPsiRight00Even = Psi_right00_even
            prevPsiRight00Odd = Psi_right00_odd
            prevPsiRight01Even = Psi_right01_even
            prevPsiRight01Odd = Psi_right01_odd
            prevPsiRight10Even = Psi_right10_even
            prevPsiRight10Odd = Psi_right10_odd
            prevPsiRight11Even = Psi_right11_even
            prevPsiRight11Odd = Psi_right11_odd
            prevPsiLeftEven = newPsiLeftEven
            prevPsiLeftOdd = newPsiLeftOdd
        end
        t_Index=0
    end

    corr00 = 0.0 + 0.0* 1.0im # Initialize the total sum
    corr01 = 0.0 + 0.0* 1.0im
    corr10 = 0.0 + 0.0* 1.0im
    corr11 = 0.0 + 0.0* 1.0im
    for i in 1:length(time_array)
        for j in 1:length(x_array)
            corr00  += correlator00Vector[i, j]
            corr01  += correlator01Vector[i, j]
            corr10  += correlator10Vector[i, j]
            corr11  += correlator11Vector[i, j]
        end
    end

    return corr00,corr01,corr10,corr11
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
    nsweeps = 10
    maxdim = 350
    cutoff = 1E-10
    totaltime = 60
    # dt = totaltime / nsweeps
    # kplus_range = 25.0:25.0:500.0
    # kplus_vec = Float64.(collect(kplus_range))

    correlatorVector00 = Vector{Complex{Float64}}(undef, totaltime,N)
    correlatorVector01 = Vector{Complex{Float64}}(undef, totaltime,N)
    correlatorVector10 = Vector{Complex{Float64}}(undef, totaltime,N)
    correlatorVector11 = Vector{Complex{Float64}}(undef, totaltime,N)
    cparameter = 1

    # cutoff = 1E-10
    maxdim = 200
    cutoff = 1E-8
    # for (counter, kplus) in enumerate(kplus_vec)
    #     correlatorVector00[counter], correlatorVector01[counter], correlatorVector10[counter], correlatorVector11[counter] = Correlator_kplus(N,MPO_Ham, groundPsi, asv; x0, cparameter, kplus, nsweeps, maxdim, cutoff, sites)
    # end

    for n0 in 1:totaltime
        for n1 in 1:N
            correlatorVector00[n0,n1], correlatorVector01[n0,n1], correlatorVector10[n0,n1], correlatorVector11[n0,n1] = Correlator_HPV(N,MPO_Ham, groundPsi, asv; x0, cparameter, n0,n1,totaltime, nsweeps, maxdim, cutoff, sites)
        end
    end
    ####################################################################################################################
    plot(kplus_vec, real(correlatorVector00),color="red",label=L"Re(C^{00}(k^+))")
    plot!(kplus_vec, real(correlatorVector01),color="blue",label=L"Re(C^{01}(k^+))")
    plot!(kplus_vec, real(correlatorVector10),color="green",label=L"Re(C^{10}(k^+))")
    plot!(kplus_vec, real(correlatorVector11),color="purple",label=L"Re(C^{11}(k^+))")
    # log scale on y-axis
    # plot!(yaxis=:log10)
    xlabel!(L"k^+ \, [{\rm GeV}] ")
    ylabel!(L"\mathcal{C}(k^+)")
    title_str = "CurrCurr_Real_Corr_HPV_N_$(N)_m_$(m)_g_$(g)_a_$(asv)"
    filename_str = "Plots/Correlators/CurrentCurrent/HPV/Julia_plots/$title_str.pdf"
    title!(title_str)
    savefig(filename_str)

    plot(kplus_vec, imag(correlatorVector00),color="red",label=L"Im(C^{00}(k^+))")
    plot!(kplus_vec, imag(correlatorVector01),color="blue",label=L"Im(C^{01}(k^+))")
    plot!(kplus_vec, imag(correlatorVector10),color="green",label=L"Im(C^{10}(k^+))")
    plot!(kplus_vec, imag(correlatorVector11),color="purple",label=L"Im(C^{11}(k^+))")
    # log scale on y-axis
    # plot!(yaxis=:log10)
    xlabel!(L"k^+ \, [{\rm GeV}] ")
    ylabel!(L"\mathcal{C}(k^+)")
    title_str = "CurrCurr_Imag_Corr_HPV_N_$(N)_m_$(m)_g_$(g)_a_$(asv)"
    filename_str = "Plots/Correlators/CurrentCurrent/HPV/Julia_plots/$title_str.pdf"
    title!(title_str)
    savefig(filename_str)

    data_dict = Dict("k_plus" => kplus_vec,
                 "correlator00Vector_imag" => imag(correlatorVector00),
                 "correlator01Vector_imag" => imag(correlatorVector01),
                 "correlator10Vector_imag" => imag(correlatorVector10),
                 "correlator11Vector_imag" => imag(correlatorVector11),
                 "correlator00Vector_real" => real(correlatorVector00),
                 "correlator01Vector_real" => real(correlatorVector01),
                 "correlator10Vector_real" => real(correlatorVector10),
                 "correlator11Vector_real" => real(correlatorVector11))
    df = DataFrame(data_dict)
    CSV.write("Plots/Correlators/CurrentCurrent/HPV/CSV/CurrCurr_Corr_HPV_N_$(N)_m_$(m)_g_$(g)_a_$(asv).csv", df)
    return
end

# Set 5
N=20
x0=5

# N=50
# x0=10
mPhys = 0.25
g = 0.5
a = 1.0
m=mPhys-g*g*a/8


main(N,m,g,a,x0)
