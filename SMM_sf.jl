using Pkg
using LinearAlgebra, CSV, DataFrames, Distributions, Parameters, Optim

@with_kw struct parameters_SMM
    θ::Float64 = -2
    α::Array = [.13 .13]
    αhat::Array = [0.33 0.33]
    β::Array = [1 1 1 1 1; 1 1 1 1 1]
    σ::Array = [1 1 1 1 1; 1 1 1 1 1]
    σε::Array = [50 50]
    λ::Array = [0.5 0.5 0.5; 0.5 0.5 0.5]
    y::Array = [100000 100000]
    ρ::Float64 = 1.1341810237625438
    γ::Float64 = 1.416
    ψ::Float64 = 0.868
    ϕ::Float64 = 5000
    κ::Float64 = 15000
end

para_SMM = parameters_SMM()

cd("/Users/bigfoot13770/Documents/UO ECON PROGRAM/FUTILITY/JMP_RCsf")
# Read in sf data
data = CSV.read("sf_data.csv",DataFrame)

function makeR(data,year,times, varR, para)
        @unpack θ, ρ = para
        R = zeros(size(data)[1], 4,times)
        for t in 1:times
                RC = zeros(size(data)[1],1,1);
                RC_C = zeros(size(data)[1],1,1);
                RX = zeros(size(data)[1],1,1);
                RO = zeros(size(data)[1],1,1);
                if year ==2010
                    for j in 1:size(data)[1]
                        # size distributions
                        dΛC = truncated(Normal(data.Csize10[j],data.Csigma10[j]),0, Inf)
                        dΛC_C = truncated(Normal(data.Csize10[j],data.Csigma10[j]),0, Inf)
                        dΛX = truncated(Normal(data.Xsize10[j],data.Xsigma10[j]),0, Inf)
                        dΛO = truncated(Normal(data.Osize10[j],data.Osigma10[j]),0, Inf)

                        # rental distributions
                        dC = truncated(Normal(data.Crent10[j],data.CsigmaR10[j]),data.Crent10[j] - varR,data.Crent10[j] + varR);
                        dC_C = truncated(Normal(data.Crent10[j],data.CsigmaR10[j]),data.Crent10[j] - varR,data.Crent10[j] + varR);
                        dX = truncated(Normal(data.Xrent10[j],data.XsigmaR10[j]),data.Xrent10[j] - varR,data.Xrent10[j] + varR);
                        dO = truncated(Normal(data.mortgage10[j],max(data.Mortgagesigma10[j])),data.mortgage10[j] - varR,data.mortgage10[j] + varR);
                            RC[j] = sum((rand(dΛC, data.Cstruct10[j]).*rand(dC,data.Cstruct10[j]).^(1+θ) ) )^(1/(1+θ));
                            RC_C[j] = ρ*sum((rand(dΛC, data.Cstruct10[j]).*rand(dC_C,data.Cstruct10[j],).^(1+θ) ) )^(1/(1+θ));
                            RX[j] = sum(rand(dΛX, data.Xstruct10[j]).*(rand(dX,min(data.Xstruct10[j],1)).^(1+θ) ) )^(1/(1+θ));
                            RO[j] = sum((rand(dΛC, data.Ostruct10[j]).*rand(dO,min(data.Ostruct10[j],1)).^(1+θ) ) )^(1/(1+θ));
                    end # j-loop
                else
                    for j in 1:size(data)[1]
                            # size distributions
                            dΛC = truncated(Normal(data.Csize10[j],data.Csigma10[j]),0, Inf)
                            dΛC_C = truncated(Normal(data.Csize10[j],data.Csigma10[j]),0, Inf)
                            dΛX = truncated(Normal(data.Xsize10[j],data.Xsigma10[j]),0, Inf)
                            dΛO = truncated(Normal(data.Osize10[j],data.Osigma10[j]),0, Inf)
                            # Rental distributions
                        dC = truncated(Normal(data.Crent19[j],data.CsigmaR10[j]),data.Crent19[j] - varR,data.Crent19[j] + varR);
                        dC_C = truncated(Normal(data.Crent10[j],data.CsigmaR10[j]),data.Crent10[j] - varR,data.Crent10[j] + varR);
                        dX = truncated(Normal(data.Xrent19[j],data.XsigmaR10[j]),data.Xrent19[j] - varR,data.Xrent19[j] + varR);
                        dO = truncated(Normal(data.mortgage19[j],max(data.Mortgagesigma19[j],1e-3)),data.mortgage19[j] - varR,data.mortgage19[j] + varR);
                            RC[j] = sum((rand(dΛC,data.Cstruct19[j]) .*rand(dC,data.Cstruct19[j]).^(1+θ) ) )^(1/(1+θ));
                            RC_C[j] = ρ*sum((rand(dΛC, data.Cstruct10[j]).*rand(dC_C,data.Cstruct10[j],).^(1+θ) ) )^(1/(1+θ));
                            RX[j] = sum(rand(dΛC,data.Xstruct19[j]) .*(rand(dX,data.Xstruct19[j]).^(1+θ) ) )^(1/(1+θ));
                            RO[j] = sum((rand(dΛC,data.Ostruct19[j]) .*rand(dO,data.Ostruct19[j]).^(1+θ) ) )^(1/(1+θ));
                    end # j loop

            end # if statement
            R[:,:,t] = [RC RX RO RC_C]
    end # t loop

        R = mean(R,dims = 3)[:,:,1]

        return R
end

# Household preference draws for neighborhood attributes
function drawν(N)
    νH = randn(N[1,2],5)
    νL = randn(N[2,2],5)
    ν = []
    push!(ν,νH)
    push!(ν,νL)
    return ν
end

function hhendow(data, R, N)
        # Set up price matrix
        rents = repeat([R[:,1] ; R[:,2] ; R[:,3]]', Int(N))
        # Loop by tract and endow some agents with RC_C by replacement
        index = 1
        poptrack = 0
        for x in 1:size(data.CshareH10)[1]
                # Find number of RCs to replace
                pop = round(Int,N[1]*data.CshareH10[x])
                # Reassign values
                rents[index:(pop + poptrack),x] .= R[x,4]
                # Update indexes
                index += pop
                poptrack += pop
        end

        return rents
end

function make_mu(consum, data, para_SMM)
        @unpack α, σ, β = para_SMM
        # Tract characteristics X
        XH = [data.cbd./1000 data.Hcommute./1000 data.parks data.park_share data.water./1000]
        XL = [data.cbd./1000 data.Lcommute./1000 data.parks data.park_share data.water./1000]

        # Assemble inner utility
        μH = -α[1].*consum[1][:,6:size(consum[1])[2]] +  repeat(consum[1][:,1:5]*Diagonal(exp.(σ[1,:]))*XH',outer = (1,3))
        μL = -α[2].*consum[2][:,6:size(consum[2])[2]] +  repeat(consum[2][:,1:5]*Diagonal(exp.(σ[2,:]))*XH',outer = (1,3))
        μ = [[μH] [μL]]
        return μ
end

function endowε(N,data,para_SMM)
    @unpack σε = para_SMM
    dH = Gumbel(0,σε[1])
    dL = Gumbel(0,σε[2])
    ε = []
    push!(ε, rand(dH,N[1,2],3*size(data)[1]))
    push!(ε, rand(dH,N[2,2],3*size(data)[1]))
    return ε
end

H = [1286 1738]
L = [1668 1770]
N = 1000 .*[H;L]

# Function for calculating household shares
function hhsort(pop0, μ, R, ξ, ε, data, para_SMM)
    @unpack θ, α, β, σ, λ = para_SMM
    J = size(data)[1]
    # Find utility maximizing j-η index
    XH = repeat([data.cbd data.Hcommute data.parks data.park_share data.water], outer = 3)
    XL = repeat([data.cbd data.Lcommute data.parks data.park_share data.water], outer = 3)
    indexH = findmax((XH*β[1,:])' + [ξ[:,1,1]; ξ[:,2,1]; ξ[:,3,1]]' .+ μ[1] .+ ε[1], dims = 2)[2]
    indexL = findmax((XL*β[1,:])' + [ξ[:,1,2]; ξ[:,2,2]; ξ[:,3,2]]' .+ μ[2] .+ ε[2], dims = 2)[2]
    # Loop by individuals to -- man this shit is getting old, I really have to loop through 3 million people just to figure out what alternative is best for them like come on
    optimlocH = zeros(size(indexH)[1])
    popH = zeros(3*J)
    for i in 1:size(indexH)[1]
        optimlocH[i] = indexH[i][2]
        popH[indexH[i][2]] = popH[indexH[i][2]] + 1
    end
    optimlocL = zeros(size(indexL)[1])
    popL = zeros(3*J)
    for i in 1:size(indexL)[1]
        optimlocL[i] = indexL[i][2]
        popL[indexL[i][2]] = popL[indexL[i][2]] + 1
    end
    # Assemble pop matrix and locs
    pop = zeros(J,3,2)
    pop[:,:,1] = [popH[1:J] popH[(J+1):2*J] popH[(2*J+1):end]]
    pop[:,:,2] = [popL[1:J] popL[(J+1):2*J] popL[(2*J+1):end]]
    # loop through locations to count movers
    moversH = zeros(3*J)
    moversL = zeros(3*J)
    pop0H = [pop0[:,1,1]; pop0[:,2,1]; pop0[:,3,1]]
    pop0L = [pop0[:,1,2]; pop0[:,2,2]; pop0[:,3,2]]
    trackH = 0
    trackH1 = 1
    trackL = 0
    trackL1 = 1
    for x in 1:3*J
        neighborhood = x÷3
        println("HH Action, neighborhood $neighborhood")
        moversH[x] = count(i -> i != x, optimlocH[Int(trackH1):Int(pop0H[x] + trackH)])
        trackH += pop0H[x]
        trackH1 += pop0H[x]
        moversL[x] = count(i -> i != x, optimlocL[Int(trackL1):Int(pop0L[x] + trackL)])
        trackL += pop0L[x]
        trackL1 += pop0L[x]
    end
    movers = zeros(J,3,2)
    movers[:,:,1] = [moversH[1:J] moversH[(J+1):2*J] moversH[(2*J+1):end]]
    movers[:,:,2] = [moversL[1:J] moversL[(J+1):2*J] moversL[(2*J+1):end]]

    return pop, movers
end

# Function to endow landowners with initial conditions
function drawr0(data, lotto, varR)
    r0 = []

    for j in 1:size(data)[1]
        dCR = truncated(Normal(data.Crent10[j],data.CsigmaR10[j]),data.Crent10[j] - varR,data.Crent10[j] + varR);
        dCR_C = truncated(Normal(data.Crent10[j],data.CsigmaR10[j]),data.Crent10[j] - varR,data.Crent10[j] + varR);
        dXR = truncated(Normal(data.Xrent10[j],data.XsigmaR10[j]),data.Xrent10[j] - varR,data.Xrent10[j] + varR);
        dOR = truncated(Normal(data.mortgage10[j],max(data.Mortgagesigma10[j],1e-3)),data.mortgage10[j] - varR,data.mortgage10[j] + varR);
        # create L0η matricies to store initial conditions by η
        r0C = rand(dCR,max(Int(round(data.Cstruct10[j]-lotto*data.win[j])),0))

        if data.Cstruct10[j]-lotto*data.win[j] < 0
            r0C_C = zeros(0)
        else
            r0C_C = rand(dCR,Int(round(lotto*data.win[j])))
        end
        r0X = rand(dXR,data.Xstruct10[j])
        r0O = rand(dOR,data.Ostruct10[j])
        # Store in L0
        push!(r0, [[r0C] [r0X] [r0O] [r0C_C]])
    end

    return r0
end

function endowL0(r0, pop0, tol, para_SMM)
    @unpack θ, αhat, y = para_SMM
    #Posit initial M0
    M0 = 1e5 .* ones(size(r0)[1],3) # M0 only 3×J
    Mnew = []
    Mstore = zeros(size(r0)[1],3)
    CONT = true
    while CONT == true
        Mnew = zeros(size(r0)[1],3)
        for j in 1:size(r0)[1], η in 1:3
            # Combine η = 1,2 for initial period (doesn't matter)
            if η == 1
                Λ0 = M0[j,1].*[r0[j][1];r0[j][2]].^θ
                R = (sum(Λ0 .* [r0[j][1];r0[j][2]].^(θ+1)))^(1/(θ+1))
                Mnew[j,η] = pop0[j,1,1].*αhat[1]*y[1]./R^(θ+1) + pop0[j,1,2].*αhat[2]*y[2]./R^(θ+1)
            else
                Λ0 = M0[j,η].*r0[j][η].^θ
                R = (sum(Λ0 .* r0[j][η].^(θ+1)))^(1/(θ+1))
                Mnew[j,η] = pop0[j,η,1].*αhat[1]*y[1]./R^(θ+1) + pop0[j,η,2].*αhat[2]*y[2]./R^(θ+1)
            end # η instructions
            # check in Inf
            if Mnew[j,η] == Inf
                Mnew[j,η] = 0
            end
        end # j-η loop
        if findmax(abs.(Mnew - M0))[1] < tol
            M0 = Mnew
            CONT = false
        else
            #println(findmax(abs.(Mnew - M0))[1])
            M0 = 0.5.*(Mnew + Mstore)
            Mstore = Mnew
        end
    end # While loop
    L0 = []
    for j in 1:size(r0)[1]
        L0C = [r0[j][1] M0[j,1].*r0[j][1].^θ]
        L0C_C = [r0[j][4] M0[j,1].*r0[j][4].^θ]
        L0X = [r0[j][2] M0[j,2].*r0[j][2].^θ]
        L0O = [r0[j][3] M0[j,3].*r0[j][3].^θ]
        push!(L0, [[L0C] [L0X] [L0O] [L0C_C]])
    end # j loop
    return L0
end

r0 = drawr0(data, lotto, varR)
L0 = endowL0(r0, pop0, 0.01, para_SMM)

replace_0(v) = map(x -> x == 0 ? 1e-100 : x, v)
replace_0(0)
# Function for characterizing landowner reactions
function losort(L0, M, R0, movers, data, ς, Vmin, para_SMM)
    # Set up parameters
    @unpack θ,αhat,y,ρ,γ,ϕ,ψ,κ = para_SMM
    ς = ς[1]
    Psi = data.Psi
    Γ = data.Gamma
    b = data.avgBO
    ϑ = 1/(1+θ - ψ*θ)
    d = Gumbel(0,ς)

    # Set up matrices to store results
    𝚲 = []
    𝐫 = []
    πV = zeros(size(data)[1])
    # Count 2010 vacant landowners Vj, number of those still vacant 2019
    Vj = max.(data.parcels - (data.Cstruct10 + data.Xstruct10 + data.Ostruct10),Vmin)
    Vbar = max.(Vj - (data.Xstruct19 + data.Ostruct19) + (data.Xstruct10 + data.Ostruct10),0)
    # Loop through by neighborhood, determine new prices, structures
    for j in 1:size(data)[1]
        println("LO Action, neighborhood $j")
        # Set up storage matrices
        ΛC = []; ΛC_C = []; ΛX = []; ΛO = [];
        rC = []; rC_C = []; rX = []; rO = [];
        # Iterate by η in each neighborhood
        for η in 1:4
            if η == 1 # Controlled C
                # Draw Gumbel shocks
                u = rand(d,size(L0[j][η])[1],3)
                Λ0 = L0[j][η][:,2]
                r0 = L0[j][η][:,1]
                # Calculate φ
                φ = max.((Λ0 .- (movers[j,1,1]*αhat[1]*y[1]/(R0[j,η].^(1+θ)).*r0.^θ + movers[j,1,2]*αhat[2]*y[2]/(R0[j,η].^(1+θ)).*r0.^θ))./Λ0, 0)
                # Set up matrix to test different adjustment possibilities
                r_test = repeat(collect(500:1:3000)',size(Λ0)[1])
                Λfree = M[j,1] .*r_test .^θ
                adj_profit = φ.*ρ.*Λ0 .+ (Γ[j]*(θ+1)/γ*θ).*((Λfree.+φ.*Λ0)./Λ0).^γ .- data.avgBO[j].*Λ0
                # Near optimal objectives
                maxi = findmax(adj_profit,dims=2)
                rstar = zeros(size(Λ0)[1],1)
                # Loop by k for cartesian indices
                for k in 1:size(Λ0)[1]
                    rstar[k] = r_test[k,Int(maxi[2][k][2])]
                end # k-loop

                # Caclulate profits under each alternative
                Stasis = φ.*ρ.*r0.*Λ0 + ((1 .-φ).*Λ0./M[j,η]).^((θ+1)/θ).*M[j,η]
                Adjustment = maxi[1][:,1]
                # Set up renovation to η = X/0 = 3/4
                RenoX = (Psi[j]*M[j,2]^(ψ-1))^(ϑ*(θ+1))*M[j,2] .- ϕ.*Λ0 .- (Psi[j]*(1+θ)/(ψ*θ)).*(Psi[j]*M[j,2]^(ψ-1))^(ϑ*ψ) .- data.avgBO[j].*Λ0
                RenoO = (Psi[j]*M[j,3]^(ψ-1))^(ϑ*(θ+1))*M[j,3] .- ϕ.*Λ0 .- (Psi[j]*(1+θ)/(ψ*θ)).*(Psi[j]*M[j,3]^(ψ-1))^(ϑ*ψ) .- data.avgBO[j].*Λ0
                profits = [Stasis+u[:,1] Adjustment+u[:,1] RenoX+u[:,2] RenoO+u[:,3]]
                # Assmeble rents
                rents = [ρ.*r0 rstar repeat([(Psi[j]*M[j,2]^(ψ-1))^ϑ], size(Λ0)[1]) repeat([(Psi[j]*M[j,3]^(ψ-1))^ϑ], size(Λ0)[1])]
                Λ = [Λ0 φ.*Λ0+M[j,1].*rstar.^θ M[j,2].*rents[:,3].^θ  M[j,3].*rents[:,4].^θ]
                # use findmax() to figure out indicies for optimal housing choice, iterate through cartesian indicies cause Julia sucks sometimes to assemble Λ and r storage matricies
                index = findmax(profits,dims = 2)[2]
                ΛC = zeros(size(profits)[1],4)
                rC = zeros(size(profits)[1],4)
                for k in 1:size(profits)[1]
                    if Int(index[k][2]) == 1
                        ΛC[k,2] = φ[k]*Λ0[k]
                        ΛC[k,1] = (1-φ[k])*Λ0[k]
                        rC[k,2] = rents[k,1]
                        rC[k,1] = rstar[k]
                    elseif Int(index[k][2]) == 2
                        ΛC[k,2] = φ[k]*Λ0[k]
                        ΛC[k,1] = rstar[k]^θ*M[j,η]
                        rC[k,2] = rents[k,1]
                        rC[k,1] = rstar[k]
                    else
                        ΛC[k,Int(index[k][2])] = Λ[k,Int(index[k][2])]
                        rC[k,Int(index[k][2])] = rents[k,Int(index[k][2])]
                    end # if statement
                end # k loop
            elseif η == 2 # Exempt X
                # Draw Gumbel shocks
                u = rand(d,size(L0[j][η])[1],3)
                # Caclulate profits under each alternative
                Λ0 = L0[j][η][:,2]
                Stasis = (Λ0./M[j,η]).^((θ+1)/θ).*M[j,η]
                Adjustment = (M[j,η]*Γ[j]./Λ0.^γ).^((θ+1)/(1-γ*θ))*M[j,η] .- (Γ[j]*(1+θ)/(γ*θ)).* ((M[j,η]*Γ[j]./Λ0.^γ).^(1/(1-γ*θ))./Λ0).^γ .- data.avgBO[j].*Λ0
                # Set up renovation to η = O = 3
                Condo = (Λ0./M[j,3]).^((θ+1)/θ) .- (κ + data.avgBO[j]).*Λ0
                Reno = (Psi[j]*M[j,3]^(ψ-1))^(ϑ*(θ+1))*M[j,3] .- ϕ.*Λ0 .- (Psi[j]*(1+θ)/(ψ*θ)).*(Psi[j]*M[j,3]^(ψ-1))^(ϑ*ψ) .- data.avgBO[j].*Λ0
                profits = [Stasis+u[:,1] Adjustment+u[:,1] Condo+u[:,2] Reno+u[:,2]]
                rents = [(Λ0./M[j,η]).^(1/θ) (M[j,η]*Γ[j]./Λ0.^γ).^(1/(1-γ*θ)) (Λ0./M[j,2]).^(1/θ) repeat([(Psi[j]*M[j,3]^(ψ-1))^ϑ], size(Λ0)[1])]
                Λ = [rents[:,1:2].^θ.*M[j,η] rents[:,3:4].^θ.*M[j,3]]
                # use findmax() to figure out indicies for optimal housing choice, iterate through cartesian indicies cause Julia sucks sometimes to assemble Λ and r storage matricies
                index = findmax(profits,dims = 2)[2][:]
                ΛX = zeros(size(profits)[1],4)
                rX = zeros(size(profits)[1],4)
                for k in 1:size(profits)[1]
                    ΛX[k,Int(index[k][2])] = Λ[k,Int(index[k][2])]
                    rX[k,Int(index[k][2])] = rents[k,Int(index[k][2])]
                end # k loop

            elseif η == 3 # Owned O
                # Draw Gumbel shocks
                u = rand(d,size(L0[j][η])[1],2)
                # Caclulate profits under each alternative
                Λ0 = L0[j][η][:,2]
                Stasis = (Λ0./M[j,η]).^((θ+1)/θ).*M[j,η]
                Adjustment = (M[j,η]*Γ[j]./Λ0.^γ).^((θ+1)/(1-γ*θ))*M[j,η] .- (Γ[j]*(1+θ)/(γ*θ)).* ((M[j,η]*Γ[j]./Λ0.^γ).^(1/(1-γ*θ))./Λ0).^γ .- data.avgBO[j].*Λ0
                # Set up renovation to η = X = 2
                Reno = (Psi[j]*M[j,2]^(ψ-1))^(ϑ*(θ+1))*M[j,2] .- ϕ.*Λ0 .- (Psi[j]*(1+θ)/(ψ*θ)).*(Psi[j]*M[j,2]^(ψ-1))^(ϑ*ψ) .- data.avgBO[j].*Λ0
                profits = [Stasis+u[:,1] Adjustment+u[:,1] Reno+u[:,2]]
                rents = [(Λ0./M[j,η]).^(1/θ) (M[j,η]*Γ[j]./Λ0.^γ).^(1/(1-γ*θ)) repeat([(Psi[j]*M[j,2]^(ψ-1))^ϑ], size(Λ0)[1])]
                Λ = [rents[:,1:2].^θ.*M[j,η] rents[:,3].^θ.*M[j,2]]
                # use findmax() to figure out indicies for optimal housing choice, iterate through cartesian indicies cause Julia sucks sometimes to assemble Λ and r storage matricies
                index = findmax(profits,dims = 2)[2][:]
                ΛO = zeros(size(profits)[1],3)
                rO = zeros(size(profits)[1],3)
                for k in 1:size(profits)[1]
                    ΛO[k,Int(index[k][2])] = Λ[k,Int(index[k][2])]
                    rO[k,Int(index[k][2])] = rents[k,Int(index[k][2])]
                end # k loop
            else  # Controlled, lotto win C⋆
                # Draw Gumbel shocks
                u = rand(d,size(L0[j][η])[1],3)
                Λ0 = L0[j][η][:,2]
                r0 = L0[j][η][:,1]
                # Calculate φ
                φ = max.((Λ0 .- (movers[j,1,1]*αhat[1]*y[1]/(R0[j,η].^(1+θ)).*r0.^θ + movers[j,1,2]*αhat[2]*y[2]/(R0[j,η].^(1+θ)).*r0.^θ))./Λ0, 0)
                # Set up matrix to test different adjustment possibilities
                r_test = repeat(collect(500:0.1:3000)',size(Λ0)[1])
                Λfree = M[j,1] .*r_test .^θ
                adj_profit = φ.*ρ.*Λ0 .+ (Γ[j]*(θ+1)/γ*θ).*((Λfree.+φ.*Λ0)./Λ0).^γ .- data.avgBO[j].*Λ0
                # Near optimal objectives
                maxi = findmax(adj_profit,dims=2)
                rstar = zeros(size(Λ0)[1],1)
                # Loop by k for cartesian indices
                for k in 1:size(Λ0)[1]
                    rstar[k] = r_test[k,Int(maxi[2][k][2])]
                end # k-loop

                # Caclulate profits under each alternative
                Stasis = φ.*ρ.*r0.*Λ0 + ((1 .-φ).*Λ0./M[j,1]).^((θ+1)/θ).*M[j,1]
                Adjustment = maxi[1][:,1]
                # Set up renovation to η = X/0 = 3/4
                Condo = (Λ0./M[j,3]).^((θ+1)/θ) .- (κ + data.avgBO[j]).*Λ0
                RenoX = (Psi[j]*M[j,2]^(ψ-1))^(ϑ*(θ+1))*M[j,2] .- ϕ.*Λ0 .- (Psi[j]*(1+θ)/(ψ*θ)).*(Psi[j]*M[j,2]^(ψ-1))^(ϑ*ψ) .- data.avgBO[j].*Λ0
                RenoO = (Psi[j]*M[j,3]^(ψ-1))^(ϑ*(θ+1))*M[j,3] .- ϕ.*Λ0 .- (Psi[j]*(1+θ)/(ψ*θ)).*(Psi[j]*M[j,3]^(ψ-1))^(ϑ*ψ) .- data.avgBO[j].*Λ0
                profits = [Stasis+u[:,1] Adjustment+u[:,1] RenoX+u[:,2] Condo+u[:,3] RenoO+u[:,3]]
                # Assmeble rents
                rents = [ρ.*r0 rstar repeat([(Psi[j]*M[j,2]^(ψ-1))^ϑ], size(Λ0)[1]) (Λ0./M[j,3]).^(1/θ) repeat([(Psi[j]*M[j,3]^(ψ-1))^ϑ], size(Λ0)[1])]
                Λ = [Λ0 φ.*Λ0+M[j,1].*rstar.^θ M[j,2].*rents[:,3].^θ Λ0 M[j,3].*rents[:,5].^θ]
                # use findmax() to figure out indicies for optimal housing choice, iterate through cartesian indicies cause Julia sucks sometimes to assemble Λ and r storage matricies
                index = findmax(profits,dims = 2)[2]
                ΛC_C = zeros(size(profits)[1],5)
                rC_C = zeros(size(profits)[1],5)
                # Available C goes in 1, occupied C goes in 2
                for k in 1:size(profits)[1]
                    if Int(index[k][2]) == 1
                        ΛC_C[k,2] = φ[k]*Λ0[k]
                        ΛC_C[k,1] = (1-φ[k])*Λ0[k]
                        rC_C[k,2] = rents[k,1]
                        rC_C[k,1] = rstar[k]
                    elseif Int(index[k][2]) == 2
                        ΛC_C[k,2] = φ[k]*Λ0[k]
                        ΛC_C[k,1] = rstar[k]^θ*M[j,η]
                        rC_C[k,2] = rents[k,1]
                        rC_C[k,1] = rstar[k]
                    else
                        ΛC_C[k,Int(index[k][2])] = Λ[Int(index[k][2])]
                        rC_C[k,Int(index[k][2])] = rents[Int(index[k][2])]
                    end # if statement
                end # k loop
            end # η check
        end # η loop

        # New Housing:
        # Calculate πj(X/O) for new builders
        πX = (Psi[j]*M[j,2]^(ψ-1))^(ϑ*(θ+1))*M[j,2] .- (Psi[j]*(1+θ)/(ψ*θ)).*(Psi[j]*M[j,2]^(ψ-1))^(ϑ*ψ)
        πO = (Psi[j]*M[j,3]^(ψ-1))^(ϑ*(θ+1))*M[j,3] .- (Psi[j]*(1+θ)/(ψ*θ)).*(Psi[j]*M[j,3]^(ψ-1))^(ϑ*ψ)
        # Calculate πj(V)
        πV[j] = log(Vbar[j]/Vj[j]*(exp(πX/ς) + exp(πO/ς))/max(1 - Vbar[j]/Vj[j],0.0000001))
        # Draw Gumbel Shocks
        uV = rand(d,Vj[j],3)
        rentsV = [0 (Psi[j]*M[j,2]^(ψ-1))^ϑ  (Psi[j]*M[j,3]^(ψ-1))^ϑ]
        Λ = [0 rentsV[2]^θ*M[j,2] rentsV[3]^θ*M[j,3]]
        index = findmax([πV[j].+uV[:,1] πX.+uV[:,2] πO.+uV[:,3]], dims = 2)[2]
        # Loop by k
        ΛV = zeros(Vj[j],3)
        rV = zeros(Vj[j],3)
        for k in 1:Vj[j]
            ΛV[k,Int(index[k][2])] = Λ[Int(index[k][2])]
            rV[k,Int(index[k][2])] = rentsV[Int(index[k][2])]
        end #
        # Push results into bf matricies
        push!(𝚲,[[ΛC] [ΛX] [ΛO] [ΛC_C] [ΛV]])
        push!(𝐫,[[rC] [rX] [rO] [rC_C] [rV]])
    end # j loop
    return 𝐫, 𝚲, πV
end # losort



testr, testΛ, testπV = losort(L0, M, R0, movers, data, ς, 100, para_SMM)

# Assemble R given landowner actions
function assembleR(r, Λ, L0, para_SMM)
    @unpack θ, ρ = para_SMM
    # set up matrices to hold results
    RC = zeros(size(r)[1])
    RC_C = zeros(size(r)[1])
    RX = zeros(size(r)[1])
    RO = zeros(size(r)[1])

    # Loop by j
    for j in 1:size(r)[1]
        RC[j] = sum(sum(Λ[j][1][:,1][Λ[j][1][:,1].>0].*r[j][1][:,1][r[j][1][:,1].>0].^(θ+1)) + sum(Λ[j][4][:,1][Λ[j][4][:,1].>0].*r[j][4][:,1][r[j][4][:,1].>0].^(θ+1)))^(1/(θ+1))
        RX[j] = sum(sum(Λ[j][1][:,3][Λ[j][1][:,3].>0].*r[j][1][:,3][r[j][1][:,3].>0].^(θ+1)) + sum(Λ[j][4][:,3][Λ[j][4][:,3].>0].*r[j][4][:,3][r[j][4][:,3].>0].^(θ+1)) + sum(Λ[j][2][:,1][Λ[j][2][:,1].>0].*r[j][2][:,1][r[j][2][:,1].>0].^(θ+1)) + sum(Λ[j][2][:,2][Λ[j][2][:,2].>0].*r[j][2][:,2][r[j][2][:,2].>0].^(θ+1)) + sum(Λ[j][3][:,3][Λ[j][3][:,3].>0].*r[j][3][:,3][r[j][3][:,3].>0].^(θ+1)))^(1/(θ+1))
        RO[j] = sum(sum(Λ[j][1][:,3][Λ[j][1][:,3].>0].*r[j][1][:,3][r[j][1][:,3].>0].^(θ+1)) + sum(Λ[j][4][:,4][Λ[j][4][:,4].>0].*r[j][4][:,4][r[j][4][:,4].>0].^(θ+1)) + sum(Λ[j][4][:,5][Λ[j][4][:,5].>0].*r[j][4][:,5][r[j][4][:,5].>0].^(θ+1)) + sum(Λ[j][2][:,3][Λ[j][2][:,3].>0].*r[j][2][:,3][r[j][2][:,3].>0].^(θ+1)) + sum(Λ[j][2][:,4][Λ[j][2][:,4].>0].*r[j][2][:,4][r[j][2][:,4].>0].^(θ+1)) + sum(Λ[j][3][:,1][Λ[j][3][:,1].>0].*r[j][3][:,1][r[j][3][:,1].>0].^(θ+1)) + sum(Λ[j][3][:,2][Λ[j][3][:,2].>0].*r[j][3][:,2][r[j][3][:,2].>0].^(θ+1)))^(1/(θ+1))
        RC_C[j] = sum(sum(L0[j][1][:,2][Λ[j][1][:,1].>0].*r[j][1][:,2][r[j][1][:,2].>0].^(θ+1)) + sum(L0[j][4][:,2][Λ[j][4][:,1].>0].*r[j][4][:,2][r[j][4][:,2].>0].^(θ+1)))^(1/(θ+1))
    end # j loop
    R = [RC RX RO RC_C]
end

# pop contains share counts of households by location and housing type, pop[j,η,e]


function assembleM(pop, R, para_SMM)
    @unpack αhat, θ, y = para_SMM
    M = pop[:,:,1].*αhat[1]*y[1]./R[:,1:3].^(θ+1) + pop[:,:,2].*αhat[2]*y[2]./R[:,1:3].^(θ+1)
    return M # M is 4×J
end

function assembleR0(r0, L0, para_SMM)
    @unpack θ = para_SMM
    # set up matrices to hold results
    RC = zeros(size(r0)[1])
    RC_C = zeros(size(r0)[1])
    RX = zeros(size(r0)[1])
    RO = zeros(size(r0)[1])
    for j in 1:size(r0)[1]
        if size(r0[j][2][:])[1] == 0
            RC[j] = (sum(L0[j][1][:,2].*L0[j][1][:,1].^(θ+1)))^(1/(θ+1))
        else
            RC[j] = (sum(L0[j][1][:,2].*L0[j][1][:,1].^(θ+1)) + sum(L0[j][2][:,2].*L0[j][2][:,1].^(θ+1))^(1/(θ+1)))^(1/(θ+1))
        end # if statement
        RX[j] = sum(L0[j][3][:,2].*L0[j][3][:,1].^(θ+1))^(1/(θ+1))
        RO[j] = sum(L0[j][4][:,2].*L0[j][4][:,1].^(θ+1))^(1/(θ+1))
    end
    R = [RC RX RO]
    return R # M0 is only 3×J
end






fata = DataFrame(
    Cstruct10 = [10000; 10000],
    Xstruct10 = [10000; 10000],
    Ostruct10 = [10000; 10000],
    Cstruct19 = 2 .*[10000; 10000],
    Xstruct19 = 2 .*[10000; 10000],
    Ostruct19 = 2 .*[10000; 10000],
    parcels = [40000;80000],
    Gamma = [1;1],
    Psi = [1;1],
    avgBO = [1000;1000]
)


# Function to find equilibrium
function equil(ς,ν,ε,tol,N,ξ,lotto,varR,data,para_SMM)
    H10 = N[1,1]
    L10 = N[2,1]
    # Calculate pop0 from 2010 data
    pop0 = zeros(size(data)[1], 3, 2)
    pop0[:,:,1] = round.(H10.*[data.CshareH10+data.CshareH10_C data.XshareH10 data.OshareH10])
    pop0[:,:,2] = round.(L10.*[data.CshareL10+data.CshareL10_C data.XshareL10 data.OshareL10])
    # Draw r0's
    r0 = drawr0(data, lotto, varR)
    # Calculate rent indicies for 2010, use to find M0 and Λ0
    L0 = endowL0(r0, pop0, .01, para_SMM)
    R0 = assembleR0(r0, L0, para_SMM)
    R = makeR(data,2019,100, 1000, para_SMM)
    CONT = true
    while CONT == true
        # endow households based on previous simulation attempt
        rentsH = hhendow(data, R, N[1,2])
        rentsL = hhendow(data, R, N[2,2])
        consum = []
        push!(consum, [ν[1] rentsH])
        push!(consum, [ν[2] rentsL])
        μ = make_mu(consum, data, para_SMM)
        # Simulate household choices conditional on R from previous loop
        pop, movers = hhsort(pop0, μ, R, ξ, ε, data, para_SMM)
        # use the results of hhsort to assemble M matrix
        M = replace_0.(assembleM(pop, R, para_SMM))
        # Simulate landowner choices conditional on M
        r, Λ, boogers = losort(L0, M, R, movers, data, ς, 100, para_SMM)
        # Calculate rent indicies
        Rnew = assembleR(r, Λ, L0, para_SMM)
        # Convergence check
        if findmax(abs.(replace!(Rnew,Inf=>0) - replace!(R,Inf=>0)))[1] < tol
            R = Rnew
            CONT = false
        elseif isnan(findmax(abs.(Rnew - R))[1]) == true
            println(findmax(abs.(Rnew - R))[1])
            R = Rnew
            CONT = false
        else
            println(findmax(abs.(replace!(Rnew - R, Inf=>0)))[1])
            R = 0.5.*(Rnew + R)
        end # Convergence check
    end # while loop
    # endow households based on previous simulation attempt
    rentsH = hhendow(data, R, N[1,2])
    rentsL = hhendow(data, R, N[2,2])
    consum = []
    push!(consum, [ν[1] rentsH])
    push!(consum, [ν[2] rentsL])
    μ = make_mu(consum, data, para_SMM)
    # Simulate household choices conditional on R from previous loop
    pop, movers = hhsort(pop0, μ, R, ξ, ε, data, para_SMM)
    # use the results of hhsort to assemble M matrix
    M = replace_0.(assembleM(pop, R, para_SMM))
    # Simulate landowner choices conditional on M
    r, Λ, πV = losort(L0, M, R, movers, data, ς, 100, para_SMM)
    return pop, R, r, Λ, πV
end


function obj_SMM(ς,ν,ε,tol,N,ξ,lotto,varR,data,para_SMM)
    J = size(data)[1]
    # Find equilibrium given current ς guess
    pop, R, r, Λ = equil(ς,ν,ε,tol,N,ξ,lotto,varR,data,para_SMM)
    # Assemble simulated moments
    # Population moments order: j, η, e
    pop_sim = [[pop[:,1,1]; pop[:,2,1]; pop[:,3,1]]./N[1,2]; [pop[:,1,2]; pop[:,2,2]; pop[:,3,2]]./N[2,2]]
    # find median rents and number of units to match on
    # Order: j, η
    r_sim = zeros(J,3)
    Λ_sim = zeros(J,3)
    units_sim = zeros(J,3)
    for j in 1:J
        # extract rents
        r_sim[j,1] = median([r[j][1][:,1][r[j][1][:,1].>0]; r[j][1][:,2][r[j][1][:,2].>0]; r[j][4][:,1][r[j][4][:,1].>0]; r[j][4][:,2][r[j][4][:,2].>0];0])
        r_sim[j,2] = median([r[j][1][:,3][r[j][1][:,3].>0]; r[j][4][:,3][r[j][4][:,3].>0]; r[j][2][:,1][r[j][2][:,1].>0]; r[j][2][:,2][r[j][2][:,2].>0]; r[j][3][:,3][r[j][3][:,3].>0];0])
        r_sim[j,3] = median([r[j][1][:,4][r[j][1][:,4].>0]; r[j][4][:,4][r[j][4][:,4].>0]; r[j][4][:,5][r[j][4][:,5].>0]; r[j][2][:,3][r[j][2][:,3].>0]; r[j][2][:,4][r[j][2][:,4].>0]; r[j][3][:,1][r[j][3][:,1].>0]; r[j][3][:,2][r[j][3][:,2].>0];0])
        # extract number of units
        Λ_sim[j,1] = sum([Λ[j][1][:,1][Λ[j][1][:,1].>0]; Λ[j][1][:,2][Λ[j][1][:,2].>0]; Λ[j][4][:,1][Λ[j][4][:,1].>0]; Λ[j][4][:,2][Λ[j][4][:,2].>0];0])
        Λ_sim[j,2] = sum([Λ[j][1][:,3][Λ[j][1][:,3].>0]; Λ[j][4][:,3][Λ[j][4][:,3].>0]; Λ[j][2][:,1][Λ[j][2][:,1].>0]; Λ[j][2][:,2][Λ[j][2][:,2].>0]; Λ[j][3][:,3][Λ[j][3][:,3].>0];0])
        Λ_sim[j,3] = sum([Λ[j][1][:,4][Λ[j][1][:,4].>0]; Λ[j][4][:,4][Λ[j][4][:,4].>0]; Λ[j][4][:,5][Λ[j][4][:,5].>0]; Λ[j][2][:,3][Λ[j][2][:,3].>0]; Λ[j][2][:,4][Λ[j][2][:,4].>0]; Λ[j][3][:,1][Λ[j][3][:,1].>0]; Λ[j][3][:,2][Λ[j][3][:,2].>0];0])
    end
    r_sim = [r_sim[:,1]; r_sim[:,2]; r_sim[:,3]]
    Λ_sim = [Λ_sim[:,1]; Λ_sim[:,2]; Λ_sim[:,3]]
    # Simulation Moments
    m_sim = [pop_sim; r_sim; Λ_sim]
    # Assemble data moments
    pop_data = [data.CshareH19 + data.CshareH19_C; data.XshareH19; data.OshareH19; data.CshareL19 + data.CshareL19_C; data.XshareL19; data.OshareL19]
    r_data = [data.Crent19; data.Xrent19; data.mortgage19]
    Λ_data = [data.Cstock19; data.Xstock19; data.Ostock19]
    m_data = [pop_data; r_data; Λ_data]
    # Assemble weight matrix
    W = Diagonal(inv((m_data .- mean(m_data))*(m_data .- mean(m_data))'))
    # Calculate objective function value
    obj = (m_data - m_sim)'*W*(m_data - m_sim)
    return obj[1]
end

obj_SMM(ς,ν,ε,tol,N,ξ,lotto,varR,data,para_SMM)

# Function to execute SMM repeatedly with different initial conditions, compile averate results
function SMM(N,guess,times)
    ς_store = zeros(times)
    for t in 1:times
        ν = drawν(N)
        ε =  endowε(N,data,para_SMM)
        results = optimize(x -> obj_SMM(x,ν,ε,tol,N,ξ,lotto,varR,data,para_SMM),guess,Newton(),Optim.Options(iterations = 10000))
        ς_store[t] = Optim.minimizer(results)
    end # t loop
    ς = mean(ς_store)
    return ς
end
guess = [200000.0]

N = Int.(1000000 .*ones(2,2))
ς = 20000

ξ = ones(39,3,2)

SMM(N,guess,1)
