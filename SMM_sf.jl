using Pkg
using LinearAlgebra, CSV, DataFrames, Distributions, Parameters

@with_kw struct parameters_SMM
    θ::Float64 = -2
    α::Array = [.13 .13]
    αhat::Array = [0.33 0.33]
    β::Array = [1 1 1 1 1; 1 1 1 1 1]
    σ::Array = [1 1 1 1 1; 1 1 1 1 1]
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
                        dC = truncated(Normal(data.Crent10[j],data.CsigmaR10[j]),data.Crent10[j] - varR,data.Crent10[j] + varR);
                        dC_C = truncated(Normal(data.Crent10[j],data.CsigmaR10[j]),data.Crent10[j] - varR,data.Crent10[j] + varR);
                        dX = truncated(Normal(data.Xrent10[j],data.XsigmaR10[j]),data.Xrent10[j] - varR,data.Xrent10[j] + varR);
                        dO = truncated(Normal(data.mortgage10[j],max(data.Mortgagesigma10[j],1e-3)),data.mortgage10[j] - varR,data.mortgage10[j] + varR);
                            RC[j] = sum((rand(dC,data.Cstruct10[j]).^(1+θ) ) )^(1/(1+θ));
                            RC_C[j] = ρ*sum((rand(dC_C,data.Cstruct10[j]).^(1+θ) ) )^(1/(1+θ));
                            RX[j] = sum((rand(dX,data.Xstruct10[j]).^(1+θ) ) )^(1/(1+θ));
                            RO[j] = sum((rand(dO,data.Ostruct10[j]).^(1+θ) ) )^(1/(1+θ));
                    end
                else

                    for j in 1:size(data)[1]
                        dC = truncated(Normal(data.Crent19[j],data.CsigmaR10[j]),data.Crent19[j] - varR,data.Crent19[j] + varR);
                        dC_C = truncated(Normal(data.Crent19[j],data.CsigmaR10[j]),data.Crent19[j] - varR,data.Crent19[j] + varR);
                        dX = truncated(Normal(data.Xrent19[j],data.XsigmaR10[j]),data.Xrent19[j] - varR,data.Xrent19[j] + varR);
                        dO = truncated(Normal(data.mortgage19[j],max(data.Mortgagesigma19[j],1e-3)),data.mortgage19[j] - varR,data.mortgage19[j] + varR);
                            RC[j] = sum((rand(dC,data.Cstruct19[j]).^(1+θ) ) )^(1/(1+θ));
                            RC_C[j] = ρ*sum((rand(dC_C,data.Cstruct10[j]).^(1+θ) ) )^(1/(1+θ));
                            RX[j] = sum((rand(dX,data.Xstruct19[j]).^(1+θ) ) )^(1/(1+θ));
                            RO[j] = sum((rand(dO,data.Ostruct19[j]).^(1+θ) ) )^(1/(1+θ));
                    end

            end
            R[:,:,t] = [RC RC_C RX RO]
        end

        R = mean(R,dims = 3)[:,:,1]

        return R
end


function endow(data, R, N)
        # ν draws:
        ν = randn(N[2],5)

        # Set up price matrix
        rents = repeat([R.RC ; R.RX ; R.RO]',N[2])
        # Loop by tract and endow some agents with RC_C by replacement
        index = 1
        poptrack = 0
        for x in 1:size(data.CshareH10)[1]
                # Find number of RCs to replace
                pop = round(Int,N[1]*data.CshareH10[x])

                # Reassign values
                rents[index:(pop + poptrack),x] .= R.RC_C[x]

                # Update indexes
                index += pop
                poptrack += pop
        end
        # Assemble output
        consum = [ν rents]

        return consum
end

# Function for calculating household shares
function hhsort(data,consumH,consumL,R,para_SMM)
    @unpack θ, α, βH, βL, σH, σL, λH, λL = para_SMM
    # Calculate δ's conditional on endowments
    XH = repeat([data.cbd data.Hcommute data.parks data.park_share data.water], outer = 3)
    XL = repeat([data.cbd data.Lcommute data.parks data.park_share data.water], outer = 3)

end

lotto = 600
varR = 200

# Function to endow landowners with initial conditions
function lo_endow(data, lotto, varR)
    L0 = []

    for j in 1:size(data)[1]

        # Set up distributions to draw from for initial conditions
        # Structure size Λ0
        dC = truncated(Normal(data.Csize10[j],max(data.Csigma10[j],.001)),1, Inf);
        dC_C = truncated(Normal(data.Csize10[j],max(data.Csigma10[j],.001)),2, 6)
        dX = truncated(Normal(data.Xsize10[j],max(data.Xsigma10[j],.001)), 1, Inf);
        dO = truncated(Normal(data.Osize10[j],max(data.Osigma10[j],.001)),1,Inf);
        # Rents r0
        dCR = truncated(Normal(data.Crent10[j],data.CsigmaR10[j]),data.Crent10[j] - varR,data.Crent10[j] + varR);
        dCR_C = truncated(Normal(data.Crent10[j],data.CsigmaR10[j]),data.Crent10[j] - varR,data.Crent10[j] + varR);
        dXR = truncated(Normal(data.Xrent10[j],data.XsigmaR10[j]),data.Xrent10[j] - varR,data.Xrent10[j] + varR);
        dOR = truncated(Normal(data.mortgage10[j],max(data.Mortgagesigma10[j],1e-3)),data.mortgage10[j] - varR,data.mortgage10[j] + varR);
        # create L0η matricies to store initial conditions by η
        L0C = [rand(dCR,max(Int(round(data.Cstruct10[j]-lotto*data.win[j])),0)) rand(dC,max(Int(round(data.Cstruct10[j] - lotto*data.win[j])),0))]
        if data.Cstruct10[j]-lotto*data.win[j] < 0
            L0C_C = [0 0]
        else
        L0C_C = [rand(dCR,Int(round(lotto*data.win[j]))) rand(dC,Int(round(lotto*data.win[j])))]
        end
        L0X = [rand(dXR,data.Xstruct10[j]) rand(dX,data.Xstruct10[j])]
        L0O = [rand(dOR,data.Ostruct10[j]) rand(dO,data.Ostruct10[j])]
        # Store in L0
        push!(L0, [[L0C] [L0C_C] [L0X] [L0O]])
    end

    return L0
end



# Function for characterizing landowner reactions
function losort(L0, M, R0, movers, data, ς, Vmin, para_SMM)
    # Set up parameters
    @unpack θ,αhat,y,ρ,γ,ϕ,ψ,κ = para_SMM
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
                φ = movers[j,1]*αhat[1]*y[1]/(R0[j,η].^(1+θ)).*r0.^θ + movers[j,2]*αhat[2]*y[2]/(R0[j,η].^(1+θ)).*r0.^θ
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
                Stasis = ρ.*r0.*Λ0
                Adjustment = maxi[1][:,1]
                # Set up renovation to η = X/0 = 3/4
                RenoX = (Psi[j]*M[j,3]^(ψ-1))^(ϑ*(θ+1))*M[j,3] .- ϕ.*Λ0 .- (Psi[j]*(1+θ)/(ψ*θ)).*(Psi[j]*M[j,4]^(ψ-1))^(ϑ*ψ) .- data.avgBO[j].*Λ0
                RenoO = (Psi[j]*M[j,4]^(ψ-1))^(ϑ*(θ+1))*M[j,4] .- ϕ.*Λ0 .- (Psi[j]*(1+θ)/(ψ*θ)).*(Psi[j]*M[j,4]^(ψ-1))^(ϑ*ψ) .- data.avgBO[j].*Λ0
                profits = [Stasis+u[:,1] Adjustment+u[:,1] RenoX+u[:,2] RenoO+u[:,3]]
                # Assmeble rents
                rents = [ρ.*r0 rstar repeat([(Psi[j]*M[j,3]^(ψ-1))^ϑ], size(Λ0)[1]) repeat([(Psi[j]*M[j,4]^(ψ-1))^ϑ], size(Λ0)[1])]
                Λ = [Λ0 φ.*Λ0+M[j,1].*rstar.^θ M[j,3].*rents[:,3].^θ  M[j,4].*rents[:,4].^θ]
                # use findmax() to figure out indicies for optimal housing choice, iterate through cartesian indicies cause Julia sucks sometimes to assemble Λ and r storage matricies
                index = findmax(profits,dims = 2)[2][:]
                ΛC = zeros(size(profits)[1],1)
                rC = zeros(size(profits)[1],1)
                for k in 1:size(profits)[1]
                    ΛC[k] = Λ[k,Int(index[k][2])]
                    rC[k] = rents[k,Int(index[k][2])]
                end # k loop
            elseif η == 2 # Controlled, lotto win C⋆
                # Draw Gumbel shocks
                u = rand(d,size(L0[j][η])[1],3)
                Λ0 = L0[j][η][:,2]
                r0 = L0[j][η][:,1]
                # Calculate φ
                φ = movers[j,1]*αhat[1]*y[1]/R0[j,η].*r0.^θ + movers[j,2]*αhat[2]*y[2]/R0[j,η].*r0.^θ
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
                Stasis = ρ.*r0.*Λ0
                Adjustment = maxi[1][:,1]
                # Set up renovation to η = X/0 = 3/4
                Condo = (Λ0./M[j,4]).^((θ+1)/θ) .- (κ + data.avgBO[j]).*Λ0
                RenoX = (Psi[j]*M[j,3]^(ψ-1))^(ϑ*(θ+1))*M[j,3] .- ϕ.*Λ0 .- (Psi[j]*(1+θ)/(ψ*θ)).*(Psi[j]*M[j,4]^(ψ-1))^(ϑ*ψ) .- data.avgBO[j].*Λ0
                RenoO = (Psi[j]*M[j,4]^(ψ-1))^(ϑ*(θ+1))*M[j,4] .- ϕ.*Λ0 .- (Psi[j]*(1+θ)/(ψ*θ)).*(Psi[j]*M[j,4]^(ψ-1))^(ϑ*ψ) .- data.avgBO[j].*Λ0
                profits = [Stasis+u[:,1] Adjustment+u[:,1] RenoX+u[:,2] Condo+u[:,3] RenoO+u[:,3]]
                # Assmeble rents
                rents = [ρ.*r0 rstar repeat([(Psi[j]*M[j,3]^(ψ-1))^ϑ], size(Λ0)[1]) (Λ0./M[j,4]).^(1/θ) repeat([(Psi[j]*M[j,4]^(ψ-1))^ϑ], size(Λ0)[1])]
                Λ = [Λ0 φ.*Λ0+M[j,1].*rstar.^θ M[j,3].*rents[:,3].^θ Λ0 M[j,4].*rents[:,5].^θ]
                # use findmax() to figure out indicies for optimal housing choice, iterate through cartesian indicies cause Julia sucks sometimes to assemble Λ and r storage matricies
                index = findmax(profits,dims = 2)[2][:]
                ΛC_C = zeros(size(profits)[1],1)
                rC_C = zeros(size(profits)[1],1)
                for k in 1:size(profits)[1]
                    ΛC_C[k] = Λ[k,Int(index[k][2])]
                    rC_C[k] = rents[k,Int(index[k][2])]
                end # k loop
            elseif η == 3 # Exempt X
                # Draw Gumbel shocks
                u = rand(d,size(L0[j][η])[1],3)
                # Caclulate profits under each alternative
                Λ0 = L0[j][η][:,2]
                Stasis = ρ.*L0[j][η][:,1].*Λ0
                Stasis = (Λ0./M[j,η]).^((θ+1)/θ).*M[j,η]
                 Adjustment = (M[j,η]*Γ[j]./Λ0.^γ).^((θ+1)/(1-γ*θ))*M[j,η] .- (Γ[j]*(1+θ)/(γ*θ)).* ((M[j,η]*Γ[j]./Λ0.^γ).^(1/(1-γ*θ))./Λ0).^γ .- data.avgBO[j].*Λ0
                # Set up renovation to η = O = 4
                Condo = (Λ0./M[j,4]).^((θ+1)/θ) .- (κ + data.avgBO[j]).*Λ0
                Reno = (Psi[j]*M[j,4]^(ψ-1))^(ϑ*(θ+1))*M[j,4] .- ϕ.*Λ0 .- (Psi[j]*(1+θ)/(ψ*θ)).*(Psi[j]*M[j,4]^(ψ-1))^(ϑ*ψ) .- data.avgBO[j].*Λ0
                profits = [Stasis+u[:,1] Adjustment+u[:,1] Condo+u[:,2] Reno+u[:,2]]
                rents = [(Λ0./M[j,η]).^(1/θ) (M[j,η]*Γ[j]./Λ0.^γ).^(1/(1-γ*θ)) (Λ0./M[j,4]).^(1/θ) repeat([(Psi[j]*M[j,4]^(ψ-1))^ϑ], size(Λ0)[1])]
                Λ = [rents[:,1:2].^θ.*M[j,η] rents[:,3:4].^θ.*M[j,4]]
                # use findmax() to figure out indicies for optimal housing choice, iterate through cartesian indicies cause Julia sucks sometimes to assemble Λ and r storage matricies
                index = findmax(profits,dims = 2)[2][:]
                ΛX = zeros(size(profits)[1],1)
                rX = zeros(size(profits)[1],1)
                for k in 1:size(profits)[1]
                    ΛX[k] = Λ[k,Int(index[k][2])]
                    rX[k] = rents[k,Int(index[k][2])]
                end # k loop

            else # Owned O
                # Draw Gumbel shocks
                u = rand(d,size(L0[j][η])[1],2)
                # Caclulate profits under each alternative
                Λ0 = L0[j][η][:,2]
                Stasis = (Λ0./M[j,η]).^((θ+1)/θ).*M[j,η]
                Adjustment = (M[j,η]*Γ[j]./Λ0.^γ).^((θ+1)/(1-γ*θ))*M[j,η] .- (Γ[j]*(1+θ)/(γ*θ)).* ((M[j,η]*Γ[j]./Λ0.^γ).^(1/(1-γ*θ))./Λ0).^γ .- data.avgBO[j].*Λ0
                # Set up renovation to η = X = 3
                Reno = (Psi[j]*M[j,3]^(ψ-1))^(ϑ*(θ+1))*M[j,3] .- ϕ.*Λ0 .- (Psi[j]*(1+θ)/(ψ*θ)).*(Psi[j]*M[j,3]^(ψ-1))^(ϑ*ψ) .- data.avgBO[j].*Λ0
                profits = [Stasis+u[:,1] Adjustment+u[:,1] Reno+u[:,2]]
                rents = [(Λ0./M[j,η]).^(1/θ) (M[j,η]*Γ[j]./Λ0.^γ).^(1/(1-γ*θ)) repeat([(Psi[j]*M[j,3]^(ψ-1))^ϑ], size(Λ0)[1])]
                Λ = [rents[:,1:2].^θ.*M[j,η] rents[:,3].^θ.*M[j,3]]
                # use findmax() to figure out indicies for optimal housing choice, iterate through cartesian indicies cause Julia sucks sometimes to assemble Λ and r storage matricies
                index = findmax(profits,dims = 2)[2][:]
                ΛO = zeros(size(profits)[1],1)
                rO = zeros(size(profits)[1],1)
                for k in 1:size(profits)[1]
                    ΛO[k] = Λ[k,Int(index[k][2])]
                    rO[k] = rents[k,Int(index[k][2])]
                end # k loop
            end # η check
        end # η loop
        # Calculate πj(X/O) for new builders
        πX = (Psi[j]*M[j,3]^(ψ-1))^(ϑ*(θ+1))*M[j,3] .- (Psi[j]*(1+θ)/(ψ*θ)).*(Psi[j]*M[j,3]^(ψ-1))^(ϑ*ψ)
        πO = (Psi[j]*M[j,4]^(ψ-1))^(ϑ*(θ+1))*M[j,4] .- (Psi[j]*(1+θ)/(ψ*θ)).*(Psi[j]*M[j,4]^(ψ-1))^(ϑ*ψ)

        # Calculate πj(V)
        πV[j] = log(Vbar[j]/Vj[j]*(exp(πX) + exp(πO))/max(1 - Vbar[j]/Vj[j],0.0000001))
        # Push results into bf matricies
        push!(𝚲,[[ΛC] [ΛC_C] [ΛX] [ΛO]])
        push!(𝐫,[[rC] [rC_C] [rX] [rO]])
    end # j loop
    return 𝐫, 𝚲, πV
end # losort

data = CSV.read("sf_data.csv",DataFrame)

R0 = makeR(data,2019,100,1000,para_SMM)
L0 = lo_endow(data, lotto, varR)
movers = 10 .*ones(size(data)[1],2)
M = 1000000 .*ones(size(data)[1],4)
ς = 2

testr, testΛ, testπV = losort(L0, M, R0, movers, data, ς, 100, para_SMM)

function assembleR(r,L0,para_SMM)
    @unpack θ, ρ = para_SMM
    RC = zeros(size(r)[1])
    RC_C = zeros(size(r)[1])
    RX = zeros(size(r)[1])
    RO = zeros(size(r)[1])
    # Loop by j
    for j in 1:size(r)[1]
        RC[j] = (sum((r[j][1]).^(θ+1)))^(1/(θ+1))
        RC_C[j] = ρ*(sum((L0[j][2][:,1]).^(θ+1)))^(1/(θ+1))
        RX[j] = (sum((r[j][3]).^(θ+1)))^(1/(θ+1))
        RO[j] = (sum((r[j][4]).^(θ+1)))^(1/(θ+1))
    end # j loop
    R = [RC RC_C RX RO]
end

assembleR(testr,L0,para_SMM)




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

R = 10 .*ones(2,4)
R0 = R
L0 = []
push!(L0,[[repeat([1000 5],5000)] [repeat([1000 5],5000)] [repeat([1000 5],5000)] [repeat([1000 5],5000)] ])

testr, testΛ, testπV = losort(L0, M, R0, movers, fata, ς, 100, para_SMM)

function SMM()

end
