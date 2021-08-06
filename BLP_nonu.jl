# A program to run BLP methodology on sf_data to obtain parameter estimates for counterfactual simulation, this time with out idiosyncratic preferences

# ******Before running, use ] activate fresh, add Optim*******
using Pkg
using LinearAlgebra, CSV, DataFrames, Distributions, Parameters, Optim

@with_kw struct parameters
    θ::Float64 = -2
    ρ::Float64 = 1.1341810237625438
end

para = parameters()

cd("/Users/bigfoot13770/Documents/UO ECON PROGRAM/FUTILITY/JMP_RCsf")
# Read in sf data

#------------------------
# Keep λ's between 0 and 1 using a logit constraint
function constrain1(x; min=0, max=1)
    (max - min) * (1/(1 + exp(-x))) + min
end


# DECLARE BLP FUNCTIONS FOR USE BELOW
# endow -- Give agents initial states and sets of prices
function endow_noν(data, R, N)

        # Set up price matrix
        rents = repeat([R[:,1] ; R[:,3] ; R[:,4]]',N[2])
        # Loop by tract and endow some agents with RC_C by replacement
        index = 1
        poptrack = 0
        for x in 1:size(data.CshareH10)[1]
                # Find number of RCs to replace
                pop = round(Int,N[1]*data.CshareH10[x])
                # Reassign values
                rents[index:(pop + poptrack),x] .= R[x,2]
                # Update indexes
                index += pop
                poptrack += pop
        end
        # Assemble output
        consum = rents

        return consum
end


# Berry Contraction loop
function innerloop(alphasigma, μ, tol, iter, skill, data)
        N = Int(size(μ)[1])
        J = Int(size(μ)[2]/3)

        # Declare lambda's
        λC = constrain1(alphasigma[2])
        λX = constrain1(alphasigma[3])
        λO = constrain1(alphasigma[4])
        # Initialize delta loop
        δC = ones(J)
        δX = ones(J)
        δO = ones(J)
        δ = 0
        count = Int(1)
        CONT = true
        while CONT == true && count < 200
                # Calculate model shares based on μ, etc.
                denom = (sum(exp.((δC' .+ reshape(μ[1:3:end],N,J))./λC),dims = 2)).^λC + (sum(exp.((δX' .+ reshape(μ[2:3:end],N,J)./λX)),dims = 2)).^λX  + (sum(exp.((δO' .+ reshape(μ[1:3:end],N,J)./λO)),dims = 2)).^λO .+ 1
                # Calculate model predicted shares
                # Note: since λ is different for each housing type, have to do this separately then compile into a single δ later
                sharem_C = sum(exp.((δC' .+ reshape(μ[1:3:end],N,J))/λC) .*(sum(exp.((δC' .+reshape(μ[1:3:end],N,J))/λC)))^(λC-1)./denom, dims = 1)./N
                sharem_X = sum(exp.((δX' .+ reshape(μ[2:3:end],N,J))/λX) .*(sum(exp.((δX' .+reshape(μ[2:3:end],N,J))/λX)))^(λX-1)./denom, dims = 1)./N
                sharem_O = sum(exp.((δO' .+ reshape(μ[3:3:end],N,J))/λO) .*(sum(exp.((δO' .+reshape(μ[3:3:end],N,J))/λO)))^(λO-1)./denom, dims = 1)./N

                # Delta contraction
                # Condition on skill level to get right moments from data
                δCnew = 0
                δXnew = 0
                δOnew = 0
                if skill == "H"
                        δCnew = δC + λC.*(log.(data.CshareH19 + data.CshareH19_C) - log.(sharem_C'))
                        δXnew = δX + λX.*(log.(data.XshareH19) - log.(sharem_X'))
                        δOnew = δO + λO.*(log.(data.OshareH19) - log.(sharem_O'))

                else
                        δCnew = δC + λC.*(log.(data.CshareL19 + data.CshareL19_C) - log.(sharem_C'))
                        δXnew = δX + λX.*(log.(data.XshareL19) - log.(sharem_X'))
                        δOnew = δO + λO.*(log.(data.OshareL19) - log.(sharem_O'))
                end
                # Check for convergence
                δ = [δC ; δX ; δO]
                δnew = [δCnew ; δXnew ; δOnew]
                if findmax(abs.(δ - δnew))[1] < tol || count > iter
                        δ = δnew
                        CONT = false
                elseif isnan(findmax(abs.(δ - δnew))[1]) == true
                        δ = δnew
                        println("Oops")
                        CONT = false

                else
                        println(findmax(abs.(δ - δnew))[1])
                        δC = δCnew
                        δX = δXnew
                        δO = δOnew
                        count += 1
                end
        end


        return δ
end

# Function to fminsearch over
function outerloop_noν(alphasigma, consum, R, skill, iter, data)
        N = size(data)[1]*3
        # set up μ
        μ = -alphasigma[1].*consum[:,1:size(consum)[2]]
        # find δ
        δ = innerloop(alphasigma, μ, 0.01, iter, skill, data)
        # Format data matrices
        X = []
        if skill =="H"
                X = repeat([data.cbd./1000 data.Hcommute./1000 data.parks data.park_share data.water./1000], outer = 3)
        else
                X = repeat([data.cbd./1000 data.Lcommute./1000 data.parks data.park_share data.water./1000], outer = 3)
       end
        # Z, matrix of zoning instruments
        Z = repeat([data.res_share data.com_share data.ind_share data.mix_share], outer = 3)
        W = inv(Z'Z)
        #
        β = inv(transpose(X)*Z*inv(W)*transpose(Z)*X)*transpose(X)*Z*W*transpose(Z)*δ
        # Solve for unobservable amenities
        ξ = δ - X*β
        # return objective function
        obj = ξ'*Z*W*Z'*ξ
        return obj[1]
end


H = Int.([1286 1738].*2);
#H = [1286088 1738678];
#L = [1668244 1770388];

alphasigma = [0.5;0.5;0.5;0.5];
β = [1;1;1;1;1];

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

R = makeR(data,2019,100, 1000, para)
consum = endow_noν(data, R, H)

test = optimize(x -> outerloop_noν(x, consum, R, "H", 1000, data),alphasigma,NelderMead(),Optim.Options(iterations = 10000))


alphasigma_noν = Optim.minimizer(test)


function outerloop_noν_return(alphasigma, consum, R, skill, iter, data)
        N = size(data)[1]*3
        # set up μ
        μ = -alphasigma[1].*consum[:,1:size(consum)[2]]
        # find δ
        δ = innerloop(alphasigma, μ, 0.01, iter, skill, data)
        # Format data matrices
        X = []
        if skill =="H"
                X = repeat([data.cbd./1000 data.Hcommute./1000 data.parks data.park_share data.water./1000], outer = 3)
        else
                X = repeat([data.cbd./1000 data.Lcommute./1000 data.parks data.park_share data.water./1000], outer = 3)
       end
        # Z, matrix of zoning instruments
        Z = repeat([data.res_share data.com_share data.ind_share data.mix_share], outer = 3)
        W = inv(Z'Z./N)
        #
        β = inv(transpose(X)*Z*inv(W)*transpose(Z)*X)*transpose(X)*Z*W*transpose(Z)*δ
        # Solve for unobservable amenities
        ξ = δ - X*β
        # return objective function
        obj = ξ'*Z*W*Z'*ξ
        return β, δ, ξ
end


β, δ, ξ = outerloop_noν_return(alphasigma, consum, R, "H", 1000, data)


function makeG(Z,diffξ,N,K)
        G = zeros(size(Z)[2],K)
        for x in 1:N
                G = G + Z[x,:] * transpose(diffξ[x,:])
        end
        G = G./N
        return G
end

function se_noν(perturb, alphasigma, consum, R, skill, iter, data)
        N = 3*size(data)[1]
        K = size(alphasigma)[1]
        # Reassemble X, Z, W
        X = []
        if skill =="H"
                X = repeat([data.cbd./1000 data.Hcommute./1000 data.parks data.park_share data.water./1000], outer = 3)
        else
               X = repeat([data.cbd./1000 data.Lcommute./1000 data.parks data.park_share data.water./1000], outer = 3)
       end
        # Z, matrix of zoning instruments
        Z = repeat([data.res_share data.com_share data.ind_share data.mix_share], outer = 3)
        W = inv(Z'Z)
        # Obtain ξ
        shitty, poop, ξ = outerloop_noν_return(alphasigma, consum, R, skill, iter, data)

        diffξ = zeros(N,K)
        # Perturb alphasigma, observe change in ξ
        for k in 1:K
                turb = zeros(K)
                turb[k] = perturb
                pert = alphasigma + turb
                crap, δalt, ξalt = outerloop_noν_return(pert, consum, R, skill, iter, data)
                diffξ[:,k] = (ξalt-ξ)./perturb
        end

        G = makeG(Z,diffξ,N,K)

        # Assemble S
        S = var(ξ, corrected = false)/N .* Z'Z
        # Compute Σ
        Σ = inv(G'W*G)*G'W*S*W*G*inv(G'W*G)

        alphasigma_var = 1/N * sqrt.(diag(Σ))

        return alphasigma_var
end

se = se_noν(1e-10, alphasigma_noν, consum, R, "H", 1000, data)
