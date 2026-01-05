# This code only defines the functions necessary to run the evolutionary simulations.
# You can run this code in the "run.ipynb" Jupyter notebook for any parameter set desired.

function rA(n, parameters)
    return n <= parameters.N_U ? parameters.r_o : (parameters.N_U * parameters.r_o + (n - parameters.N_U) * parameters.r_f)/n
end

function rB(n, parameters)
    return n <= parameters.N_L ? parameters.r_o : (parameters.N_L * parameters.r_o + (n - parameters.N_L) * parameters.r_f)/n
end

function pA(parameters)
    return [i >= parameters.N_S ? 0.0 : 1.0 - parameters.mixing  for i in 1:parameters.N_max ]
end

function pB(parameters)
    return ones(parameters.N_max) - pA(parameters)
end

function r(parameters)
    return [pA(parameters)[n] * rA(n, parameters) + pB(parameters)[n] * rB(n, parameters) for n in 1:parameters.N_max]
end

# We consider three life cycles: N+1, Nx1, N+N. The life cycles are implemented with the same fragmentation size. 

function dynamicalModel!(u, p, t)
    LCtype = p.LCtype
    LC = p.lifeCyclesPresent
    m = length(LC)
    m_max = maximum(LC)
    # u will be a matrix with m rows and m_max columns

    #competition 
    T_A = transpose(ones(m)) * u * (pA(p)[1:m_max] .* (1:m_max))
    T_B = transpose(ones(m)) * u * (pB(p)[1:m_max] .* (1:m_max))

    # reproduction
    births = u .* (ones(m) * transpose(r(p)[1:m_max] .* (1:m_max)))
    gain = zeros(size(u))
    gain[:,2:m_max] = births[:,1:(m_max - 1)]
    if LCtype == "N+1"
        gain[:,1] += [births[i,LC[i]] for i in 1:m]
        for i in 1:m
            births[i,LC[i]] = 0.0
        end
    elseif LCtype == "Nx1"
        gain[:,1] += [(1+LC[i]) * births[i,LC[i]] for i in 1:m]
    elseif LCtype == "N+N"
        for i in 1:m
            gain[i,Int(floor(0.5*(1+LC[i])))] += births[i,LC[i]]
            gain[i,Int(ceil(0.5*(1+LC[i])))] += births[i,LC[i]]
        end
    end
    for i in 1:m
        for j in 1:m_max
            if j > LC[i]
                gain[i,j] = 0
                births[i,j] = 0
            end
        end
    end
    du = gain - births - p.gamma * T_A * (ones(m) * transpose(pA(p)[1:m_max])).* u - p.gamma * T_B * (ones(m) * transpose(pB(p)[1:m_max])) .* u 
end

logR0 = function(x,n,parameters)
    TA = sum(x * ((1:parameters.N_max) .* pA(parameters)))
    TB = sum(x * ((1:parameters.N_max) .* pB(parameters)))
    r = function(m)
        pA(parameters)[m] * rA(m, parameters) + pB(parameters)[m] * rB(m, parameters)
    end
    c = function(m)
        parameters.gamma * (pA(parameters)[m] * TA + pB(parameters)[m] * TB)
    end
    if n == 1
        logR = log(r(1)/c(1))
    elseif n > 1
        if parameters.LCtype == "N+1"
            logR = log(n * r(n) / c(n))
            lower = 1
            upper = n-1
        elseif parameters.LCtype == "Nx1"
            logR = log(n+1)
            lower = 1
            upper = n
        elseif parameters.LCtype == "N+N"
            if iseven(n+1)
                logR = log(2)
                lower = Int(round((n+1)/2))
            else
                i = Int(round(n/2)) 
                lower = Int(round((n/2)+1))
                logR = log(1.0+(i * r(i))/ (i*r(i) + c(i)))
            end
            upper = n
        end
        for i in lower:upper
            logR = logR + log((i * r(i))/ (i*r(i) + c(i)))
        end
    end
    return(logR)
end

evolutionarySimulation = function(parameters, invasionAttempts = 1000, R0cutoff = 1e-8)
    Random.seed!(parameters.seed)
    N_max = parameters.N_max
    x = zeros(N_max,N_max)
    x[1,1] = r(parameters)[1]/(parameters.gamma * (pA(parameters)[1] .^2 + pB(parameters)[1] .^2))
    lifeCyclesPresent = [1]

    for z in 1:invasionAttempts
        println("The life cycles present are " * string(lifeCyclesPresent))
        logR0s = [logR0(x,i,parameters) for i in 1:N_max]
        
        if maximum(logR0s) <= R0cutoff
            neutralInvaders = filter(i -> logR0(x, i, parameters) >= -R0cutoff, setdiff(1:N_max, lifeCyclesPresent))
            lifeCyclesPresent = sort([lifeCyclesPresent; neutralInvaders])
            if length(neutralInvaders) > 0
                println("The life cycles " * string(neutralInvaders) * " were allowed to invade neutrally!")
            end
        end

        potentialInvaders = filter(i->logR0(x,i,parameters) > R0cutoff, setdiff(1:N_max, lifeCyclesPresent))

        if length(potentialInvaders) > 0
            # select invader
            invader = rand(potentialInvaders) 
            println(string(invader) * " is invading with an R0 of " * string(exp(logR0(x,invader,parameters))) * " ...")

            # sprinkle in invader with 1% of the total cell population
            x[invader, 1] = 0.01 * sum(x * (1:N_max))
            lifeCyclesPresent = sort([lifeCyclesPresent; invader])
            
            # run simulation 
            names = (:N_max, :gamma, :N_U, :N_L, :N_S, :r_o, :r_f, :mixing, :LCtype, :lifeCyclesPresent)
            values = (parameters.N_max, parameters.gamma, parameters.N_U, parameters.N_L, parameters.N_S, parameters.r_o, parameters.r_f, parameters.mixing, parameters.LCtype, lifeCyclesPresent)
            p = NamedTuple{names}(values)
            u0 = x[lifeCyclesPresent,1:maximum(lifeCyclesPresent)]
            prob = ODEProblem(dynamicalModel!, u0, (0,5e7), p)
            sol =  solve(prob, lsoda(), reltol = 1e-12, abstol = 1e-12)
            u = sol.u[length(sol.t)]

            x = zeros(N_max,N_max)
            for i in 1:(length(lifeCyclesPresent))
                x[lifeCyclesPresent[i],1:(lifeCyclesPresent[i])] = u[i,1:(lifeCyclesPresent[i])]
            end

            # exclude life cycles
            excludedLifeCycles = filter(i->((logR0(x,i,parameters) < - R0cutoff) && (sum((x * (1:N_max))[i]) < 1e-6)), 1:N_max)
            for n in excludedLifeCycles
                for i in 1:n
                    x[n,i] = 0.0
                end
            end
            x = tril(x) 
            lifeCyclesPresent = sort(setdiff(lifeCyclesPresent, excludedLifeCycles))
        else
            println("The life cycles present are " * string(lifeCyclesPresent))
            println("ESC reached!")
            return lifeCyclesPresent
        end

    end
end