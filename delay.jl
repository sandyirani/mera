

function approx1(L,tau)
    d = zeros(L+1)
    d[1] = 1
    for j=1:L
        d[j+1] = d[j]*(L-j+1)*(tau+j-1)/j/(tau+L+j)
    end
    return(d)
end

function apply(d,z)
    sum = d[1]
    pow = 1
    for j=2:length(d)
        pow = pow*z
        sum = sum + d[j]*pow
    end
    return(sum)
end

function test1(L)
    tau = 1/2
    vals = [1,1.1,1.2,1.3]
    d = approx1(L,tau)
    for j = 1:length(vals)
        out = apply(d,vals[j])/apply(d,1/vals[j])
        diff = out - (1/(vals[j]^tau))
        @show(diff)
    end
end

function approx2(L,tau)
    d = ones(L+1)
    for j = 2:L+1
      d[j] = d[j-1]*(L-j+2)*(L-j+2-tau)/(j-1)/(j-1+tau)
    end
    @show(d)
    return(d)
end

function test2(L)
    tau = 1/2
    vals = [.7,.8,.9,1,1.1,1.2,1.3,1.9]
    d = approx2(L,tau)
    d2 = [d[length(d)-j+1] for j = 1:length(d)]
    for j = 1:length(vals)
        out = apply(d,vals[j])/apply(d2,vals[j])

        out = apply(d2,vals[j])
        #out = apply(d,vals[j])
        diff = out - (1/(vals[j]^tau))
        @show(out)
    end
end
