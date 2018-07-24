L = 4
tau = 1/2
n = 10

function makeD()
    d = ones(L+1)
    for j = 2:L+1
      d[j] = d[j-1]*(L-j+2)*(L-j+2-tau)/(j-1)/(j-1+tau)
    end
    @show(d)
    return(d)
end

function makeQ(d)
    d2 = copy(d)
    sq = [j^2 for j = 0:length(d)-1]
    a = ones[n]
    for j = 2:n
        a[j] = -1 * a[j-1] * tau^2 / (2*j) / (2*j-1)
    end
    b = ones(n)
    for j = 1:n
        b[j] = sum(d2) * (-1)^(j-1) / factorial(2*j-2)
        d2 = d2 .* sq
    end
    return(conv(a,b))
end

function makeP(d)
    d2 = copy(d)
    s = [j+tau for j = 0:length(d)-1]
    d2 = d2 .* s
    p = zeros(n)
    for j = 1:n
        p[j] = sum(d2) * (-1)^(j-1) / factorial(2*j-1)
        d2 = d2 .* s
        d2 = d2 .* s
    end
    return(p)
end

function makeR()
    d = makeD()
    p = makeP(d)
    q = makeQ(d)
    r = zeros(n)
end
