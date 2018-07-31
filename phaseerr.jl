L = 2
tau = 1/2
taup = -(L-tau)/2
n = 20

function makeD()
    d = ones(L+1)
    for j = 2:L+1
      d[j] = d[j-1]*(L-j+2)*(L-j+2-tau)/(j-1)/(j-1+tau)
    end
    @show(d)
    return(d)
end

function testD(d)
    d2 = copy(d)
    s = [j+taup for j = 0:length(d)-1]
    d2 = d2 .* s
    p = zeros(n)
    for j = 1:n
        p[j] = sum(d2)
        d2 = d2 .* s
        d2 = d2 .* s
    end
    return(p)
end

function makeQ(d)
    d2 = copy(d)
    sq = [j^2 for j = 0:length(d)-1]
    a = ones(n)
    for j = 2:n
        a[j] = -1 * a[j-1] * taup^2 / (2*j) / (2*j-1)
    end
    b = ones(n)
    for j = 0:n-1
        b[j+1] = sum(d2) * (-1)^(j)
        d2 = d2 .* sq
        d2 = d2/(j+1)/(j+2)
    end
    return(conv(a,b))
end

function makeP(d)
    d2 = copy(d)
    s = [j+taup for j = 0:length(d)-1]
    d2 = d2 .* s
    p = zeros(n)
    for k = 0:n-1
      for j = 0:L
        c = 1
        for m = 2*k+1:-1:1
          c = c * (j+taup)/m
        end
        p[k+1] = p[k+1] + d[j+1] * c
      end
      p[k+1] = p[k+1] * (-1)^k
    end
    return(p)
end

function makeR()
    d = makeD()
    p = makeP(d)
    q = makeQ(d)
    r = zeros(n)
    for k = 0:n-1
        r[k+1] = p[k+1]
        for j = 1:k
            r[k+1] = r[k+1] - (r[k-j+1]*q[j+1])
        end
        r[k+1] = r[k+1]/q[1]
    end
    return(r)
end
