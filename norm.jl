tau = .5
function makeD(L)
    d = ones(L+1)
    for j = 2:L+1
      d[j] = d[j-1]*(L-j+2)*(L-j+2-tau)/(j-1)/(j-1+tau)
    end
    return(d)
end

function testPi()
    n = 20
    res = zeros(n)
    for j = 1:n
        d = makeD(j)
        p = [(-1)^k for k = 1:length(d)]
        res[j] = abs.(d'*p)
    end
    return(res)
end
