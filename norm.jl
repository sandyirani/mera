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

function test2(L)
  d = makeD(L)
  p = [(-1)^(k-1) for k = 1:length(d)]
  a = d .* p
  h = Int8(floor((L+1)/2))
  b = [a[j] + a[L+2-j] for j = h:-1:1]
  r = [b[j]/b[j-1] for j = 2:length(b)]
  c = zeros(length(b))
  c[1] = b[1]
  for j = 2:length(c)
    c[j] = c[j-1] + b[j]
  end
  done = false
  change = 0
  for j = 2:length(c)
      if (!done && c[j]*c[j-1] > 0)
          change = j
          done = true
      end
  end
  #@show(change/h)
  @show(c[length(c)]^2,c[1])
  return(b,c,r)
end
