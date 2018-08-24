tau = .5
function makeD(L)
    d = ones(L+1)
    for j = 2:L+1
      d[j] = d[j-1]*(L-j+2)*(L-j+2-tau)/(j-1)/(j-1+tau)
    end
    return(d)
end

function makeR(L)
    h = Int8(floor((L+1)/2))
    r = ones(h-1)
    for j = 1:h-1
      r[j] = r[j]*(1 - 4*j/(L+2*j+1))
      r[j] = r[j]*(1 - 4*j/(L+2*j+3))
      r[j] = r[j]*(1 - 2/(L+2*j+2))
      #r[j] = (1 - (8*j+2)/(L+2*j+1))
      r[j] = r[j]*(1 + 2/(2*j-1))
    end
    return(r)
end

function makeDsimple(L)
    d = ones(L+1)
    for j = 2:L+1
      d[j] = d[j-1]*(L-j+2)*(L-j+2)/(j-1)/(j-1)
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
  r2 = [b[j]/c[j-1] for j = 2:length(b)]
  #@show(change/h)
  @show(c[length(c)]^2,c[1])
  return(b,c,r,r2)
end

function testR(L)
    r = makeR(L)
    len = Int8(floor((length(r)+1)/2))
    res = zeros(len)
    res2 = zeros(len-1)
    res[1] = 1 - 1/r[1]
    res[2] = r[2]*(1-r[3])
    cum = r[2]
    res2[1] = res[2]
    for j = 4:2:length(r)-1
        cum = cum*r[j-1]*r[j]
        idx = Int8(floor(j/2))
        res[idx+1] = res[idx] + cum*(1-r[j+1])
        res2[idx] = cum*(1-r[j+1])
    end
    return(r,res,res2)
end
