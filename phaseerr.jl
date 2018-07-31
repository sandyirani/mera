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

function genAngles(a,b,m)
  startIdx = Int32(ceil(a*m))
  finalIdx = Int32(ceil(b*m)-1)
  numIdx = (finalIdx - startIdx + 1)
  angles = [pi*j/m for j = 0:numIdx-1]
end

function phaseErr(a,b,m)
     r = makeR()
     ang = genAngles(a,b,m)
     err = zeros(length(ang))
     for j = 1:length(ang)
       v = [(ang[j])^k for k = 0:length(r)-1]
       err[j] = v'*r
     end
     err
end

function phase(a,b,m)
  ang = genAngles(a,b,m)
  p = zeros(length(ang))
  p2 = zeros(length(ang))
  diff = zeros(length(ang))
  d = makeD()
  for j = 1:length(ang)
    s = [sin(ang[j]*k) for k = 0:length(d)-1]
    c = [cos(ang[j]*k) for k = 0:length(d)-1]
    p[j] = 2 * atan((d'*s) / (d'*c)) - ((L+.5)*ang[j])
    p2[j] = (d'*s) / (d'*c) - tan((L+.5)*ang[j]/2)
    diff[j] = abs2.(exp(im*p[j]) - exp(-im*ang[j]/2))
    while (p[j] >= 2*pi)
      p[j] = p[j] - 2*pi
    end
    while (p[j] < 0)
      p[j] = p[j] + 2*pi
    end
  end
  return(p,p2,diff)
end
