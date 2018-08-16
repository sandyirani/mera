L = 12
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

function pBreakdown(d)
    d2 = copy(d)
    s = [j+taup for j = 0:length(d)-1]
    d2 = d2 .* s
    aSum = zeros(n)
    a = zeros(length(d),n)
    len = Int8(floor(length(d)/2))
    a2 = zeros(len,n)
    for k = 0:n-1
      a[:,k+1] = d2[:]
      a2[:,k+1] = [d2[j]+d2[length(d)-j+1] for j = 1:len]
      aSum[k+1] = sum(d2)
      d2 = d2 .* s
      d2 = d2 .* s
    end
    return(aSum,a,a2)
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
  ph = zeros(length(ang))
  del = zeros(length(ang))
  eps = zeros(length(ang))
  eps2 = zeros(length(ang))
  diff = zeros(length(ang))
  d = makeD()
  for j = 1:length(ang)
    s = [sin(ang[j]*k) for k = 0:length(d)-1]
    c = [cos(ang[j]*k) for k = 0:length(d)-1]
    ph[j] = 2 * atan((d'*s) / (d'*c)) - (L*ang[j])
    diff[j] = abs2.(exp(im*ph[j]) - exp(-im*ang[j]/2))
    del[j] = ph[j] + .5*ang[j]
    eps[j] = (d'*s) / (d'*c) - tan((L-.5)*ang[j]/2)
    eps2[j] = eps[j] * (cos((L-.5)*ang[j]/2))^2
    while (del[j] >= 2*pi)
      del[j] = del[j] - 2*pi
    end
    while (del[j] < 0)
      del[j] = del[j] + 2*pi
    end
  end
  return(ph,del,eps,eps2,diff)
end

function phase2(a,b,m)
  ang = genAngles(a,b,m)
  num = zeros(length(ang))
  den = zeros(length(ang))
  co = zeros(length(ang))
  total = zeros(length(ang))
  d = makeD()
  for j = 1:length(ang)
    s = [sin(ang[j]*(k+taup)) for k = 0:length(d)-1]
    c = [cos(ang[j]*k) for k = 0:length(d)-1]
    num[j] = d' * s
    den[j] = d' * c
    co[j] = cos(ang[j]*taup)
    total[j] = co[j] / (den[j])
  end
  output(num,den,co,total)
  return(num,den,co,total)
end

function output(a,b,c,d)
  @show(length(a))
  @show(length(d))
  f = open("output.txt","w")
  for j = 1:length(a)
    write(f, string(a[j]), ", ", string(b[j]), ", ", string(c[j]), ", ", string(d[j]),"\n")
  end
  close(f)
end

function test(n)
  d = makeD()
  p = makeP(d)
  ang = genAngles(0,1,n)
  @show(length(ang))
  ap = zeros(length(ang))
  for w = 1:length(ang)
    sum = 0
    for k = 1:length(p)
      sum = sum + p[k]*ang[w]^(2*k-1)
    end
    ap[w] = sum
  end
  bp = zeros(length(ang))
  for w = 1:length(ang)
    bp[w] = p[L+1]*ang[w]^(2*L+1)
  end
  cp = zeros(length(ang))
  for w = 1:length(ang)
    cp[w] = p[L+1]*ang[w]^(2*L+1) + p[L+2]*ang[w]^(2*L+3)
  end
  dp = zeros(length(ang))
  (asum,a,a2) = pBreakdown(d)
  plen = length(p)
  num = zeros(plen)
  num[L+1:plen] = a[L+1,L+1:plen]
  for m = 2*L+1:-1:1
      num = num/m
  end
  for k = L+2:n
    num[k:plen] = num[k:plen]/(2*k-1)
    num[k:plen] = num[k:plen]/(2*k-2)
  end
  @show(num)
  @show(p)
  for w = 1:length(ang)
    #d[w] = p[L+1]*ang[w]^(2*L+1) + p[L+2]*ang[w]^(2*L+3) + p[L+3]*ang[w]^(2*L+5)
      sum = 0
      for k = 1:length(p)
        sum = sum + num[k]*ang[w]^(2*k-1)
      end
      dp[w] = sum
  end
  output(ap,bp,cp,dp)
end
