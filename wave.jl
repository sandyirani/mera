using PolynomialRoots
using Plots
pyplot()


function applyPoly(p,x)
  n = length(p)
  f = im*zeros(n)
  h = Int8(floor(length(p)/2))
  f[h+1] = 1
  for j = 1:h
    f[h+1+j] = f[h+j]*x
    f[h+1-j] = f[h+2-j]/x
  end
  return(f'*p)
end

function testPoly(p,eps)

  for x = 0:eps:2*pi
    w = exp(im*x)
    val = applyPoly(p,w)
    if (real(val) < 0)
      @show(x, val)
    end
    @show(val)
  end
end

function factorPoly(p)
  t = roots(p)
  h = Int8(floor(length(p)/2))
  newT = im*zeros(h)
  #This assumes that the roots are given in sorted order according to magnitude
  newT = t[1:h]
  q = createPolyFromRootsIm(newT)
  @show(prod(newT))
  #q = q * sqrt(p[1]/prod(newT))
  #q = q * sqrt(p[1])
  if (sum(abs.(imag(q))) > .00001)
      println("Polynomial should have real coefficients")
  end
  return(real(q))
end

function createPolyFromRoots(r)
    p = zeros(length(r)+1)
    a = zeros(length(r)+1)
    b = zeros(length(r)+1)
    p[2] = 1
    p[1] = -r[1]
    for j = 2:length(r)
        a[2:j+1] = p[1:j]
        b[1:j] = -r[j]*p[1:j]
        p = a + b
    end
    return(p)
end

function createPolyFromRootsIm(r)
    p = im*zeros(length(r)+1)
    a = im*zeros(length(r)+1)
    b = im*zeros(length(r)+1)
    p[2] = 1
    p[1] = -r[1]
    for j = 2:length(r)
        a[2:j+1] = p[1:j]
        b[1:j] = -r[j]*p[1:j]
        p = a + b
    end
    return(p)
end

function makeFourierMatrix(l, n, start, final)
   startIdx = Int32(ceil(start*n))
   finalIdx = Int32(ceil(final*n)-1)
   numIdx = (finalIdx - startIdx + 1)
   F = im*zeros(numIdx,l)
   norm = 1/sqrt(2*n)
   norm = 1
   for j = startIdx:finalIdx
   	for k = 1:l
   		F[j-startIdx+1,k] = norm*exp(im*pi*j*(k-1)/n)
   	end
   end
   F
end

function getFourier(v, a, b, n)
  vNorm = v/sqrt(v'*v)
  vNorm = v
  F = makeFourierMatrix(length(v),n, a, b)
  f = F*vNorm
  f = abs.(f)
  f
end

function getTargetPhase(start, final, n)
  startIdx = Int32(ceil(start*n))
  finalIdx = Int32(ceil(final*n)-1)
  numIdx = (finalIdx - startIdx + 1)
  F = [exp(-im*pi*j/(2*n)) for j = startIdx:finalIdx]
  return(F)
end

function getPhaseDiff(v, w, a, b, n)
  vNorm = v/sqrt(v'*v)
  wNorm = w/sqrt(w'*w)
  F = makeFourierMatrix(length(v),n, a, b)
  fv = F*vNorm
  fw = F*wNorm
  t = getTargetPhase(a, b, n)
  return(abs.(fv./fw - t))
end

function checkOffset(v,m)
  sum = 0
  for j = m+1:length(v)
    sum = sum+v[j]*v[j-m]
  end
  return(sum)
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

function levels(g,gp)
  n = length(g)
  h = [(-1)^j*g[j] for j = n:-1:1]
  hp = [(-1)^j*gp[j] for j = n:-1:1]

  hdub = zeros(2*n)
  hpdub = zeros(2*n)
  for j = 1:n
    hdub[2*j-1] = h[j]
    hpdub[2*j-1] = hp[j]
  end
  big = conv(hdub,h)
  bigp = conv(hpdub,hp)
  w = zeros(length(big)+length(bigp))
  for j = 1:length(big)
      w[2*j-1] = big[j]/sqrt(2)
      w[2*j] = bigp[j]/sqrt(2)
  end
  Fbig = getFourier(big,0,1,100)
  Fbigp = getFourier(bigp,0,1,100)
  Fw = getFourier(w,0,1,100)
  output(Fbig,Fbigp,Fw,Fbigp)

end


L = 9
K = 6
tau = 1/2
d = ones(L+1)
for j = 2:L+1
  d[j] = d[j-1]*(L-j+2)*(L-j+2-tau)/(j-1)/(j-1+tau)
end
drev = [d[j] for j = length(d):-1:1]
s1 = [binomial(2*K,j) for j = 0:2*K]
s2 = conv(d,drev)
s = conv(s1,s2)
M = K+L
C = zeros(4*M-1,2*M-1)
for j = 1:4*M-1
  last = min(j,2*M-1)
  first = max(1,j-2*M)
  ffirst = max(1,j-2*M+2)
  flast = last-first+ffirst
  C[j,first:last] = s[flast:-1:ffirst]
end
C = C[2:2:size(C,1),:]
b = zeros(2*M-1)
b[M] = 1

r = C\b

q = factorPoly(r)
#@show(q)
bk = [binomial(K,j) for j = 0:K]
bk2 = [binomial(K,j)^2 for j = 0:K]
f = conv(q,bk)
g = conv(f,drev)
h = conv(f,d)
q = q/sqrt(g'*g)
g = g/sqrt(g'*g)
h = h/sqrt(h'*h)
@show(conv(q,d))

w = zeros(length(h)+length(g))
for j = 1:length(h)
    w[2*j-1] = g[j]/sqrt(2)
    w[2*j] = h[j]/sqrt(2)
end

Fq = getFourier(q,0,1,100)
Fd = getFourier(d,0,1,100)
Fqd = getFourier(conv(q,d),0,1,100)
Fb = getFourier(bk,0,1,100)
Fb2 = getFourier(bk2,0,1,100)
Ff = getFourier(f,0,1,100)
Fw = getFourier(w,0,1,100)
Fg = getFourier(g,0,1,100)
phaseDiff = getPhaseDiff(h, g, 0, 1, 100)
phaseDiff2 = phaseDiff .* Fd
output(Fq,Fd,Fg,Fw)
