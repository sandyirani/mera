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
  num = 0
  while (num < h)
    j = 1
    while(t[j]==0)
      j = j+1
    end
    newT[num+1] = t[j]
    num = num+1
    pair = t[j]
    t[j] = 0
    j = j+1
    while(abs(inv(conj(pair))-t[j]) < .00001)
      j = j+1
    end
    t[j] = 0
  end
  newT = real(newT)
  q = createPolyFromRoots(newT)
  q = q * sqrt(p[h+1])
  return(q)
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

function makeFourierMatrix(l, n, start, final)

   startIdx = Int32(ceil(start*n))
   finalIdx = Int32(floor(final*n))
   numIdx = (finalIdx - startIdx + 1)

   F = im*zeros(numIdx,l)

   norm = 1/sqrt(2*n)

   for j = startIdx:finalIdx
   	for k = 1:l
   		F[j-startIdx+1,k] = norm*exp(im*pi*j*(k-1)/n)
   	end
   end

   F

end



function getFourier(v, a, b)
  n = 2000
  F = makeFourierMatrix(length(v),n, a, b)
  f = F*v
  f = abs2.(f)
  f

end



L = 4
K = 3
tau = 1/2
d = ones(L+1)
for j = 2:L+1
  d[j] = d[j-1]*(L-j+2)*(L-j+2-tau)/(j-1)/(j-1+tau)
end
s1 = [binomial(2*K,j) for j = 0:2*K]
s2 = conv(d,[d[j] for j = length(d):-1:1])
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
bk = [binomial(K,j) for j = 0:K]
f = conv(q,bk)
h = conv(f,d)
g = conv(f,[d[j] for j = length(d):-1:1])
w = zeros(length(h)+length(g))
for j = 1:length(h)
    w[2*j-1] = g[j]
    w[2*j] = h[j]
end
w = w/sqrt(w'*w)
@show(w'*w)
mag = getFourier(w,0,2);
