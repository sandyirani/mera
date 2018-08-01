function taylorSin(w,n)
  s = zeros(n)
  s[1] = w
  last = w
  for j = 1:n-1
    last = -1 * last * w^2 / (2*j) / (2*j+1)
    s[j+1] = s[j] + last
  end
  return(s)
end

function taylorCos(w,n)
  c = zeros(n)
  c[1] = 1
  last = 1
  for j = 1:n-1
    last = -1 * last * w^2 / (2*j) / (2*j-1)
    c[j+1] = c[j] + last
  end
  return(c)
end
