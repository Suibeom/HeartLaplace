n=90
function circd(i)
 x, y = div(i, n), rem(i,n)
 if (x-50.5)^2 + (y-50.5)^2 < 1600 && y < 85
  return true
 else
  return false
 end
end
function cardioid(x,y,a)
  return (x^2+y^2)^2 - 4*a*y*(x^2 + y^2) - 4*a^2*x^2
end
function cardd(i)
 x, y = div(i,n),rem(i,n)
 if cardioid(1.3*(x-(n/2)+0.5), 1.3*(y-(n/4)+.5),13) < 0
   return true
 else
   return false
 end
end
function heart(x,y)
  if ((x-10-n/2)^2 + (y+8-n/2)^2 < 169 || ((x+10-n/2)^2 + (y+8-n/2)^2 < 169) ) || (y>(n/2) && (x-n/2)  -(y-n/2) > -21 && -(x-n/2)  -(y-n/2)>-21 )
    return true
  end
  return false
end

function heart90(x,y)
  if ((x-15-n/2)^2 + (y+12-n/2)^2 < 380.25 || ((x+15-n/2)^2 + (y+12-n/2)^2 < 380.25) ) || (y>(n/2) && (x-n/2)  -(y-n/2) > -31.5 && -(x-n/2)  -(y-n/2)>-31.5 )
    return true
  end
  return false
end

function heartxn(x,y,k)
  if ((x-10*k-n/2)^2 + (y+8*k-n/2)^2 < (13*k)^2 || ((x+10*k-n/2)^2 + (y+8*k-n/2)^2 < (13*k)^2) ) || (y>(n/2) && (x-n/2)  -(y-n/2) > -21*k && -(x-n/2)  -(y-n/2)>-21*k )
    return true
  end
  return false
end

function mkhrtxn(k)
  x,y = i->div(i,n),i->rem(i,n)
  return i -> heartxn(x(i),y(i),k)
end


function heart480(x,y)
  if ((x-80-n/2)^2 + (y+64-n/2)^2 < 10816 || ((x+80-n/2)^2 + (y+64-n/2)^2 < 10816) ) || (y>(n/2) && (x-n/2)  -(y-n/2) > -168 && -(x-n/2)  -(y-n/2)>-168 )
    return true
  end
  return false
  end

function hrt480(i)
  x,y = div(i,n),rem(i,n)
  return heart480(x,y)
end

function hrt90(i)
  x,y = div(i,n),rem(i,n)
  return heart90(x,y)
end
function hrt(i)
  x, y = div(i,n),rem(i,n)
  return heart(x,y)
end

function b(i,j)
 if !d(i) || !d(j)
  return 0
 end
 if abs(i-j)==n || abs(i-j) == 1
  return 0.25
end
 if i == j
  return sum(!d(i+k) ? 0.25 : 0 for k in [1,-1,n,-n])
 end
 return 0
end
function bd(i)
  if i == 1 || i == n
    return false
  end
  return true
end

function l1d(i, j,dx)
 b = 0
 if !d(i) || !d(j) # boundary forcing mask
   return b
 end
 if abs(i-j) == 1
  b = 1/dx^2
 end
 if i == j
  b = -2/dx^2
 end
 return b

end

function l1d2(i,j,dx)
  b = 0
  if !d(i) || !d(j) # boundary forcing mask
    return b
  end
  if abs(i-j) == 2
   b = -1/(12*dx^2)
  end
  if abs(i-j) == 1
   b = 1/(3*dx^2)
  end
  if i == j
   b = -5/(2*dx^2)
  end
  return b
end
#This part salvaged from the REPL! Be gentle.
using Arpack, Images, LinearAlgebra, SparseArrays, ProgressMeter
L = spzeros(n^2,n^2)
# This is where the domain can be changed
d = hrt90
#
#This takes a massive amount of time...
#Should parallelize better...
@showprogress for i = 1:n^2
for j in [1,-1,n,-n,0]
if i+j>0 && i+j<=n^2
L[i,i+j] = b(i,i+j)
end
end
end
#Use Arpack's eigs function to do this super quick
e, v= eigs(L, nev=1000, ritzvec=true)
#This will run through the eigenfunctions,
#normalize them, reshape them, and save them!
nameheader = "./imggsssss/hrttt"
for i in 1:1000
V1 = v[:,i]
V2 = (V1 .- min(V1...))/(max(V1...) - min(V1...))
save(nameheader*lpad(string(i),5,string(0))*".jpg", colorview(Gray,reshape(V2,(n,n))))
end
