# congruent in python for use in SageMath

import re

def is_congruent(i):
   return i%8 in [5,6,7] or ((i%2 == 1) and odd(i)) or ((i%2 == 0) and even(i))

# print the sides of congruent numbers from lo to hi, inclusive
def compute_sides(lo, hi, init_second_lim = 12, heegner = True):
    todo = [False]*(hi+1)
    for i in range(lo, hi+1):
        if squarefree(i) == 1 and is_congruent(i):
            todo[i] = True

    while sum(todo) > 0:
      print('*** second_lim=', init_second_lim, ' ', sum(todo), 'left to do', flush = True)
      print([i for i in range(len(todo)) if todo[i]])
      for i in range(lo, hi+1):
         if todo[i]:
           x = compute_sides_one(i, init_second_lim, heegner)
           if x == 1:
              todo[i] = False
      init_second_lim += 1

# a simpler version of compute_sides that only tries one value of second_lim
def compute_sides_lim(lo, hi, lim=12, heegner = True):
    for i in range(lo, hi+1):
        if squarefree(i) == 1 and is_congruent(i):
            compute_sides_one(i, lim, heegner)

# return 1 if n squarefree
def squarefree(n):
    if n % 4 == 0:
        return 0
    rt = floor(sqrt(n))
    for i in range(3, rt+1, 2):
       if n % (i*i) == 0:
           return 0
    return 1

# print sides of triangle with area N.  'lim' controls how much searching is done
def compute_sides_one(N, lim=12, heegner = True):

   def print_g(g):
       print_points(g[0].numerator(), g[0].denominator(), g[1].numerator(), g[1].denominator(), N)

   e = EllipticCurve([-N*N, 0])

   if heegner:
     try:
       x,y = pari(e).ellheegner()
     except (RuntimeError) as msg:
       x,y = (0,0)
   else:
     x = 0
   if x == 0:
     try:
         g = e.gens(descent_second_limit = lim)
     except (RuntimeError) as msg:
         #  msg might be of the form "generators found=[[323144640, -6423052104, 274625]])."
         match = re.search(r"\[([0-9-]*),(.*),([0-9- ]*)", str(msg))
         if match is not None:
            preparser(False)
            a = int(match.group(1))
            b = int(match.group(2))
            c = int(match.group(3))
            preparser(True)
            print('      mwrank')
            print_points(a, c, b, c, N)
            return 1
         else:
            print (N, "Cannot find Solution")
            return 0
     else:
        print('      mwrank', len(g), 'generators')
        for i in range(0, len(g)):
            print_g(g[i])
        if len(g) > 1:
            for i in range(len(g)):
              for j in range(i+1, len(g)):
                  print_g(g[i] + g[j])
                  print_g(g[i] - g[j])
        return 1
   else:
     x1 = Rational(x)
     y1 = Rational(y)
     print('      Heegner')
     print_points1(x1.numerator(), x1.denominator(), y1.numerator(), y1.denominator(), N)
     return 1

# for heegner.  Try to find a generator
def print_points1(s,t,u,v, N):
   e = EllipticCurve([-N*N, 0])
   heeg = e(s/t, u/v)
   T = e.torsion_subgroup()

   for i in range(10, 1, -1): # XXX
      for t1 in T:
        tmp = heeg + t1
        x = tmp.division_points(i)
        if len(x) != 0:
           x = x[0]
           print('Heegner divided by', i)
           print_points(x[0].numerator(), x[0].denominator(), x[1].numerator(),
                      x[1].denominator(), N)
           return
   print_points(s, t, u, v, N)

# Get the side lengths alpha = a/d, beta = b/e, compute
# the corresponding values of P,Q, write them as a square
# times a factor of N,  and print them
# use method of https://math.okstate.edu/people/pyan/CongruentNumbersAndEllipticCurves.pdf
def print_points(s,t,u,v, N):  # pts x=s/t, y=u/v on elliptic curve

   x = abs(s/t)
   y = abs(u/v)
   a = abs(x**2 - N**2)   # alpha = (x**2 - N**2)/y
   b = 2*N*x              # beta = 2*N*x/y
   g_a = gcd(a,y);
   g_b = gcd(b, y);
   a = a/g_a;
   d = y/g_a;
   b = b/g_b;
   e = y/g_b;
   if a/d > b/e:
       a,d,b,e = b,e,a,d
   print (N, "   alpha=", a, "/", d)
   print ("      beta=", b, "/", e)

   #  compute A,B,C so that A^2 + B^2 = C^2
   A = a*e
   B = b*d
   g = gcd(A,B)
   A = A/g
   B = B/g
   C = sqrt(A*A + B*B)

   # compute P,Q so that A = 2PQ, B = P^2 - Q^2
   if B % 2 == 1:
      P = sqrt((C+B)/2)
      Q = sqrt((C-B)/2)
      print ("     P=", P, " Q=", Q)
   else:
      P = sqrt((C+A)/2)
      Q = sqrt((C-A)/2)
      print ("     P=", P, " Q=",  Q)

   # factor P = i0*P0^2,  Q = j0*Q0^2
   for i in range(1, N+1):
       if ( P%i == 0 and floor(sqrt(P/i))**2 == P/i):
            P0 = sqrt(P/i)
            i0 = i
            break
   for j in range(1, N+1):
       if (Q%j == 0 and  floor(sqrt(Q/j))**2 == Q/j):
            Q0 = sqrt(Q/j)
            j0 = j
            break
   print ("    P= ", i0, "*", P0, "^2, Q=", j0, "*", Q0, "^2", flush = True)

# test for n congruent when n is odd
def odd(n):
    cnt1 = 0
    cnt2 = 0
    mx = floor(sqrt(n/2))
    my = floor(sqrt(n))
    mz = floor(sqrt(n/32))
    for x in range(-mx, mx+1):
      for y in range(-my, my+1):
        for z in range(-mz, mz+1):
           if 2*x*x + y*y + 32*z*z == n:
               cnt1 = cnt1 + 1
    mz = floor(sqrt(n/8))
    for x in range(-mx, mx+1):
      for y in range(-my, my+1):
        for z in range(-mz, mz+1):
           if 2*x*x + y*y + 8*z*z == n:
               cnt2 = cnt2 + 1
    return(2*cnt1 == cnt2)

# test for n congruent when n is even
def even(n):
    cnt1 = 0
    cnt2 = 0
    n = n/2
    mx = floor(sqrt(n/4))
    my = floor(sqrt(n))
    mz = floor(sqrt(n/32))
    for x in range(-mx, mx+1):
      for y in range(-my, my+1):
        for z in range(-mz, mz+1):
           if 4*x*x + y*y + 32*z*z == n:
               cnt1 = cnt1 + 1
    mz = floor(sqrt(n/8))
    for x in range(-mx, mx+1):
      for y in range(-my, my+1):
        for z in range(-mz, mz+1):
           if 4*x*x + y*y + 8*z*z == n:
               cnt2 = cnt2 + 1
    return(2*cnt1 == cnt2)

