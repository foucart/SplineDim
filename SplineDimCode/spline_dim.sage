## Created by SF on 02/22/2013
## Last modified on 04/18/2013
###########################################################################
## This file contains routines for the computations of spline dimensions ##
## It is mostly based on the code that Patrick Clarke wrote in Feb 2013  ##
###########################################################################


###########################################################################
## A partition Delta can be visualized by calling 
#  type_triangulation(Delta)
## which returns a triangulation object from an iterable of triangles

def type_triangulation(triangles):
   long_points = []
   for t in triangles:
      long_points = long_points + list(t)
   points = set(long_points)
   pc = PointConfiguration(points)
   pc_triangulations = list(pc.triangulations())
   list_triangles = [[list(p) for p in t] for t in triangles]
   for t in list_triangles:
      t.sort()
   list_triangles.sort()
   for it_triang in pc_triangulations:
      list_it_triang = [[list(pc[p]) for p in t] for t in it_triang]
      for t in list_it_triang:
         t.sort()
      list_it_triang.sort()
      if list_it_triang ==  list_triangles:		
         return it_triang
   print "no match"

###########################################################################
## To get the dimension of splines of degree d and smoothness over Delta: 
## spline_dim(Delta,r,d)

def spline_dim(Delta,r,d):
   (OY, s, J) = triangle_splines(Delta, r)
   HS=spline_gen_fun(OY, s, J)
   t = HS.parent().gens()[0]
   taylor_exp=taylor(HS,t,0,d)
   return taylor_exp.coefficients()[d][0]

## To get all the dimensions up to a degree d:
## spline_dims(Delta,r,d)

def spline_dims(Delta,r,d):
   (OY, s, J) = triangle_splines(Delta, r)
   HS=spline_gen_fun(OY, s, J)
   t = HS.parent().gens()[0]
   taylor_exp=taylor(HS,t,0,d)
   return [ c[0] for c in taylor_exp.coefficients() ]

###########################################################################
## In fact, everything is based on the computation of the generating function
## The main routine is called as
## (OY, s, J) = triangle_splines((t0, t1, t2, t3), r)

def triangle_splines(triangle_list, r):
    t_list = list(triangle_list)
    n = len(t_list[0][0])
    s = len(t_list)
    OY = PolynomialRing(QQ,'y',n)
    def J(j,k):
        outJ =  hyperplane_ideal(t_list[j], t_list[k])
        return outJ^(r+1)
    return (OY,s,J)

## this routine requires the following subroutine

def hyperplane_ideal(t1, t2):
   n = len(t1[0])
   OY = PolynomialRing(QQ,'y',n)
   y = OY.gens()
   if t1 == t2:
      return 0*OY
   t1_cap_t2 = [p for p in t1 if p in t2]
   if t1_cap_t2 == []:
      return 1*OY
   int_vects = [vector(p) for p in t1_cap_t2]
   based_vects = [v-int_vects[0] for v in int_vects[1:]]
   M = matrix( based_vects + [y])
   minor_size = min(M.nrows(), M.ncols())
   linear_terms = M.minors(minor_size)
   ideal_generators = [l - l(list(t1_cap_t2[0])) for l in linear_terms]
   h_ideal = ideal(ideal_generators)
   return h_ideal


###########################################################################
## The generating function is obtained as a rational function from data using

def spline_gen_fun(OY, s, J):
   G = G_ideal(OY, s, J)
   Ggb = G.groebner_basis()
   LTgens = [f.lt() for f in Ggb]
   M = M_ideal(OY, s)
   LT = ideal(LTgens)+M
   LThs = LT.hilbert_series()
   t = LThs.parent().gens()[0]
   Mhs = M.hilbert_series()
   shs = (Mhs - LThs)/t
   return shs/(1-t)
## the last line is a modification of Patrick's code, which returned the Hilbert series

# this routines relies on the following routines

def M_ideal(OY, s):
   (OA, e, y) = OAey(OY, s)
   M = ideal(e)^2
   return M

def G_ideal(OY, s, J):
   (OA, e, y) = OAey(OY, s)
   zero_list = [0 for i in range(s)]
   OY_to_OA = OY.hom(y, OA)
   #this corresponds to the projection from AsY to Y
   M = M_ideal(OY, s)
   def jkfilter(j,k):
       return lambda x:((x !=e[j]) and (x !=e[k]))
   def GM(j,k):
      Jjk_in_OA = ideal([ OY_to_OA(h) for h in J(j,k).gens()])
      if s<3:
         return (M + Jjk_in_OA*ideal([e[j], e[k]]) + ideal(e[j]+e[k]))
      return (M + Jjk_in_OA*ideal([e[j], e[k]]) + ideal(e[j]+e[k]) + ideal(filter(jkfilter(j,k), e)))
   def intersect(I,J):
        return I.intersection(J)
   if s < 3:
      return GM(0,1)
   G = reduce(intersect,
                 [GM(j,k) for k in range(1,s) for j in range(k)])
   return G

def OAey(OY, s):
   n = len(OY.gens())  
   OA = PolynomialRing(QQ, 'a', s+n) 
   return (OA, OA.gens()[:s], OA.gens()[s:])

###########################################################################
## The dimensions can be obtained from the generating function by calling

def gf_to_dim(HS,d):
   t = HS.parent().gens()[0]
   taylor_exp=taylor(HS,t,0,d)
   return taylor_exp.coefficients()[d][0]

## The dimensions up to degree d can be obtained by calling

def gf_to_dims(HS,d):
   t = HS.parent().gens()[0]
   taylor_exp=taylor(HS,t,0,d)
   return [ c[0] for c in taylor_exp.coefficients() ]


###########################################################################
## The Hilbert series can be computed from Delta and r using

def spline_gf(Delta,r):
   (OY, s, J) = triangle_splines(Delta, r)
   HS=spline_gen_fun(OY, s, J)
   return HS 

###########################################################################
## The smallest d for which dimension agrees with Hilbert polynomial

def gf_to_dsup(HS):
   HSnum=HS.numerator()
   HSden=HS.denominator()
   ksup=HSnum.degree()
   n=HSden.degree()-1
   return ksup - n

def dsup(Delta,r):
   HS=spline_gf(Delta,r)
   return gf_to_dsup(HS)

###########################################################################
## The largest d for which the spline space is just the polynomial space

def gf_to_dsub(HS):
   HSnum=HS.numerator()
   c=HSnum.coeffs()
   positions=[i for i in range(len(c)) if c[i]!=0]
   return positions[1]-1

def dsub(Delta,r):
   HS=spline_gf(Delta,r)
   return gf_to_dsub(HS)

   
