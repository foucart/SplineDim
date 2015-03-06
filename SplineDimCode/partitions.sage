## Created by SF on 03/08/2013
##########
# This file contains the definition of classical n-dimensional partitions
##########

##################################
## Alfeld split in n dimensions ##
##################################

# define a simplicial partition in R^n by entering 
# the vertices v_i as a list of points defined by their coordinates
# the simplices t_i as a list of (n+1)-tuples of vertices

def alfeld(n):
   V = [tuple([0 for i in range(n)])]+identity_matrix(n).columns()+[tuple([-1 for i in range(n)])]
   T = [tuple(V[0:i+1]+V[i+2:n+2]) for i in range(n+1)]
   return tuple(T)

