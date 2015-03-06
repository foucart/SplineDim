## Created by SF on 02/22/2013
# Modified om 03/29/2013
###############################################################################
# This file contains the definitions of some triangulations (=2D-partitions) ##
###############################################################################


########################
## Alfeld split in 2D ##
########################

def alfeld2d():
   (v0, v1, v2, v3) = ( (0,0), (1,0), (0,1), (-1,-1) )
   (t0, t1, t2) = ( (v0, v1, v2), (v0, v2, v3), (v0, v3, v1 ) )
   return (t0, t1, t2)

################################
## Morgan-Scott triangulation ##
################################

## the "regular" version

def morgan_scott_reg():
    (v0, v1, v2, v3, v4, v5) = ( (0,-1), (1,1), (-1,1), (0,4), (-4,-4), (4,-4) )
    (t0, t1, t2, t3, t4, t5, t6) = ( (v0,v1,v2), (v0,v4,v5), (v0,v1,v5), (v1,v3,v5), (v1,v2,v3), (v2,v3,v4), (v0,v4,v2) )
    return (t0, t1, t2, t3, t4, t5, t6)

## the "irregular" version

def morgan_scott_irreg():
    (v0, v1, v2, v3, v4, v5) = ( (0,-1), (1,1), (-1,1), (0,4), (-4,-4), (5,-4) )
    (t0, t1, t2, t3, t4, t5, t6) = ( (v0,v1,v2), (v0,v4,v5), (v0,v1,v5), (v1,v3,v5), (v1,v2,v3), (v2,v3,v4), (v0,v4,v2) )
    return (t0, t1, t2, t3, t4, t5, t6)




