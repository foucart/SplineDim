## Created by SF on 02/22/2013
# Modified on 03/08/2013
#################################################################################
## This file contains the definitions of some tetrahedrations (=3D-partitions) ##
## They are the ones given in Alfeld's applet
#################################################################################


####################
## Two tetrahedra ##
####################

def two_tetrahedra():
   (v0,v1,v2,v3,v4) = (  (-3,6,0), (3,0,3), (-3,-6,0), (0,0,6), (-2,0,-4) )
   (t0,t1) = ( (v0,v1,v2,v3), (v0,v1,v2,v4) )
   return (t0,t1)

##############
## 3-Orange ##
##############

def three_orange():
   (v0,v1,v2,v3,v4) = ( (-1,-1,0), (0,1,0), (1,-1,0), (0,0,1), (0,0,-1) )
   (t0,t1,t2) = ( (v0,v1,v3,v4), (v1,v2,v3,v4), (v0,v2,v3,v4) ) 
   return (t0,t1,t2)

############################
## Type-I split of a cube ##
############################

def type1():
   (v0, v1, v2, v3, v4, v5, v6, v7) = ( (-1,-1,-1), (1,-1,-1), (1,1,-1), (-1,1,-1), (-1,1,1), (-1,-1,1), (1,-1,1), (1,1,1) )
   (t0, t1, t2, t3, t4, t5) = ( (v0,v1,v2,v7), (v0,v1,v6,v7), (v0,v2,v3,v7), (v0,v3,v4,v7), (v0,v4,v5,v7), (v0,v5,v6,v7) )
   return (t0, t1, t2, t3, t4, t5)

########################################
## Alfeld split (Clough-Tocher in 3D) ##
########################################

def alfeld3d():
    (v0, v1, v2, v3, v4) = ( (0,0,0), (1,0,0), (0,1,0), (0,0,1), (-1,-1,-1) )
    (t0, t1, t2, t3) = ( (v0, v1, v2, v4), (v0, v1, v2, v3), (v0, v2, v3, v4), (v0, v1, v3, v4) )
    return (t0, t1, t2, t3)

#############################################
## Morgan-Scott partition of a tetrahedron ##
#############################################

def Morgan_Scott():
   (v0,v1,v2,v3,v4,v5,v6,v7) = ( (14,42,42), (14,0,42), (-14,21,21), (0,21,0), (2,18,24), (2,24,24), (6,21,27), (4,21,30) )
   (t0,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14) = ( (v4,v5,v6,v7), (v3,v4,v5,v6), (v2,v4,v5,v7), (v1,v4,v6,v7), (v0,v5,v6,v7), (v2,v3,v4,v5), (v1,v3,v4,v6), (v1,v2,v4,v7), (v0,v3,v5,v6), (v0,v2,v5,v7), (v0,v1,v6,v7), (v0,v1,v2,v7), (v0,v1,v3,v6), (v0,v2,v3,v5), (v1,v2,v3,v4) )
   return (t0,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14)


########################
## Regular octahedron ##
########################

def regular_octahedron():
   (v0, v1, v2, v3, v4, v5, v6) = ( (0,0,0), (1,0,0), (0,1,0), (-1,0,0), (0,-1,0), (0,0,1), (0,0,-1) )
   (t0, t1, t2, t3, t4, t5, t6, t7) = ( (v0,v1,v2,v5), (v0,v2,v3,v5), (v0,v3,v4,v5), (v0,v1,v4,v5), (v0,v1,v2,v6), (v0,v2,v3,v6), (v0,v3,v4,v6), (v0,v1,v4,v6) )
   return (t0, t1, t2, t3, t4, t5, t6, t7)


##########################################
## Worsey--Farin split of a tetrahedron ##
##########################################

def Worsey_Farin():
   (v0, v1, v2, v3, v4, v5, v6, v7, v8) = ( (0,0,0), (12,0,0), (0,12,0), (0,0,12), (3,3,3), (4,4,4), (0,4,4), (4,0,4), (4,4,0) )
   (t0, t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11) = ( (v1,v2,v4,v8), (v0,v2,v4,v8), (v0,v1,v4,v8), (v1,v3,v4,v7), (v0,v3,v4,v7), (v0,v1,v4,v7), (v2,v3,v4,v6), (v0,v3,v4,v6), (v0,v2,v4,v6), (v2,v3,v4,v5), (v1,v3,v4,v5), (v1,v2,v4,v5) )
   return (t0, t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11)

###########################
## Double Clough--Tocher ##
###########################

def double_CT():
   (v0,v1,v2,v3,v4,v5,v6,v7,v8) = ( (0,0,0), (16,0,0), (0,16,0), (0,0,16), (4,4,4), (5,5,5), (1,5,5), (5,1,5), (5,5,1) )
   (t0,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,15) = ( (v0,v1,v2,v8), (v1,v2,v4,v8), (v0,v2,v4,v8), (v0,v1,v4,v8), (v0,v1,v3,v7), (v1,v3,v4,v7), (v0,v3,v4,v7), (v0,v1,v4,v7), (v0,v2,v3,v6), (v2,v3,v4,v6), (v0,v3,v4,v6), (v0,v2,v4,v6), (v1,v2,v3,v5), (v2,v3,v4,v5), (v1,v3,v4,v5), (v1,v2,v4,v5) )
   return (t0,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,15)

#############################
## Type-IV split of a cube ##
#############################
# for the visualization, it is stuck at the construction of the triangulation

def type4():
   (v0,v1,v2,v3,v4,v5,v6,v7,v8,v9,v10,v11,v12,v13,v14) = ( (0,0,0), (2,0,0), (0,2,0), (0,0,2), (1,1,1), (0,2,2), (2,2,2), (2,0,2), (2,2,0), (1,1,2), (0,1,1), (1,0,1), (2,1,1), (1,2,1), (1,1,0) )
   (t0,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23) = ( (v3,v4,v5,v9), (v4,v5,v6,v9), (v4,v6,v7,v9), (v3,v4,v7,v9), (v3,v4,v5,v10), (v2,v4,v5,v10), (v0,v2,v4,v10), (v0,v3,v4,v10), (v0,v1,v4,v11), (v1,v4,v7,v11), (v3,v4,v7,v11), (v0,v3,v4,v11), (v4,v6,v7,v12), (v1,v4,v7,v12), (v1,v4,v8,v12), (v4,v6,v8,v12), (v4,v5,v6,v13), (v4,v6,v8,v13), (v2,v4,v8,v13), (v2,v4,v5,v13), (v2,v4,v8,v14), (v1,v4,v8,v14), (v0,v1,v4,v14), (v0,v2,v4,v14) )
   return (t0,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23)


###########################
## Generic Morgan--Scott ##
###########################

def generic_Morgan_Scott():
   (v0,v1,v2,v3,v4,v5,v6,v7) = ( (140,420,420), (140,0,420), (-140,210,210), (0,21,0), (20,127,242), (20,185,238), (60,157,274), (39,181,297) )
   (t0,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14) = ( (v4,v5,v6,v7), (v3,v4,v5,v6), (v2,v4,v5,v7), (v1,v4,v6,v7), (v0,v5,v6,v7), (v2,v3,v4,v5), (v1,v3,v4,v6), (v1,v2,v4,v7), (v0,v3,v5,v6), (v0,v2,v5,v7), (v0,v1,v6,v7), (v0,v1,v2,v7), (v0,v1,v3,v6), (v0,v2,v3,v5), (v1,v2,v3,v4) )
   return (t0,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14)

########################
## Generic octahedron ##
########################

def generic_octahedron():
   (v0, v1, v2, v3, v4, v5, v6) = ( (1,2,3), (10,0,1), (3,11,2), (-12,0,0), (-1,-13,2), (1,-2,14), (2,-3,-15) )
   (t0, t1, t2, t3, t4, t5, t6, t7) = ( (v0,v1,v2,v5), (v0,v2,v3,v5), (v0,v3,v4,v5), (v0,v1,v4,v5), (v0,v1,v2,v6), (v0,v2,v3,v6), (v0,v3,v4,v6), (v0,v1,v4,v6) )
   return (t0, t1, t2, t3, t4, t5, t6, t7)

###########################
## Generic Worsey--Farin ##
###########################

def generic_Worsey_Farin():
   (v0,v1,v2,v3,v4,v5,v6,v7,v8) = ( (0,0,0), (24,0,0), (0,24,0), (0,0,24), (6,7,8), (8,8,8), (0,8,8), (8,0,8), (8,8,0) )
   (t0,t1,t2,t3,t4,t5,t6,t7,t8,t9,10,t11) = ( (v1,v2,v4,v8), (v0,v2,v4,v8), (v0,v1,v4,v8), (v1,v3,v4,v7), (v0,v3,v4,v7), (v0,v1,v4,v7), (v2,v3,v4,v6), (v0,v3,v4,v6), (v0,v2,v4,v6), (v2,v3,v4,v5), (v1,v3,v4,v5), (v1,v2,v4,v5) )
   return (t0,t1,t2,t3,t4,t5,t6,t7,t8,t9,10,t11)

###################################
## Generic double Clough--Tocher ##
###################################

def generic_double_CT():
   (v0,v1,v2,v3,v4,v5,v6,v7,v8) = ( (0,0,0), (32,0,0), (0,32,0), (0,0,32), (8,9,10), (11,10,9), (1,10,11), (12,3,10), (9,9,1) )
   (t0,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15) = ( (v0,v1,v2,v8), (v1,v2,v4,v8), (v0,v2,v4,v8), (v0,v1,v4,v8), (v0,v1,v3,v7), (v1,v3,v4,v7), (v0,v3,v4,v7), (v0,v1,v4,v7), (v0,v2,v3,v6), (v2,v3,v4,v6), (v0,v3,v4,v6), (v0,v2,v4,v6), (v1,v2,v3,v5), (v2,v3,v4,v5), (v1,v3,v4,v5), (v1,v2,v4,v5) )
   return (t0,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15)

#####################################
## Generic type-IV split of a cube ##
#####################################
# for the visualization, it is stuck at the construction of the triangulation

def generic_type4():
   (v0,v1,v2,v3,v4,v5,v6,v7,v8,v9,v10,v11,v12,v13,v14) = ( (0,0,0), (30,1,-2), (-1,30,3), (2,-3,30), (15,15,17), (1,26,31), (32,30,34), (32,-3,24), (29,32,5), (15,12,31), (-2,12,16), (18,0,14), (29,14,17), (13,28,17), (17,17,-1) )
   (t0,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23) = ( (v3,v4,v5,v9), (v4,v5,v6,v9), (v4,v6,v7,v9), (v3,v4,v7,v9), (v3,v4,v5,v10), (v2,v4,v5,v10), (v0,v2,v4,v10), (v0,v3,v4,v10), (v0,v1,v4,v11), (v1,v4,v7,v11), (v3,v4,v7,v11), (v0,v3,v4,v11), (v4,v6,v7,v12), (v1,v4,v7,v12), (v1,v4,v8,v12), (v4,v6,v8,v12), (v4,v5,v6,v13), (v4,v6,v8,v13), (v2,v4,v8,v13), (v2,v4,v5,v13), (v2,v4,v8,v14), (v1,v4,v8,v14), (v0,v1,v4,v14), (v0,v2,v4,v14) )
   return (t0,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23)


####################
## Generic 8-cell ##
####################

def generic_8cell():
   (v0, v1, v2, v3, v4, v5, v6) = ( (22,24,31), (7,8,9), (72,0,24), (0,72,0), (0,0,72), (26,24,33), (-18,28,23) )
   (t0, t1, t2, t3, t4, t5, t6, t7) = ( (v0,v1,v2,v3), (v0,v1,v2,v4), (v0,v1,v3,v6), (v0,v1,v4,v6), (v0,v3,v4,v6),  (v0,v2,v4,v5), (v0,v2,v3,v5), (v0,v3,v4,v5) )
   return (t0, t1, t2, t3, t4, t5, t6, t7)


##############
## Inverted ##
##############

def inverted():
   (v0,v1,v2,v3,v4,v5,v6,v7) = ( (0,0,0), (3,0,0), (0,3,0), (0,0,3), (1,1,1), (0,1,1), (1,0,1), (1,1,0) )
   (t0,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10) = ( (v4,v5,v6,v7), (v0,v5,v6,v7), (v1,v4,v6,v7), (v2,v4,v5,v7), (v3,v4,v5,v6), (v0,v1,v6,v7), (v0,v2,v5,v7), (v1,v2,v4,v7), (v0,v3,v5,v6), (v1,v3,v4,v6), (v2,v3,v4,v5) )
   return (t0,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10)


#######################
## Symmetric MS cone ##
#######################

def symmetric_MS_cone():
   (v0,v1,v2,v3,v4,v5,v6) = ( (5,5,-5), (0,0,0), (0,15,0), (15,0,0), (6,6,0), (6,3,0), (3,6,0) )
   (t0,t1,t2,t3,t4,t5,t6) = ( (v0,v4,v5,v6), (v0,v1,v5,v6), (v0,v2,v4,v6), (v0,v3,v4,v5), (v0,v1,v3,v5), (v0,v1,v2,v6), (v0,v2,v3,v4) )
   return (t0,t1,t2,t3,t4,t5,t6)

#####################
## Generic MS cone ##
#####################

def generic_MS_cone():
   (v0,v1,v2,v3,v4,v5,v6) = ( (50,50,-50), (0,0,0), (0,150,0), (150,0,0), (60,63,0), (57,33,0), (30,60,0) )
   (t0,t1,t2,t3,t4,t5,t6) = ( (v0,v4,v5,v6), (v0,v1,v5,v6), (v0,v2,v4,v6), (v0,v3,v4,v5), (v0,v1,v3,v5), (v0,v1,v2,v6), (v0,v2,v3,v4) )
   return (t0,t1,t2,t3,t4,t5,t6)






### Other triangulations given on Peter Alfeld's web page ###


##########################################################
## A triangular prism partitioned into three tetrahedra ##
##########################################################

def prism():
   (v0, v1, v2, v3, v4, v5) = ( (0, 0, 0), (0, 1, 0), (0, 0, 1), (1, 0, 0), (1, 1, 0), (1, 0, 1) )
   (t0, t1, t2) = ( (v0, v1, v2, v3), (v1, v2, v3, v4), (v2, v3, v4, v5) )
   return (t0, t1, t2)
	







   


