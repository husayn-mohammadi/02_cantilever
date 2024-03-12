import sys, time
# import pandas as pd
import numpy as np
import functions.FuncAnalysis  as fa
exec(open("Input/unitsSI.py").read()) 
# exec(open("Input/inputData.py").read()) 
recordToLogDesign = True
fa.replace_line('MAIN.py', 27, "recordToLog     = False                      # True, False")
if recordToLogDesign == True:
    sys.stdout = open("logDesign.txt", 'w') 

#||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
#                   Step 1: Input Data
#||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

# Seismic Coefficients
# Cd      = 5.5
# Ie      = 1.
# R       = 8
# Omega0  = 2.5
# Rho     = 1.

# Composite wall resistance factor
# Fi_v = Fi_b = Fi_c = Fi_t = 0.9

# Building Geometry
# Lf      = 200   *ft
# Wf      = 120   *ft
# h_typ   = 14    *ft
# h_1     = 17    *ft
# n_story = 10

# Floor Loads
# DL      = 0.12  *ksf
# LL      = 0

# Material
# Es      = 29000 *ksi
# Gs      = 11200 *ksi
# Fy      = 50    *ksi
# Fu      = 65    *ksi
# Ry      = 1.1

# fpc     = 6     *ksi
# Ec      = 4500  *ksi
# Gc      = 1800  *ksi
# Rc      = 1.3

# Wall and Coupling Beam Sections
# Lw      = 132   *inch
# tw      = 9/16  *inch
# t_sc    = 24    *inch
# btie    = 12    *inch
# Stie    = 12    *inch
# dtie    = 1     *inch

# L_CB    = 96    *inch
# b_CB    = 24    *inch
# h_CB    = 24    *inch
# t_pfCB  = 0.5   *inch
# t_pwCB  = 3/8   *inch


# Coupling Beam Properties
# As_CB   = 2 *t_pwCB *h_CB + 2 *t_pfCB *(b_CB -2 *t_pwCB)
# Ac_CB   = (b_CB -2 *t_pwCB) *(h_CB -2 *t_pfCB)
# Asw_CB  = 2 *h_CB *t_pwCB
# Is_CB   = 2 *(1/12 *t_pwCB *h_CB **3 + 1/12 *(b_CB -2 *t_pwCB) *t_pfCB **3 +
#             t_pfCB *(b_CB -2 *t_pwCB) *((h_CB -t_pfCB) /2) **2)
# Ic_CB   = 1/12 *(b_CB -2 *t_pwCB) *(h_CB -2 *t_pfCB) **3

# EAuncrCB= 0.8 *(Es *As_CB +Ec *Ac_CB)
# C3      = min(0.9, 0.45 +3 *(As_CB /b_CB/ h_CB))
# EIeff_CB= 0.64 *(Es *Is_CB +C3 *Ec *Ic_CB)
# GAv_CB  = Gs *Asw_CB + Gc *Ac_CB

# Planar SpeedCore Wall Properties
# As      = tw *(2 *(Lw -2 *tw) + 2 *t_sc)
# Ac      = Lw *t_sc -As
# Asw     = 2 *Lw *tw
# Is      = (2 *(1/12 *tw *(Lw -2 *tw) **3) + 
#            2 *(1/12 *t_sc *tw **3 + t_sc*tw *((Lw -tw) /2) **2))
# Ic      = 1/12 * (t_sc -2 *tw) * (Lw -2 *tw) **3

# EAeff   = Es *As  +0.45 *Ec *Ac
# EIeff   = Es *Is  +0.35 *Ec *Ic
# GAveff  = Gs *Asw +      Gc *Ac

#______________________________________________________________________________
#                   Check with Thresholds
"""'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''"""

# 1.    Check Minimum and Maximum Area of Steel
"""'''''''''''''''''''''''''''''''''''''''''"""
# rho_min  = 0.01
# rho_max  = 0.1

# 1.1.  Coupling Beams
# As_CBmin = rho_min *b_CB *h_CB
# As_CBmax = rho_max *b_CB *h_CB
# if As_CB < As_CBmin: 
#     print("Area of Steel is less than minimum for Coupling Beams!!!")
# elif As_CB > As_CBmax:
#     print("Area of Steel is greater than maximum for Coupling Beams!!!")

# 1.2.  SpeedCore Walls
# Asmin = rho_min *t_sc *Lw
# Asmax = rho_max *t_sc *Lw
# if As < Asmin: 
#     print("Area of Steel is less than minimum for Composite Walls!!!")
# elif As > Asmax:
#     print("Area of Steel is greater than maximum for Composite Walls!!!")


# 2. Check Plate Slenderness Ratios
"""'''''''''''''''''''''''''''''"""
# 1.1.  Coupling Beams
# lambdaP_fCB     = 2.37 *(Es /Ry /Fy) **0.5
# lambdaP_wCB     = 2.66 *(Es /Ry /Fy) **0.5

# b_cCB           = b_CB -2 *t_pwCB
# lambda_fCB      = b_cCB /t_pfCB

# h_cCB           = h_CB -2 *t_pfCB
# lambda_wCB      = h_cCB /t_pwCB

# if lambda_fCB/lambdaP_fCB > 1.0: 
#     print(f"Coupling Beam Flange Plate Slenderness OVERRATED!!!\n===>>>Lam_fCB/LamP_fCB = {lambda_fCB/lambdaP_fCB:.2f}")
# if lambda_wCB/lambdaP_wCB > 1.0: 
#     print(   f"Coupling Beam Web Plate Slenderness OVERRATED!!!\n===>>>Lam_wCB/LamP_wCB = {lambda_wCB/lambdaP_wCB:.2f}")

# 1.2.  SpeedCore Walls
# lambda_btie = btie/tw
# lambdaP_btie= 1.05 *(Es /Ry /Fy) **0.5
# lambdaP_btieTop = 1.2 *(Es /Fy) **0.5

# if lambda_btie/lambdaP_btie > 1.0:    
#     print(f"Wall Faceplate Slenderness in flexural yielding zone not CHECKED!!!\n===>>>Lam_fW/LamP_fW = {lambda_btie/lambdaP_btie:.2f}") 
# if lambda_btie/lambdaP_btieTop > 1.0: 
#     print(f"Wall Faceplate Slenderness above flexural yielding zone not CHECKED!!!\n===>>>Lam_fW/LamP_fW = {lambda_btie/lambdaP_btieTop:.2f}") 

# Alpha       = 1.7 *(t_sc /tw -2) *(tw /dtie) **4
# lambdaP_tbs = (Es /(2 *Alpha +1)) **0.5
# lambda_tbs  = Stie /tw

# if lambda_tbs/lambdaP_tbs > 1.0: 
#     print(f"Wall Tie bar spacing not CHECKED!!!\n===>>>Lam_tbs/LamP_tbs = {lambda_tbs/lambdaP_tbs:.2f}") 


# 3. Check Shear Criticality of Coupling Beams
"""''''''''''''''''''''''''''''''''''''''''"""
# t_cCB   = b_CB -2 *t_pwCB

# C_CBexp = ((2 *t_pwCB *Ry *Fy*h_CB  +0.85 *Rc *fpc *t_cCB *t_pfCB)/
#            (4 *t_pwCB *Ry *Fy       +0.85 *Rc *fpc *t_cCB))
# print(f"C_CBexp = {C_CBexp *1000:.1f} mm")
# C_1exp  = (b_CB -2 *t_pwCB) *t_pfCB *Ry *Fy
# C_2exp  = 2 *t_pwCB *C_CBexp *Ry *Fy
# C_3exp  = 0.85 *Rc *fpc *t_cCB *(C_CBexp -t_pfCB)
# T_1exp  = (b_CB -2 *t_pwCB) *t_pfCB *Ry *Fy
# T_2exp  = 2 *t_pwCB *(h_CB -C_CBexp) *Ry *Fy
# M_exp = (C_1exp *(C_CBexp -t_pfCB /2) +
#             C_2exp *(C_CBexp /2) +
#             C_3exp *((C_CBexp -t_pfCB) /2) +
#             T_1exp *(h_CB -C_CBexp -t_pfCB /2) +
#             T_2exp *((h_CB -C_CBexp) /2))
# V_exp    = 0.6 *Ry *Fy *Asw_CB +0.06 *Ac_CB *(Rc *fpc /MPa) **0.5 #Take care NOT to put fpc in other units except MPa

# print(f"M_exp = {M_exp /1000:.1f} kN.m")
# print(f"V_exp = {V_exp /1000:.1f} kN")

# # Check Flexure-Criticality Condition
# if V_exp *L_CB /M_exp >= 2.4: 
#     print("Flexure-Criticality Confirmed!")

numSign = 65
#||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
#                   Step 2: Analysis for Design
#||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
fa.replace_line('MAIN.py', 31, "linearity       = True")
fa.replace_line('MAIN.py', 32, "typeBuild       = 'coupledWalls'            # 'CantileverColumn', 'coupledWalls', 'buildBeam', 'ShearCritBeam'")
fa.replace_line('MAIN.py', 78, "plot_MomCurv    = False")
t_EAna_i    = time.time()
print(f"{'='*numSign}\nElastic Analysis Started.\n{'='*numSign}\n")
exec(open("MAIN.py").read()) 
t_EAna_f    = time.time()
dur_EAna    = (t_EAna_f - t_EAna_i)/60
print(f"{'='*numSign}\nElastic Analysis Finished in {dur_EAna:.2f} mins.\n{'='*numSign}\n")

# Effective distance between wall centroids
L_eff   = L_CB +Lw

#---------------------------------------------
# Check Maximum Inter-Story Drift here using amplified displacement. ASCE7 12.12.1
ID_allowable = 0.02 # for example here for Risk Category II. 
IDe     = driftMaximum
ID = IDe *Cd
if ID > ID_allowable:
    print(f"ID = {ID*100:.5f}% which is greater than {ID_allowable*100}%\nProgram exits here!"); sys.exit()
else:
    print(f"ID = {ID*100:.5f}%")


# Loads 
# SDloads = pd.read_excel(io="Input/Excel/Example.xlsx", sheet_name="Loads", 
#                     usecols      = "A:D",
#                     skiprows     = 5,
#                     nrows        = 11)

# Base Shear
# V_base  = 879 *kip
Cvx     = fa.verDistFact(We, T1, h_1, h_typ, n_story)

OTM     = 0
for i in range(1, n_story+1):
    def h(n):
        if n == 1:
            return h_1
        else:
            return (h_1 + (n-1) *h_typ)
    OTM += h(i) * Cvx[i] *V_base
    # print(f"OTM = {OTM/kip/inch:.1f} kip-in")
    # pass

# OTM     = 890000 *kip*inch
print(f"OTM = {OTM /1000:.1f} kN.m = {OTM/kip/inch:.1f} kip-in")
print(f"V_base = {V_base /1000:.1f} kN = {V_base /kip:.1f} kip")
# Coupling Beams Required Strength
# Vr_CB   = 324 *kip # This is an average


#||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
#                   Step 3: Design of Coupling Beams
#||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
# 3-1   Calculate Required Strength
"""'''''''''''''''''''''''''''''"""
# a)    Shear Strength
Vu_CB   = Vr_CB
print(f"Vu_CB = {Vu_CB /1000:.1f} kN")
# b)    Flexural Strength
Mu_CB   = 1/2 *Vu_CB *L_CB
print(f"Mu_CB = {Mu_CB /1000:.1f} kN.m")

# 3-2   Calculate Available Strength
"""''''''''''''''''''''''''''''''"""
# a)    Shear Strength
Kc          = 1 # AISC 360, Section I4.2 (Usually equals to 1)
Vn_CB       = 0.6 *Fy *Asw_CB +0.06 *Kc *Ac_CB *(fpc /MPa) **0.5
Vn_CB_Fi_v  = Fi_v *Vn_CB
print(f"Vn_CB*Fi_v = {Vn_CB_Fi_v /1000:.1f} kN")

# b)    Flexural Strength
C_CB        = ((2 *t_pwCB *h_CB *Fy +0.85 *fpc *t_cCB *t_pfCB)/
               (4 *t_pwCB *Fy       +0.85 *fpc *t_cCB))
print(f"C_CB = {C_CB *1000:.1f} mm")
C1          = (b_CB -2 *t_pwCB) *t_pfCB *Fy
C2          = 2 *t_pwCB *C_CB *Fy
C3          = 0.85 *fpc *t_cCB *(C_CB -t_pfCB)
T1          = (b_CB -2 *t_pwCB) *t_pfCB *Fy
T2          = 2 *t_pwCB *(h_CB -C_CB) *Fy

Mpn_CB      = (C1 *(C_CB -t_pfCB /2) + 
               C2 *(C_CB /2) + 
               C3 *((C_CB -t_pfCB) /2) + 
               T1 *(h_CB -C_CB -t_pfCB /2) + 
               T2 *((h_CB -C_CB) /2))
Mn_CB       = Mpn_CB
Mn_CB_Fi_b  = Fi_b *Mn_CB
print(f"Mn_CB*Fi_b = {Mn_CB_Fi_b /1000:.1f} kN.m")

# 3-3   Check Strength Ratios
"""'''''''''''''''''''''''"""
# a)    Shear Strength
R__V_CB     = Vu_CB /Vn_CB_Fi_v
if R__V_CB >= 1.0:
    print("The Available Shear Strength of Coupling Beam is NOT SUFFICIENT!!!")
elif 0.95 < R__V_CB < 1.0:
    print("The Available Shear Strength of Coupling Beam is OK")
else:
    print(f"The Available Shear Strength of Coupling Beam is OK but NOT OPTIMUM: \n===>>>Ratio = {R__V_CB:.2f}")

# b)    Flexural Strength
R__M_CB     = Mu_CB /Mn_CB_Fi_b
if R__M_CB > 1.0:
    print("The Available Flexural Strength of Coupling Beam is NOT SUFFICIENT!!!")
elif 0.95 < R__M_CB <= 1.0:
    print("The Available Flexural Strength of Coupling Beam is OK")
else:
    print(f"The Available Flexural Strength of Coupling Beam is OK but NOT OPTIMUM: \n===>>>Ratio = {R__M_CB:.2f}")


#||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
#                   Step 4: Design of Composite Walls
#||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

# 4-1   Calculate Required Strength
"""'''''''''''''''''''''''''''''"""
# a)    Required Axial Strength
Vn_Mp_exp   = 2 *(1.2 *M_exp) /L_CB
print(f"Vn_Mp_exp = {Vn_Mp_exp /1000:.1f} kN")
Pu          = n_story *Vn_Mp_exp
print(f"Pu = {Pu /1000:.1f} kN")

# b)    Required Shear Strength
V_amp   = 4 *V_base
print(f"V_amp = {V_amp /1000:.1f} kN")
Vu      = V_amp /2
print(f"Vu = {Vu /1000:.1f} kN")

# c)    Required Flexural Strength All Walls
gamma1      = (n_story *(1.2 *M_exp)) /(n_story *Mu_CB)
Mu_Both     = gamma1 *OTM -Pu *L_eff
print(f"Mu_Both = {Mu_Both /1000:.1f} kN.m")

# EIeff_Ten   = 7.7e9   *kip*inch **2
# EIeff_Com   = 1.81e10 *kip*inch **2
fa.replace_line('MAIN.py', 31, "linearity       = False")
fa.replace_line('MAIN.py', 32, "typeBuild       = 'CantileverColumn'            # 'CantileverColumn', 'coupledWalls', 'buildBeam', 'ShearCritBeam'")
fa.replace_line('MAIN.py', 78, "plot_MomCurv    = True")
t_IEAna_i   = time.time()
print(f"{'='*numSign}\nInelastic Analysis Started.\n{'='*numSign}\n")
exec(open("MAIN.py").read()) 
t_IEAna_f   = time.time()
dur_IEAna   = (t_IEAna_f - t_IEAna_i)/60
print(f"{'='*numSign}\nInelastic Analysis Finished in {dur_IEAna:.2f} mins.\n{'='*numSign}\n")
EIeff_Ten   = EIeff_walls[0]
EIeff_Com   = EIeff_walls[1]

# The portion of the overturning moment resisted by the individual wall is
Mu_T        = (EIeff_Ten /(EIeff_Ten +EIeff_Com)) *Mu_Both
print(f"Mu_T = {Mu_T /1000:.1f} kN.m")
Mu_C        = (EIeff_Com /(EIeff_Ten +EIeff_Com)) *Mu_Both
print(f"Mu_C = {Mu_C /1000:.1f} kN.m")



# 4-2   Calculate Available Strength 
"""''''''''''''''''''''''''''''''"""
# a)    Axial Strength
# a.1)  Tensile Strength
Pn_T        = Fy *As
Pn_T_Fi_t   = Fi_t *Pn_T
print(f"Pn_T*Fi_t = {Pn_T_Fi_t /1000:.1f} kN")

# a.2)    Compressive Strength
Pno         = Fy *As +0.85 *fpc *Ac
Is_minor    = (2 *(1/12 *(Lw -2 *tw) *tw **3 + (Lw -2 *tw) *tw *((t_sc -tw) /2) **2) + 2 *(1/12 *tw *t_sc **3))
Ic_minor    = 1/12 *(Lw -2 *tw) *(t_sc -2 *tw) **3
EIeff_minor = Es *Is_minor +0.35 *Ec *Ic_minor
L_cr        = max(h_1, h_typ)
Pe          = np.pi **2 *EIeff_minor /L_cr **2

if Pno/Pe < 2.25:
    Pn_C    = Pno *(0.658 **(Pno /Pe))
else:
    Pn_C    = 0.877 *Pe

Pn_C_Fi_c   = Fi_c *Pn_C
print(f"Pn_C*Fi_c = {Pn_C_Fi_c /1000:.1f} kN")


# b)    Shear Strength
Ks          = Gs *Asw
Ksc         = ((0.7 *(Es *Asw) *(Ec *Ac))/
               (  4 *(Es *Asw) +(Ec *Ac)))

Vn          = ((Ks +Ksc)
               /((3 *Ks **2 +Ksc  **2)  **0.5)
               *(Fy *Asw))

Vn_Fi_v     = Fi_v *Vn
print(f"Vn*Fi_v = {Vn_Fi_v /1000:.1f} kN")


# c)    Flexural Strength
# c.1)  Flexural Strength of Tension SpeedCore Wall
C_T         = (( 0.85 *fpc *(t_sc -2 *tw) *tw +2 *tw *Fy *Lw -Pu)/
                (0.85 *fpc *(t_sc -2 *tw)      +4 *tw *Fy))
print(f"C_T = {C_T *1000:.1f} mm")
C1_T        = (t_sc -2 *tw) *tw *Fy
C2_T        = 2 *tw *C_T *Fy
C3_T        = 0.85 *fpc *(t_sc -2 *tw) *(C_T -tw)
T1_T        = (t_sc -2 *tw) *tw *Fy
T2_T        = 2 *tw *(Lw -C_T) *Fy

Mp_T        = (C1_T *(C_T -tw /2) +
               C2_T *(C_T /2) +
               C3_T *((C_T -tw) /2) +
               T1_T *(Lw -C_T -tw /2) +
               T2_T *((Lw -C_T) /2) +
               Pu *(C_T -Lw /2))

Mn_T        = Mp_T
Mn_T_Fi_b   = Fi_b *Mn_T
print(f"Mn_T*Fi_b = {Mn_T_Fi_b /1000:.1f} kN.m")
    
# c.2)  Flexural Strength of Compression SpeedCore Wall
C_C         = ((0.85 *fpc *(t_sc-2 *tw) *tw +2 *tw *Fy *Lw +Pu)/
               (0.85 *fpc *(t_sc-2 *tw)      +4 *tw *Fy))
print(f"C_C = {C_C *1000:.1f} mm")

C1_C        = (t_sc -2 *tw) *tw *Fy
C2_C        = 2 *tw *C_C *Fy
C3_C        = 0.85 *fpc *(t_sc -2 *tw) *(C_C -tw)
T1_C        = (t_sc -2 *tw) *tw *Fy
T2_C        = 2 *tw *(Lw -C_C) *Fy

Mp_C        = (C1_C *(C_C -tw /2) +
               C2_C *(C_C /2) +
               C3_C *((C_C -tw) /2) +
               T1_C *(Lw -C_C -tw /2) +
               T2_C *((Lw -C_C) /2) +
               Pu* (Lw /2 -C_C)) 

Mn_C        = Mp_C
Mn_C_Fi_b   = Fi_b *Mn_C
print(f"Mn_C*Fi_b = {Mn_C_Fi_b /1000:.1f} kN.m")
print("\n")


# 4-3   Check Strength Ratios
"""'''''''''''''''''''''''"""
# a)    Axial Strength
# a.1)  Tensile Strength
R__P_Twall  = Pu/Pn_T_Fi_t
print(f"\nR__P_Twall = {R__P_Twall*100:.1f}%")
if R__P_Twall > 1.0:
    print("The Available Tensile Strength of Tension Wall is NOT OK!!!")
elif 0.95 < R__P_Twall <= 1.0:
    print("The Available Tensile Strength of Tension Wall is OK")
else:
    print("The Available Tensile Strength of Tension Wall is OK, but NOT OPTIMUM!")# \n===>>>Ratio = {R__P_Twall*100:.1f}%")


# a.2)  Compressive Strength
R__P_Cwall  = Pu/Pn_C_Fi_c
print(f"\nR__P_Cwall = {R__P_Cwall*100:.1f}%")
if R__P_Cwall > 1.0:
    print("The Available Compressive Strength of Compression Wall is NOT OK!!!")
elif 0.95 < R__P_Cwall <= 1.0:
    print("The Available Compressive Strength of Compression Wall is OK")
else:
    print("The Available Compressive Strength of Compression Wall is OK, but NOT OPTIMUM!")# \n===>>>Ratio = {R__P_Cwall*100:.1f}%")

# b)    Shear Strength
R__V_wall   = Vu /Vn_Fi_v
print(f"\nR__V_wall = {R__V_wall*100:.1f}%")
if R__V_wall > 1.0:
    print("The Available Shear Strength of Walls is NOT OK!!!")
elif 0.95 < R__V_wall <= 1.0:
    print("The Available Shear Strength of Walls is OK")
else:
    print("The Available Shear Strength of Walls is OK, but NOT OPTIMUM!")# \n===>>>Ratio = {R__V_wall*100:.1f}%")

# c)    Flextural Strength
# c.1)  Flexural Strength of Tension SpeedCore Wall
R_M_Twall   = Mu_T/Mn_T_Fi_b
print(f"\nR_M_Twall = {R_M_Twall*100:.1f}%")
if R_M_Twall > 1.0:
    print("The Available Flexural Strength of Tension Wall is NOT OK!!!")
elif 0.95 < R_M_Twall <= 1.0:
    print("The Available Flexural Strength of Tension Wall is OK")
else:
    print("The Available Flexural Strength of Tension Wall is OK, but NOT OPTIMUM!")# \n===>>>Ratio = {R_M_Twall*100:.1f}%")

# c.2)  Flexural Strength of Compression SpeedCore Wall
R__M_Cwall   = Mu_C/Mn_C_Fi_b
print(f"\nR__M_Cwall = {R__M_Cwall*100:.1f}%")
if R__M_Cwall > 1.0:
    print("The Available Flexural Strength of Compression Wall is NOT OK!!!")
elif 0.95 < R__M_Cwall <= 1.0:
    print("The Available Flexural Strength of Compression Wall is OK")
else:
    print("The Available Flexural Strength of Compression Wall is OK, but NOT OPTIMUM!")# \n===>>>Ratio = {R__M_Cwall*100:.1f}%")




#||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
#                   Coupling Ratio
#||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

M_couplingBeams     = L_eff * Pu
M_all               = Mu_Both + M_couplingBeams
R__coupling         = M_couplingBeams /M_all

print(f"\n\nCoupling Ratio = {R__coupling*100:.1f}%")


if recordToLogDesign == True:
    sys.stdout.close()
    sys.stdout = sys.__stdout__
    
fa.replace_line('MAIN.py', 27, "recordToLog     = True                      # True, False")