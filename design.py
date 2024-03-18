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

def calc_C(fpc, Fy, tc, tw, tf, h, P=0, Ry=1, Rc=1):
    if   P < 0:
        print("Section is under tension.")    
    C = ((0.85 *Rc  *fpc *tc *tf +2 *tw *Ry *Fy *h -P)/
         (0.85 *Rc  *fpc *tc     +4 *tw *Ry *Fy))
    if C < 0:
        print(f"C = {C:.3f} < 0\nThe program exits in calc_C()!!!"); sys.exit()
    if C > h:
        print(f"C = {C:.3f} > h = {h:.3f}\nThe program exits in calc_C()!!!"); sys.exit()
    return C

def calc_Mp(fpc, Fy, tc, tw, tf, h, bf, P=0, Ry=1, Rc=1):
    y       = calc_C(fpc, Fy, tc, tw, tf, h, P, Ry, Rc)
    C1 = T1 = bf *tf            *Ry *Fy
    C2      = 2 *tw *(y -tf)    *Ry *Fy
    C3      = tc    *(y -tf)    *Rc *(0.85 *fpc)
    T2      = 2 *tw *(h -y -tf) *Ry *Fy

    Mp = ( C1 *( y -tf     /2) 
          +C2 *((y -tf)    /2) 
          +C3 *((y -tf)    /2) 
          +T1 *( h -y -tf  /2) 
          +T2 *((h -y -tf) /2) 
          -P *(h /2 -y))
    return Mp

#||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
#                   Step 1: Input Data
#||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
# Passed to inputData.py 



numSign = 65
#||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
#                   Step 2: Analysis for Design
#||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
fa.replace_line('MAIN.py', 31, "linearity       = True")
fa.replace_line('MAIN.py', 32, "typeBuild       = 'coupledWalls'            # 'CantileverColumn', 'coupledWalls', 'buildBeam', 'ShearCritBeam'")
fa.replace_line('MAIN.py', 43, "incrMono        = 0.5*((H_typical*n_story)/4000)")
fa.replace_line('MAIN.py', 46, "dispTarget      = drift*(H_typical*n_story)")
fa.replace_line('MAIN.py', 78, "plot_MomCurv    = False")
fa.replace_line('MAIN.py', 83, "Pu_1wall        = -load['wallG']")
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
    print(f"ID = {ID*100:.3f}%")


# Loads 

# Base Shear
Cvx     = fa.verDistFact(We, T1, h_1, h_typ, n_story)

OTM     = 0
for i in range(1, n_story+1):
    def h(n):
        if n == 1:
            return h_1
        else:
            return (h_1 + (n-1) *h_typ)
    OTM += h(i) * Cvx[i] *V_base

print(f"OTM = {OTM /1000:.1f} kN.m = {OTM/kip/inch:.1f} kip-in")
print(f"V_base = {V_base /1000:.1f} kN = {V_base /kip:.1f} kip")
# sys.exit()

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
C_CB        = calc_C(fpc, Fy, t_cCB, t_pwCB, t_pfCB, h_CB)
print(f"C_CB = {C_CB *1000:.1f} mm")
Mp_CB       = calc_Mp(fpc, Fy, tc, t_pwCB, t_pfCB, h_CB, bf)
Mn_CB       = Mp_CB # Assuming that the slenderness ratios in coupling beams are satisfied as per AISC
Mn_CB_Fi_b  = Fi_b *Mn_CB
print(f"Mn_CB*Fi_b = {Mn_CB_Fi_b /1000:.1f} kN.m")

# 3-3   Check Strength Ratios
"""'''''''''''''''''''''''"""
print(f"{'-'*numSign}\nCoupling Beams Strength Check\n{'-'*numSign}\n")
# a)    Shear Strength
R__V_CB     = Vu_CB /Vn_CB_Fi_v
print(f"===>>>Ratio = {R__V_CB:.2f}")
if R__V_CB >= 1.0:
    print("The Available Shear Strength of Coupling Beam is NOT SUFFICIENT!!!")
elif 0.95 < R__V_CB < 1.0:
    print("The Available Shear Strength of Coupling Beam is OK")
else:
    print("The Available Shear Strength of Coupling Beam is OK but NOT OPTIMUM!")

# b)    Flexural Strength
R__M_CB     = Mu_CB /Mn_CB_Fi_b
print(f"===>>>Ratio = {R__M_CB:.2f}")
if R__M_CB > 1.0:
    print("The Available Flexural Strength of Coupling Beam is NOT SUFFICIENT!!!")
elif 0.95 < R__M_CB <= 1.0:
    print("The Available Flexural Strength of Coupling Beam is OK")
else:
    print("The Available Flexural Strength of Coupling Beam is OK but NOT OPTIMUM!")


#||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
#                   Step 4: Design of Composite Walls
#||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

# 4-1   Calculate Required Strength
"""'''''''''''''''''''''''''''''"""
print(f"{'-'*numSign}\nCalculate Required Strength of Composite Walls\n{'-'*numSign}\n")
# a)    Required Axial Strength
Vn_Mp_exp   = 2 *(1.2 *M_exp) /L_CB
print(f"Vn_Mp_exp = {Vn_Mp_exp /1000:.1f} kN")
Pu_exp_CB   = n_story *Vn_Mp_exp
print(f"Pu_exp_CB = {Pu_exp_CB /1000:.1f} kN")
print(f"Pu_T = {Pu_T /1000:.1f} kN")
print(f"Pu_C = {Pu_C /1000:.1f} kN")

# b)    Required Shear Strength
V_amp   = 4 *V_base
print(f"V_amp = {V_amp /1000:.1f} kN")
nWalls  = 2
Vu_wall      = V_amp /nWalls
print(f"Vu_wall = {Vu_wall /1000:.1f} kN")

# c)    Required Flexural Strength for All Walls
gamma1      = (n_story *(1.2 *M_exp)) /(n_story *Mu_CB)
print(f"{gamma1 = }")
Mu_Both     = gamma1 *OTM -Pu_exp_CB *L_eff
print(f"Mu_Both = {Mu_Both /1000:.1f} kN.m")

# EIeff_Ten   = 7.7e9   *kip*inch **2
# EIeff_Com   = 1.81e10 *kip*inch **2
fa.replace_line('MAIN.py', 31, "linearity       = False")
fa.replace_line('MAIN.py', 32, "typeBuild       = 'CantileverColumn'            # 'CantileverColumn', 'coupledWalls', 'buildBeam', 'ShearCritBeam'")
fa.replace_line('MAIN.py', 43, "incrMono        = 0.5*((H_typical)/4000)")
fa.replace_line('MAIN.py', 46, "dispTarget      = drift*(H_typical)")
fa.replace_line('MAIN.py', 78, "plot_MomCurv    = True")

# Run Nonlinear Cross-ectional Analysis on Compression Wall
t_IEAna_i   = time.time()
print(f"{'='*numSign}\nNonlinear Cross-ectional Analysis on Compression Wall Started.\n{'='*numSign}\n")
fa.replace_line('MAIN.py', 83, "Pu_1wall        = Pu_C")
exec(open("MAIN.py").read()) 
EIeff_Com   = EIeff_walls[0]
t_IEAna_f   = time.time()
dur_IEAna   = (t_IEAna_f - t_IEAna_i)/60
print(f"{'='*numSign}\nNonlinear Cross-ectional Analysis on Compression Wall Finished in {dur_IEAna:.2f} mins.\n{'='*numSign}\n")

# Run Nonlinear Cross-ectional Analysis on Tension Wall
t_IEAna_i   = time.time()
print(f"{'='*numSign}\nNonlinear Cross-ectional Analysis on Tension Wall Started.\n{'='*numSign}\n")
fa.replace_line('MAIN.py', 83, "Pu_1wall        = Pu_T")
exec(open("MAIN.py").read()) 
EIeff_Ten   = EIeff_walls[0]
t_IEAna_f   = time.time()
dur_IEAna   = (t_IEAna_f - t_IEAna_i)/60
print(f"{'='*numSign}\nNonlinear Cross-ectional Analysis on Tension Wall  Finished in {dur_IEAna:.2f} mins.\n{'='*numSign}\n")
fa.replace_line('MAIN.py', 43, "incrMono        = 0.5*((H_typical*n_story)/4000)")
fa.replace_line('MAIN.py', 46, "dispTarget      = drift*(H_typical*n_story)")
fa.replace_line('MAIN.py', 83, "Pu_1wall        = -load['wallG']")

print(f"{'-'*numSign}\nCalculate Required Strength of Composite Walls\n{'-'*numSign}\n")
print(f"Vu_CB = {Vu_CB /1000:.1f} kN")
print(f"Vn_Mp_exp = {Vn_Mp_exp /1000:.1f} kN")
print(f"Pu_exp_CB = {Pu_exp_CB /1000:.1f} kN")
print(f"V_amp = {V_amp /1000:.1f} kN")
print(f"Vu_wall = {Vu_wall /1000:.1f} kN")
print(f"Mu_Both = {Mu_Both /1000:.1f} kN.m")

# The portion of the overturning moment resisted by the individual wall is
Mu_T        = (EIeff_Ten /(EIeff_Ten +EIeff_Com)) *Mu_Both
print(f"Mu_T = {Mu_T /1000:.1f} kN.m")
Mu_C        = (EIeff_Com /(EIeff_Ten +EIeff_Com)) *Mu_Both
print(f"Mu_C = {Mu_C /1000:.1f} kN.m")



# 4-2   Calculate Available Strength 
"""''''''''''''''''''''''''''''''"""
print(f"{'-'*numSign}\nCalculate Available Strength of Composite Walls\n{'-'*numSign}\n")
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
C_T         = calc_C(fpc, Fy, tc, tw, tf, Lw, Pu_T)
print(f"C_T = {C_T *1000:.1f} mm")
Mp_T        = calc_Mp(fpc, Fy, tc, tw, tf, Lw, bf, Pu_T)
Mn_T        = Mp_T
Mn_T_Fi_b   = Fi_b *Mn_T
print(f"Mn_T*Fi_b = {Mn_T_Fi_b /1000:.1f} kN.m")
    
# c.2)  Flexural Strength of Compression SpeedCore Wall
C_C         = calc_C(fpc, Fy, tc, tw, tf, Lw, Pu_C)
print(f"C_C = {C_C *1000:.1f} mm")
Mp_C        = calc_Mp(fpc, Fy, tc, tw, tf, Lw, bf, Pu_C)
Mn_C        = Mp_C
Mn_C_Fi_b   = Fi_b *Mn_C
print(f"Mn_C*Fi_b = {Mn_C_Fi_b /1000:.1f} kN.m")
print("\n")


# 4-3   Check Strength Ratios
"""'''''''''''''''''''''''"""
print(f"{'-'*numSign}\nCheck Strength Ratios of Composite Walls\n{'-'*numSign}\n")
# a)    Axial Strength
# a.1)  Tensile Strength
R__P_Twall  = Pu_T /Pn_T_Fi_t
print(f"\nR__P_Twall = {R__P_Twall *100:.1f}%")
if R__P_Twall > 1.0:
    print("The Available Tensile Strength of Tension Wall is NOT OK!!!")
elif 0.95 < R__P_Twall <= 1.0:
    print("The Available Tensile Strength of Tension Wall is OK")
else:
    print("The Available Tensile Strength of Tension Wall is OK, but NOT OPTIMUM!")


# a.2)  Compressive Strength
R__P_Cwall  = Pu_C /Pn_C_Fi_c
print(f"\nR__P_Cwall = {R__P_Cwall *100:.1f}%")
if R__P_Cwall > 1.0:
    print("The Available Compressive Strength of Compression Wall is NOT OK!!!")
elif 0.95 < R__P_Cwall <= 1.0:
    print("The Available Compressive Strength of Compression Wall is OK")
else:
    print("The Available Compressive Strength of Compression Wall is OK, but NOT OPTIMUM!")

# b)    Shear Strength
R__V_wall   = Vu_wall /Vn_Fi_v
print(f"\nR__V_wall = {R__V_wall *100:.1f}%")
if R__V_wall > 1.0:
    print("The Available Shear Strength of Walls is NOT OK!!!")
elif 0.95 < R__V_wall <= 1.0:
    print("The Available Shear Strength of Walls is OK")
else:
    print("The Available Shear Strength of Walls is OK, but NOT OPTIMUM!")

# c)    Flextural Strength
# c.1)  Flexural Strength of Tension SpeedCore Wall
R_M_Twall   = Mu_T/Mn_T_Fi_b
print(f"\nR_M_Twall = {R_M_Twall *100:.1f}%")
if R_M_Twall > 1.0:
    print("The Available Flexural Strength of Tension Wall is NOT OK!!!")
elif 0.95 < R_M_Twall <= 1.0:
    print("The Available Flexural Strength of Tension Wall is OK")
else:
    print("The Available Flexural Strength of Tension Wall is OK, but NOT OPTIMUM!")

# c.2)  Flexural Strength of Compression SpeedCore Wall
R__M_Cwall   = Mu_C/Mn_C_Fi_b
print(f"\nR__M_Cwall = {R__M_Cwall *100:.1f}%")
if R__M_Cwall > 1.0:
    print("The Available Flexural Strength of Compression Wall is NOT OK!!!")
elif 0.95 < R__M_Cwall <= 1.0:
    print("The Available Flexural Strength of Compression Wall is OK")
else:
    print("The Available Flexural Strength of Compression Wall is OK, but NOT OPTIMUM!")




#||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
#                   Coupling Ratio
#||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

M_couplingBeams     = L_eff * Pu_exp_CB
M_all               = Mu_Both + M_couplingBeams
R__coupling         = M_couplingBeams /M_all


#||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
#                   Results all at once
#||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
print(f"{'='*numSign}\n\t\t\tResults All at Once\n{'='*numSign}\n")
#______________________________________________________________________________
# 01) SFRS Performance
print(f"01) SFRS Performance\n{'_'*numSign}\n")

# 01.01) Allowable Interstory Drift of the Structure
print(f"ID              = {ID *100:.3f}%")

# 01.02) Base Shear 
print(f"V_base          = {V_base /1000:.1f} kN = {V_base /kip:.1f} kip")

# 01.03) Overturning Moment
print(f"OTM             = {OTM /1000:.1f} kN.m = {OTM/kip/inch:.1f} kip-in")

# 01.04) Coupling Ratio
print(f"Coupling Ratio  = {R__coupling *100:.1f}%")


#______________________________________________________________________________
# 02) Coupling Beam Demand/Capacity Info
print(f"\n\n02) Coupling Beam Demand/Capacity Info\n{'_'*numSign}\n")

# 02.01) Shear
print(f"Vu_CB           = {Vu_CB /1000:.1f} kN")
print(f"Vn_CB*Fi_v      = {Vn_CB_Fi_v /1000:.1f} kN")
print(f"===>>>DCR_V_CB  = {R__V_CB *100:.1f}%")
if R__V_CB >= 1.0:
    print("The Available Shear Strength of Coupling Beam is NOT SUFFICIENT!!!")
elif 0.95 < R__V_CB < 1.0:
    print("The Available Shear Strength of Coupling Beam is OK")
else:
    print("The Available Shear Strength of Coupling Beam is OK but NOT OPTIMUM!")
print(f"Vn_Mp_exp       = {Vn_Mp_exp /1000:.1f} kN")

# 02.02) Flexure
print(f"C_CB            = {C_CB *1000:.1f} mm")
print(f"Mu_CB           = {Mu_CB /1000:.1f} kN.m")
print(f"Mn_CB*Fi_b      = {Mn_CB_Fi_b /1000:.1f} kN.m")
print(f"===>>>DCR_M_CB  = {R__M_CB *100:.1f}%")
if R__M_CB > 1.0:
    print("The Available Flexural Strength of Coupling Beam is NOT SUFFICIENT!!!")
elif 0.95 < R__M_CB <= 1.0:
    print("The Available Flexural Strength of Coupling Beam is OK")
else:
    print("The Available Flexural Strength of Coupling Beam is OK but NOT OPTIMUM!")


#______________________________________________________________________________
# 03) Composite Walls Demand/Capacity Info
print(f"\n\n03) Composite Walls Demand/Capacity Info\n{'_'*numSign}\n")

# 03.01) Axial Force
print(f"Pu_T            = {Pu_T /1000:.1f} kN")
print(f"Pn_T*Fi_t       = {Pn_T_Fi_t /1000:.1f} kN")
print(f"===>>>DCR_Twall = {R__P_Twall *100:.1f}%")
if R__P_Twall > 1.0:
    print("The Available Tensile Strength of Tension Wall is NOT OK!!!")
elif 0.95 < R__P_Twall <= 1.0:
    print("The Available Tensile Strength of Tension Wall is OK")
else:
    print("The Available Tensile Strength of Tension Wall is OK, but NOT OPTIMUM!")

print(f"Pu_C            = {Pu_C /1000:.1f} kN")
print(f"Pn_C*Fi_c       = {Pn_C_Fi_c /1000:.1f} kN")
print(f"===>>>DCR_Cwall = {abs(R__P_Cwall) *100:.1f}%")
if R__P_Cwall > 1.0:
    print("The Available Compressive Strength of Compression Wall is NOT OK!!!")
elif 0.95 < R__P_Cwall <= 1.0:
    print("The Available Compressive Strength of Compression Wall is OK")
else:
    print("The Available Compressive Strength of Compression Wall is OK, but NOT OPTIMUM!")

# 03.02) Shear
print(f"V_amp           = {V_amp /1000:.1f} kN")
print(f"Vu_wall         = {Vu_wall /1000:.1f} kN")
print(f"Vn*Fi_v         = {Vn_Fi_v /1000:.1f} kN")
print(f"===>>>DCR_Vwall = {R__V_wall *100:.1f}%")
if R__V_wall > 1.0:
    print("The Available Shear Strength of Walls is NOT OK!!!")
elif 0.95 < R__V_wall <= 1.0:
    print("The Available Shear Strength of Walls is OK")
else:
    print("The Available Shear Strength of Walls is OK, but NOT OPTIMUM!")

# 03.03) Flexure
print(f"gamma1          = {gamma1:.3f}")
print(f"Mu_Both         = {Mu_Both /1000:.1f} kN.m")
print(f"EIeff_Com       = {EIeff_Com /1000:.1f} kN.m^2")
print(f"EIeff_Ten       = {EIeff_Ten /1000:.1f} kN.m^2")
print("\n")

print(f"Tension Wall     Share of Moment = {(EIeff_Ten /(EIeff_Ten +EIeff_Com))*100:.2f}%")
print(f"C_T             = {C_T *1000:.1f} mm")
print(f"Mu_T            = {Mu_T /1000:.1f} kN.m")
print(f"Mn_T*Fi_b       = {Mn_T_Fi_b /1000:.1f} kN.m")
print(f"===>>>DCR_MTwall= {R_M_Twall *100:.1f}%")
if R_M_Twall > 1.0:
    print("The Available Flexural Strength of Tension Wall is NOT OK!!!")
elif 0.95 < R_M_Twall <= 1.0:
    print("The Available Flexural Strength of Tension Wall is OK")
else:
    print("The Available Flexural Strength of Tension Wall is OK, but NOT OPTIMUM!")
print("\n")

print(f"Compression Wall Share of Moment = {(EIeff_Com /(EIeff_Ten +EIeff_Com))*100:.2f}%")
print(f"C_C             = {C_C *1000:.1f} mm")
print(f"Mu_C            = {Mu_C /1000:.1f} kN.m")
print(f"Mn_C*Fi_b       = {Mn_C_Fi_b /1000:.1f} kN.m")
print(f"===>>>DCR_MCwall= {R__M_Cwall *100:.1f}%")
if R__M_Cwall > 1.0:
    print("The Available Flexural Strength of Compression Wall is NOT OK!!!")
elif 0.95 < R__M_Cwall <= 1.0:
    print("The Available Flexural Strength of Compression Wall is OK")
else:
    print("The Available Flexural Strength of Compression Wall is OK, but NOT OPTIMUM!")




if recordToLogDesign == True:
    sys.stdout.close()
    sys.stdout = sys.__stdout__
    
fa.replace_line('MAIN.py', 27, "recordToLog     = True                      # True, False")





















