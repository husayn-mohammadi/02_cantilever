"C-PSW/CF Section: Shafaei PP=136"
# import sys
# from functions.ClassComposite import compo
exec(open("MAIN.py").readlines()[18]) # It SHOULD read and execute exec(open(f"Input/units{'US'}.py").read())
# exec(open("../Input/unitsSI.py").read()) # It SHOULD read and execute exec(open("Input/units    .py").read())
# exec(open("MAIN.py").readlines()[20]) # It SHOULD read and execute exec(open("Input/materialParameters.py").read())

# Seismic Coefficients
Cd      = 5.5
Ie      = 1.
R       = 8
Omega0  = 2.5
Rho     = 1.

# Composite wall resistance factor
Fi_v = Fi_b = Fi_c = Fi_t = 0.9

# Material
Es      = 29000 *ksi
Gs      = 11200 *ksi
Fy      = 50    *ksi
Fu      = 65    *ksi
Ry      = 1.1

fpc     = 6     *ksi
Ec      = 4500  *ksi
Gc      = 1800  *ksi
Rc      = 1.3

linearity = 1

#=============================================================================
#    Elements
#=============================================================================
Hw          = 132   *inch
H_CB        = 24    *inch
bf          = 24    *inch
tw          = 9/16  *inch
RhoW        = 2 *tw /bf
tf          = tw
tc          = bf - 2*tw
t_sc        = tc + 2*tw
btie        = 12    *inch
Stie        = 12    *inch
dtie        = 1     *inch
lsr         = btie/tw
t_pfCB      = 0.5   *inch
t_pwCB      = 0.5   *inch

b           = 114*mm
NfibeY      = 10

Section = {
    'wall': { # C-PSW/CF Wall Section
        #tags       = [tagSec, tagMatStFlange, tagMatStWeb, tagMatCtUnconf, tagMatCtConf]
        'tags'      : [1,      1,              2,           3,              4           ],
        #propStPart = [B,      H,         Es,      Fy,      Fu,      eps_sh, eps_ult, nu,   alpha, beta, gamma, Cf,  a1,  limit] 
        'propWeb'   : [tw,     Hw,        Es,      Fy,      Fu,      0.007,  0.12,    0.28, 0.65,  1.0,  1.0,   0.5, 4.3, 0.01],
        'propFlange': [bf,     tf,        Es,      Fy,      Fu,      0.007,  0.12,    0.28, 0.65,  1.0,  1.0,   0.5, 4.3, 0.01],
        #propCore   = [tc,     fpc,       wc,      lamConf, lamUnconf]
        'propCore'  : [tc,     fpc,       0.2*mm,  0.05,     0.25    ]
    },
    'beam': { # Composite Beam Section
        #tags       = [tagSec, tagMatStFlange, tagMatStWeb, tagMatCtUnconf, tagMatCtConf]
        'tags'      : [2,      5,              6,           7,              8           ],
        #propStPart = [B,      H,         Es,      Fy,      Fu,      eps_sh, eps_ult, nu,   alpha, beta, gamma, Cf,  a1,  limit] 
        'propWeb'   : [t_pwCB, H_CB,      Es,      Fy,      Fu,      0.007,  0.12,    0.28, 0.65,  1.0,  1.0,   0.5, 4.3, 0.01],
        'propFlange': [bf,     t_pfCB,    Es,      Fy,      Fu,      0.007,  0.12,    0.28, 0.65,  1.0,  1.0,   0.5, 4.3, 0.01],
        #propCore   = [tc,     fpc,       wc,      lamConf, lamUnconf]
        'propCore'  : [tc,     fpc,       0.2*mm,  0.05,     0.25     ]
    },
    }

#=============================================================================
#    Frame Data:
#=============================================================================
n_story         = 2
H_typical       = 14    *ft
H_first         = 17    *ft
LDR_CB          = 4
L_CB            = LDR_CB * H_CB
L_Bay           = Hw + L_CB #(Hw+2*tf) + L_CB
H_story_List    = [H_first, *((n_story-1)*[H_typical])]       # [Hstory1, *((numStories-1)*[HstoryTypical])]
L_Bay_List      = 2*[L_Bay]#, 5.*m, 5.*m, 5.*m]        # [*LBays]

L               = H_typical

# Building Geometry
Lf              = 200   *ft
Wf              = 120   *ft
h_typ           = H_typical
h_1             = H_first
A_SFRS          = Lf *Wf /2
#=============================================================================
#    Loading
#=============================================================================
#   Cantilever Loads
Pno             = 12776238.963599999 *N
# or:
ALR             = 0.02  # Axial Load Ratio
Py              = ALR * Pno
#   Frame Loads
load={}
DL_Floor        = 0.12  *ksf #90 *psf
DL_PWalls       = 0 #25 *psf
LL_Floor        = 0 #50 *psf
LL_Roof         = 0 #20 *psf

##  Tributary Loading
L_Bay_y         = (42.5 +35) *ft
L_Bay_x         = (42.5 +30) /2 *ft
A_Tributary     = 0.5*L_Bay_y * L_Bay_x
DL_Tributary    = A_Tributary * DL_Floor
LL_Tributary    = A_Tributary * LL_Floor
load["wall"]    = 1.0*DL_Tributary + 0.25*LL_Tributary
# LoadG           = 72
# load["wall"]    = LoadG * kip

##  Loading the Leaning Columns
n_Bay_x         = len(L_Bay_List)
# A_SFRS          = (1.5 * L_Bay_y) * ((n_Bay_x+1) * L_Bay_x)
A_Leaning       = A_SFRS - A_Tributary*n_Bay_x
L_PWall         = L_Bay_y + ((n_Bay_x+1) * L_Bay_x) - n_Bay_x*Hw
DL_Leaning      = A_Leaning * DL_Floor + L_PWall*H_typical * DL_PWalls
LL_Leaning      = A_Leaning * LL_Floor
load["leaningColumn"] = 1.0*DL_Leaning + 0.25*LL_Leaning
# load["leaningColumn"] = 0 * kip


#||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
#                   Step 1: Input Data
#||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
# Floor Loads
DL      = DL_Floor +DL_PWalls
LL      = LL_Floor +LL_Roof

We      = (1.0 *DL  +0.25 *LL) *A_SFRS *n_story

# Wall and Coupling Beam Sections
Lw      = Hw
b_CB    = bf
h_CB    = H_CB


# Coupling Beam Properties
As_CB   = 2 *t_pwCB *h_CB + 2 *t_pfCB *(b_CB -2 *t_pwCB)
Ac_CB   = (b_CB -2 *t_pwCB) *(h_CB -2 *t_pfCB)
Asw_CB  = 2 *h_CB *t_pwCB
# Is_CB   = 2 *(1/12 *t_pwCB *h_CB **3 + 1/12 *(b_CB -2 *t_pwCB) *t_pfCB **3 +
#             t_pfCB *(b_CB -2 *t_pwCB) *((h_CB -t_pfCB) /2) **2)
# Ic_CB   = 1/12 *(b_CB -2 *t_pwCB) *(h_CB -2 *t_pfCB) **3

# EAuncrCB= 0.8 *(Es *As_CB +Ec *Ac_CB)
# C3      = min(0.9, 0.45 +3 *(As_CB /b_CB/ h_CB))
# EIeff_CB= 0.64 *(Es *Is_CB +C3 *Ec *Ic_CB)
# GAv_CB  = Gs *Asw_CB + Gc *Ac_CB

# Planar SpeedCore Wall Properties
As      = tw *(2 *(Lw -2 *tw) + 2 *t_sc)
Ac      = Lw *t_sc -As
Asw     = 2 *Lw *tw
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
rho_min  = 0.01
rho_max  = 0.1

# 1.1.  Coupling Beams
As_CBmin = rho_min *b_CB *h_CB
As_CBmax = rho_max *b_CB *h_CB
if As_CB < As_CBmin: 
    print("Area of Steel is less than minimum for Coupling Beams!!!")
elif As_CB > As_CBmax:
    print("Area of Steel is greater than maximum for Coupling Beams!!!")

# 1.2.  SpeedCore Walls
Asmin = rho_min *t_sc *Lw
Asmax = rho_max *t_sc *Lw
if As < Asmin: 
    print("Area of Steel is less than minimum for Composite Walls!!!")
elif As > Asmax:
    print("Area of Steel is greater than maximum for Composite Walls!!!")


# 2. Check Plate Slenderness Ratios
"""'''''''''''''''''''''''''''''"""
# 1.1.  Coupling Beams
lambdaP_fCB     = 2.37 *(Es /Ry /Fy) **0.5
lambdaP_wCB     = 2.66 *(Es /Ry /Fy) **0.5

b_cCB           = b_CB -2 *t_pwCB
lambda_fCB      = b_cCB /t_pfCB
R_lambda_fCB    = lambda_fCB/lambdaP_fCB

h_cCB           = h_CB -2 *t_pfCB # Clear height of the web plate
lambda_wCB      = h_cCB /t_pwCB
R_lambda_wCB    = lambda_wCB/lambdaP_wCB

if R_lambda_fCB > 1.0: 
    print(f"Coupling Beam Flange Plate Slenderness OVERRATED!!!\n===>>>Lam_fCB/LamP_fCB = {R_lambda_fCB:.2f}")
if R_lambda_wCB > 1.0: 
    print(   f"Coupling Beam Web Plate Slenderness OVERRATED!!!\n===>>>Lam_wCB/LamP_wCB = {R_lambda_wCB:.2f}")

# 1.2.  SpeedCore Walls
lambda_btie     = btie/tw
lambdaP_btie    = 1.05 *(Es /Ry /Fy) **0.5
R_lambda_btie   = lambda_btie/lambdaP_btie

lambdaP_btieTop = 1.2 *(Es /Fy) **0.5
R_lambda_btieTop= lambda_btie/lambdaP_btieTop

if R_lambda_btie > 1.0:    
    print(f"Wall Faceplate Slenderness in flexural yielding zone not CHECKED!!!\n===>>>Lam_fW/LamP_fW = {R_lambda_btie:.2f}") 
if R_lambda_btieTop > 1.0: 
    print(f"Wall Faceplate Slenderness above flexural yielding zone not CHECKED!!!\n===>>>Lam_fW/LamP_fW = {R_lambda_btieTop:.2f}") 

Alpha           = 1.7 *(t_sc /tw -2) *(tw /dtie) **4
lambdaP_tbs     = (Es /(2 *Alpha +1)) **0.5
lambda_tbs      = Stie /tw
R_lambda_tbs    = lambda_tbs/lambdaP_tbs
if R_lambda_tbs > 1.0: 
    print(f"Wall Tie bar spacing not CHECKED!!!\n===>>>Lam_tbs/LamP_tbs = {R_lambda_tbs:.2f}") 


# 3. Check Shear Criticality of Coupling Beams
"""''''''''''''''''''''''''''''''''''''''''"""
t_cCB   = b_CB -2 *t_pwCB

C_CBexp = ((2 *t_pwCB *Ry *Fy*h_CB  +0.85 *Rc *fpc *t_cCB *t_pfCB)/
           (4 *t_pwCB *Ry *Fy       +0.85 *Rc *fpc *t_cCB))
print(f"C_CBexp = {C_CBexp *1000:.1f} mm")
C_1exp  = (b_CB -2 *t_pwCB) *t_pfCB *Ry *Fy
C_2exp  = 2 *t_pwCB *C_CBexp *Ry *Fy
C_3exp  = 0.85 *Rc *fpc *t_cCB *(C_CBexp -t_pfCB)
T_1exp  = (b_CB -2 *t_pwCB) *t_pfCB *Ry *Fy
T_2exp  = 2 *t_pwCB *(h_CB -C_CBexp) *Ry *Fy
M_exp = (C_1exp *(C_CBexp -t_pfCB /2) +
            C_2exp *(C_CBexp /2) +
            C_3exp *((C_CBexp -t_pfCB) /2) +
            T_1exp *(h_CB -C_CBexp -t_pfCB /2) +
            T_2exp *((h_CB -C_CBexp) /2))
V_exp    = 0.6 *Ry *Fy *Asw_CB +0.06 *Ac_CB *(Rc *fpc /MPa) **0.5 #!!! Take care NOT to put fpc in other units except MPa

print(f"M_exp = {M_exp /1000:.1f} kN.m")
print(f"V_exp = {V_exp /1000:.1f} kN")

# Check Flexure-Criticality Condition
if V_exp *L_CB /M_exp >= 2.4: 
    print("Flexure-Criticality Confirmed!")



