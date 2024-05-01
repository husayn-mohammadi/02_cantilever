import matplotlib.pyplot as plt
import os
import numpy             as np
# from scipy.interpolate import interp1d
# from matplotlib.animation import FuncAnimation
exec(open("MAIN.py").readlines()[18]) # It SHOULD read and execute exec(open("Input/units    .py").read())
unitForce   = "kN" if kN==1 else "N"   if N==1  else "kip" if kip==1 else "lb"  
unitLength  = "m"  if m==1  else "cm"  if cm==1 else "mm"  if mm==1  else "in." if inch==1 else "ft" 

def plotPushoverX(outputDir, types="monotonic"):
    
    disp    = np.loadtxt(f"{outputDir}/top_disp.txt", delimiter= ' ')
    # reac    = np.loadtxt(f"{outputDir}/reaction.txt", delimiter= ' ')
    
    x1       =  disp[:, 1]
    # Vx1      = -reac[:, 1]
    Vx1      =  disp[:, 0]
    
    x0      = np.array([0.])
    Vx0     = np.array([0.])
    
    x       = np.append(x0, x1)
    Vx      = np.append(Vx0, Vx1)
    
    x_Vx = np.column_stack((x, Vx))
    np.savetxt(f"{outputDir}/Pushover.txt", x_Vx)
    
    Vpeak       = max(Vx)
    V20peak     = 0.2 *Vpeak
    xAtV20      = interpolate(V20peak, Vx, x)
    stiffness   = V20peak /xAtV20
    xVpeak      = 1/stiffness *Vpeak
    
    fig, ax = plt.subplots(figsize=(10, 7), dpi=200)
    fig.suptitle(f"Pushover Curve: {outputDir[16:-4]}", fontsize=16)
    ax.set_xlabel(f'Displacement ({unitLength})')
    if unitForce=="kN":
        ax.set_ylabel('Shear (kN)')
        plt.plot(x, Vx, linewidth=0.8)
    elif unitForce=="N":
        ax.set_ylabel('Shear (kN)')
        plt.plot(x, Vx/1e3, linewidth=0.8)
    elif unitForce=="kip":
        ax.set_ylabel('Shear (kip)')
        plt.plot(x, Vx, linewidth=0.8)
    elif unitForce=="lb":
        ax.set_ylabel('Shear (kip)')
        plt.plot(x, Vx/1e3, linewidth=0.8)
    if types == "monotonic":
        plt.plot([xAtV20, xAtV20], [0,            V20peak*N/kN], 'r--')        # Vertical Line
        plt.plot([0,      xAtV20], [V20peak*N/kN, V20peak*N/kN], 'r--')   # Horizontal Line
        plt.plot([0,      xVpeak], [0,            Vpeak*N/kN],   'g--', label = f" stiffness = {stiffness/(kN*m):.1f} kN/m")
        plt.legend()
    plt.tight_layout()
    plt.show()
    
    return x, Vx

def plotStressStrain(outputDir,tagEleList): 
    listFiberMat = ['fiberSt', 'fiberCt1', 'fiberCt2'] 
    n       = len(tagEleList)
    zero    = np.array([0.])
    
    fig, ax = plt.subplots(3, n, figsize=(8*n, 16))
    fig.suptitle("Stress-Strain Curve", fontsize=16)
    if unitForce == "N" or unitForce == "kN":
        [a.set_ylabel('Stress (MPa)') for a in ax] if n==1 else [a.set_ylabel('Stress (MPa)') for a in ax[:,0]]
    elif unitForce == "kip" or unitForce == "lb":
        [a.set_ylabel('Stress (ksi)') for a in ax] if n==1 else [a.set_ylabel('Stress (ksi)') for a in ax[:,0]]
        
    SS      = {}; Stress = {}; Strain = {}
    for j, tagEle in enumerate(tagEleList):
        
        ax[2].set_xlabel('Strain') if n==1 else ax[2,j].set_xlabel('Strain')
        SS[j] = {}; Stress[j] = {}; Strain[j] = {}
        for i, fiberMat in enumerate(listFiberMat):
            for k, item in enumerate(['top', 'mid', 'bot']):
                if fiberMat == "fiberCt2" and item == "mid":
                    continue 
                SSS             = np.loadtxt(f"{outputDir}/{fiberMat}_{item}.txt", delimiter= ' ')
                Stress[j][k]    = np.append(zero, SSS[:,2*j+0])
                Strain[j][k]    = np.append(zero, SSS[:,2*j+1])
                SS[j][k]        = np.column_stack((Strain[j][k], Stress[j][k]))
                np.savetxt(f"{outputDir}/StressStrain{j}_{fiberMat}_{item}.txt", SS[j][k])
                            
            if unitForce == "N":
                ax[i].plot(Strain[j][0], Stress[j][0]/1e6, color='blue', label='top', linewidth=0.8) if n==1 else ax[i,j].plot(Strain[j][0], Stress[j][0]/1e6, color='blue', label='top', linewidth=0.8) 
                if fiberMat != "fiberCt2" and item != "mid":
                    ax[i].plot(Strain[j][1], Stress[j][1]/1e6, color='gray', label='mid', linewidth=0.8) if n==1 else ax[i,j].plot(Strain[j][1], Stress[j][1]/1e6, color='gray', label='mid', linewidth=0.8) 
                ax[i].plot(Strain[j][2], Stress[j][2]/1e6, color='red',  label='bot', linewidth=0.8) if n==1 else ax[i,j].plot(Strain[j][2], Stress[j][2]/1e6, color='red',  label='bot', linewidth=0.8) 
            elif unitForce == "kN":
                ax[i].plot(Strain[j][0], Stress[j][0]/1e3, color='blue', label='top', linewidth=0.8) if n==1 else ax[i,j].plot(Strain[j][0], Stress[j][0]/1e3, color='blue', label='top', linewidth=0.8)
                if fiberMat != "fiberCt2" and item != "mid":
                    ax[i].plot(Strain[j][1], Stress[j][1]/1e3, color='gray', label='mid', linewidth=0.8) if n==1 else ax[i,j].plot(Strain[j][1], Stress[j][1]/1e3, color='gray', label='mid', linewidth=0.8)
                ax[i].plot(Strain[j][2], Stress[j][2]/1e3, color='red',  label='bot', linewidth=0.8) if n==1 else ax[i,j].plot(Strain[j][2], Stress[j][2]/1e3, color='red',  label='bot', linewidth=0.8)
            elif unitForce == "lb":
                ax[i].plot(Strain[j][0], Stress[j][0]/1e3, color='blue', label='top', linewidth=0.8) if n==1 else ax[i,j].plot(Strain[j][0], Stress[j][0]/1e3, color='blue', label='top', linewidth=0.8)
                if fiberMat != "fiberCt2" and item != "mid":
                    ax[i].plot(Strain[j][1], Stress[j][1]/1e3, color='gray', label='mid', linewidth=0.8) if n==1 else ax[i,j].plot(Strain[j][1], Stress[j][1]/1e3, color='gray', label='mid', linewidth=0.8)
                ax[i].plot(Strain[j][2], Stress[j][2]/1e3, color='red',  label='bot', linewidth=0.8) if n==1 else ax[i,j].plot(Strain[j][2], Stress[j][2]/1e3, color='red',  label='bot', linewidth=0.8)
            elif unitForce == "kip":
                ax[i].plot(Strain[j][0], Stress[j][0], color='blue',     label='top', linewidth=0.8) if n==1 else ax[i,j].plot(Strain[j][0], Stress[j][0], color='blue',     label='top', linewidth=0.8) 
                if fiberMat != "fiberCt2" and item != "mid":
                    ax[i].plot(Strain[j][1], Stress[j][1], color='gray',     label='mid', linewidth=0.8) if n==1 else ax[i,j].plot(Strain[j][1], Stress[j][1], color='gray',     label='mid', linewidth=0.8) 
                ax[i].plot(Strain[j][2], Stress[j][2], color='red',      label='bot', linewidth=0.8) if n==1 else ax[i,j].plot(Strain[j][2], Stress[j][2], color='red',      label='bot', linewidth=0.8) 
            ax[i].set_title(f"tagEle: {tagEle} | Material: {fiberMat}") if n==1 else ax[i,j].set_title(f"tagEle: {tagEle} | Material: {fiberMat}")
            ax[i].legend() if n==1 else ax[i,j].legend()
    plt.tight_layout()
    plt.show()
    
def plotNTHA(H_typical, H_first, nFloors, outputDir, ta, tag, scaleFactor, SaGM, rec):
    os.makedirs(f"{outputDir}/{rec}", exist_ok=True)
    n       = nFloors -1
    fig, ax = plt.subplots(n+1, 1, figsize=(10, 5*n), dpi=100)
    ax[n].set_xlabel('time (s)')
    
    # Ground Motion
    gmt     = ta[:, 0]
    gma     = ta[:, 1] *scaleFactor
    igmaMax = np.argmax(np.abs(gma))
    gmaMax  = gma[igmaMax]
    gmtMax  = gmt[igmaMax]
    ax[0].set_ylabel('GMA [g]')
    ax[0].plot(gmt, gma, linewidth=0.8)
    ax[0].plot(gmtMax, gmaMax, color='red', marker='o', label=f"PGA = {abs(gmaMax):.4f}g\nS_T = {scaleFactor *SaGM:.4f} g")
    ax[0].legend()
    
    disp    = np.loadtxt(f"{outputDir}/disp{tag}.txt", delimiter= ' ')
    driftMAX= []
    for i in range(n):
        if i == 0:
            H = H_first
        else:
            H = H_typical
        
        # Structue Motion
        t           = np.array(disp[:, 0])
        drift       = np.array((disp[:, i+2] - disp[:, i+1])/H)
        iDriftMax   = np.argmax(np.abs(drift))
        tDriftMax   = t[iDriftMax]
        driftMax    = drift[iDriftMax]
        driftMAX.append(abs(driftMax))        
        # velo    = np.loadtxt(f"{outputDir}/velo{tag}.txt", delimiter= ' ')
        # acce    = np.loadtxt(f"{outputDir}/acce{tag}.txt", delimiter= ' ')
        # reac    = np.loadtxt(f"{outputDir}/R{tag}.txt", delimiter= ' ')
        
        # v       = velo[:, 1]
        # a       = np.array(acce[:, 1])
        # Vx      = reac[:, 1]; Vy = reac[:, 2]; Mz = reac[:, 3]
        
        # iaMax   = np.argmax(np.abs(a))
        # aMax    = a[iaMax]
        # taMax   = t[iaMax]
        
        ax[i+1].set_ylabel(f'drift[{i+1}]')
        ax[i+1].plot(t, drift, linewidth=0.8)
        ax[i+1].plot(tDriftMax, driftMax, color='red', marker='o', label=f"driftMax = {abs(driftMax) *100} %")
        ax[i+1].legend()
        
        # ax[1].set_ylabel('MMA (m/s^2)');    ax[1].plot(t, a, linewidth=0.8);        ax[1].plot(taMax, aMax, color='red', marker='o')
        # ax[1].set_ylabel('velocity     (m/s)');     ax[1].plot(t, v, linewidth=0.8)
        # ax[3].set_ylabel('reaction     (N)');       ax[3].plot(t, Vx, color='blue', label='Vx', linewidth=0.8); ax[3].plot(t, Vy, color='black', label='Vy', linewidth=0.8); ax[3].plot(t, Mz, color='red', label='Mz', linewidth=0.8); ax[3].legend()
    
    fig.suptitle(f"Dynamic Analysis Curves: {rec}", fontsize=16)
    plt.tight_layout()
    plt.savefig(f"{outputDir}/{rec}/GM-{rec}-{scaleFactor *SaGM:.5f}g.png")
    # plt.show()
    driftMax = max(driftMAX)
    return driftMax
    
def interpolate(Xp, x, y):

    results = []
    if len(x) > 2:
        for i in range(len(x) - 1):
            if x[i] <= Xp <= x[i + 1]:
                Yp = y[i] + (y[i + 1] - y[i]) * (Xp - x[i]) / (x[i + 1] - x[i])
                results.append(Yp)
            elif x[i] > x[i + 1] and x[i] >= Xp >= x[i + 1]:
                Yp = y[i] + (y[i + 1] - y[i]) * (Xp - x[i]) / (x[i + 1] - x[i])
                results.append(Yp)
                
    elif len(x) == 1:
        Yp = y[0]
        results.append(Yp)
    else:
        Yp = y[0] + (y[1] - y[0]) * (Xp - x[0]) / (x[1] - x[0])
        results.append(Yp)

    result_min = min(results)
    return result_min

def plotIDA(x, y, outputDirIDA, rec, mins, Threshold=False):
    os.makedirs(f"{outputDirIDA}/{rec[:-4]}", exist_ok=True)
    if Threshold == True:
        fig, ax = plt.subplots(1, 1, figsize=(10, 7), dpi=200)
        ax.set_xlabel('Drift (%)')
        ax.set_ylabel('Spectral Acceleration (g)')
        ax.scatter(x, y, color='red', marker='o', label=f"current drif={x[-1]} %\ncurrent S_T={y[-1]} g\nElapsed Time = {mins} min")
        ax.plot(x, y, linewidth=0.8)
        if len(x) >=2:
            Sa  = interpolate(5, x[-2:], y[-2:])
        else:
            Sa  = interpolate(5, x, y)
            
        plt.plot([5, 5], [0,  Sa], 'r--', linewidth=3)
        plt.plot([0, 5], [Sa, Sa], 'r--', label=f" Sa_5% = {Sa:.3f} g", linewidth=3)
        ax.legend()
        fig.suptitle(f"IDA Curve: {rec[:-4]}", fontsize=16)
        plt.tight_layout()
        plt.savefig(f"{outputDirIDA}/IDA-{rec[:-4]}-{Sa:.5f}g.png")
        plt.show()
        return Sa
    else:
        fig, ax = plt.subplots(1, 1, figsize=(10, 7), dpi=50)
        ax.set_xlabel('Drift (%)')
        ax.set_ylabel('Spectral Acceleration (g)')
        ax.scatter(x, y, color='red', marker='o', label=f"current drif={x[-1]} %\ncurrent S_T={y[-1]} g\nElapsed Time = {mins} min")
        ax.plot(x, y, linewidth=0.8)
        ax.legend()
        fig.suptitle(f"IDA Curve: {rec}", fontsize=16)
        plt.tight_layout()
        plt.savefig(f"{outputDirIDA}/{rec[:-4]}/IDA-{rec[:-4]}-{mins}.png")
        # plt.show()

def plotMomCurv(outputDir, tagEle, section, typeBuild):
    if typeBuild == "CantileverColumn":
        momenIndex = 1
    else:
        momenIndex = 2
    zero        = np.array([0.])
    # Store the section's Moment in an array
    momentTXT   = np.loadtxt(f"{outputDir}/moment{tagEle}.txt", delimiter= ' ')
    moment      = np.append(zero, momentTXT[:, momenIndex])
    Mpeak       = max(moment)
    Mpeak60perc = 0.6*Mpeak
    # Store the section's top and bottom strains in an array
    SSListTop   = np.loadtxt(f"{outputDir}/SS_top{tagEle}.txt", delimiter= ' ')
    SSListBot   = np.loadtxt(f"{outputDir}/SS_bot{tagEle}.txt", delimiter= ' ')
    # Calculate the curvature 
    StrainTop   = np.append(zero, SSListTop[:,1])
    StrainBot   = np.append(zero, SSListBot[:,1])
    h           = section.Hw + section.tw
    curvature   = ((StrainTop -StrainBot)
                   /h)
    # curAtM60per = np.interp(Mpeak60perc, moment, curvature)
    curAtM60per = interpolate(Mpeak60perc, moment, curvature)
    EI          = Mpeak60perc /curAtM60per
    curAtMpeakE = 1/EI *Mpeak
    MomCurv     = np.column_stack((curvature, moment))
    np.savetxt(f"{outputDir}/MomentCurvature_{tagEle}.txt", MomCurv)
    
    fig, ax = plt.subplots(figsize=(10, 7), dpi=200)
    fig.suptitle(f"Momemnt-Curvature: {tagEle}")
    ax.set_xlabel(f'Curvature (m^-1)')
    ax.set_ylabel('Moment (kN.m)')
    plt.plot(curvature, moment*N/kN, linewidth=0.8, label=f"M-C: Mpeak={Mpeak*N/kN:.1f} kN.m")
    plt.plot([curAtM60per, curAtM60per], [0, Mpeak60perc*N/kN], 'r--')        # Vertical Line
    plt.plot([0, curAtM60per], [Mpeak60perc*N/kN, Mpeak60perc*N/kN], 'r--')   # Horizontal Line
    plt.plot([0, curAtMpeakE], [0, Mpeak*N/kN], 'g--', label = f" EI = {EI/(kN*m**2):.1f} kN.m^2")                        # Horizontal Line
    plt.legend()
    plt.tight_layout()
    plt.show()
    
    return(EI)
    


















