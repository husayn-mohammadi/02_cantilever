import os, sys, time
import openseespy.opensees     as ops
import functions.FuncPlot      as fp
import functions.FuncAnalysis  as fa
import numpy                   as np
# import matplotlib.pyplot       as plt


#=============================================================================
#    Options
#=============================================================================
exec(open("Input/inputData.py").readlines()[16])  # Assigns SDC from input
# SDC     = "Dmax"  # "Dmax", "Dmin"
recordToLogIDA  = True

AlgorithmIDA    = "list__SF_CLP"       # Hunt_Fill, list__SF_CLP, manual__SF_CLP
SF_CLP          = 20                # If "manual__SF_CLP" is selected above

approach        = 3                 # "1" for scaling the Average RAS to MCE_AS at T1, 
                                    # "2" for scaling GMR RAS to MCE_AS at T1, 
                                    # "3" for scaling the Average RAS to MCE_AS between 0.2T1 to 1.5T1

extraTime       = 0

#=============================================================================
#    Pre-Settings for IDA
#=============================================================================
g = 9.80665
start_timeIDA = time.time()
fa.replace_line('MAIN.py', 27, "recordToLog     = False                      # True, False")
fa.replace_line('MAIN.py', 31, "linearity       = False")
fa.replace_line('MAIN.py', 32, "typeBuild       = 'coupledWalls'            # 'CantileverColumn', 'coupledWalls', 'buildBeam', 'ShearCritBeam'")
fa.replace_line('MAIN.py', 34, "typeAnalysis    = ['NTHA']             # 'monotonic', 'cyclic', 'NTHA'")
fa.replace_line('MAIN.py', 78, "plot_MomCurv    = False")
fa.replace_line('MAIN.py', 83, "Pu_1wall        = -load['wallG']")

recList         = fa.get_file_names("Input/GM")
RAS_average     = np.loadtxt("Input/GM/RAS_average.txt")
outputDirIDA    = "Output/IDA"; os.makedirs(outputDirIDA, exist_ok=True)


#=============================================================================
#    Start IDA
#=============================================================================
# Find T1
with open("MAIN.py") as file:
    lines = file.readlines()[:126]
code_to_exec = ''.join(lines)
exec(code_to_exec)
T1 = Periods[0]

# Record to log option check
if recordToLogIDA == True: sys.stdout = open('logIDA.txt', 'w') 

# Begin
recList     = recList[0:44]
numRecords  = len(recList)
list_S_CT   = []
for i_rec, rec in enumerate(recList):
    t_begIDA1   = time.time()
    filePath    = f"Input/GM/{rec}" 
    dtGM        = float(rec[3:10])
    NPTS        = int(rec[11:16])
    duration    = NPTS *dtGM +extraTime
    S_MT,SF_MCE = fa.get_spectral_acceleration(filePath, dtGM, T1, SDC, RAS_average, outputDirIDA, approach) #Here it should scale the RAS per SDC (Dmax or Dmin) and T1
    
    if AlgorithmIDA == "Hunt_Fill":
        SF_CLP          = 20
        tag             = 1
        list_driftMax   = []
        list_SCTtest    = []
        while True:
            S_CTtest = SF_CLP *S_MT
            list_SCTtest.append(S_CTtest)
            print(f"\n\n\n\n\n{'#'*65}")
            print(f"Running record {i_rec+1:02}/{numRecords:02}: {rec}")
            print(f"Sa = {S_MT:.3f}*g SF = {SF_CLP:.3f} ==> S_CTtest = {S_CTtest:.3f}")
            print(f"{'#'*65}\n\n")
            
            ops.wipe();  exec(open("MAIN.py").read())
            list_driftMax.append(100 *abs(driftMax))
            
            durIDA1     = time.time() - t_begIDA1; mins = int(durIDA1 /60)
            fp.plotIDA(list_driftMax, list_SCTtest, outputDirIDA, rec, mins)
            
            # Hunt & Fill Algorithm starts from here
            Times = (0.05/abs(driftMax))
            print(f"times = {Times}")
            
            if abs(driftMax) < 0.05:
                if "sign" not in globals():
                    sign    = "f"
                if sign == "b" and 0.04 <= abs(driftMax) <= 0.05:
                    break
                times = 0.3 *Times
                if times < 1.05:
                    times = 0.5 *Times
                if times < 1.05:
                    times = 0.75 *Times
                if times < 1.05:
                    times = 0.85 *Times
                if times < 1.05:
                    times = 0.95 *Times
                if times < 1.0:
                    times = Times
                SF_CLP *= times
                sign        = "f"
            else:
                if "sign" not in globals():
                    sign    = "b"
                if sign == "f" and 0.05 <= abs(driftMax) <= 0.06:
                    break
                times = 2 * Times
                if times > 0.95:
                    times = 1.5 *Times
                if times > 0.95:
                    times = 1.25 *Times
                if times > 0.95:
                    times = 1.15 *Times
                if times > 0.95:
                    times = 1.05 *Times
                if times > 1.0:
                    times = Times
                SF_CLP *= times
                sign        = "b"
            tag += 1
        # Now Interpolate/Extrapolate to get S_CT
        durIDA1 = time.time() - t_begIDA1; mins = int(durIDA1 /60)
        S_CT    = fp.plotIDA(list_driftMax, list_SCTtest, outputDirIDA, rec, mins, True); print(f"S_CT = {S_CT}")
        list_S_CT.append([rec, S_CT])
        
    elif AlgorithmIDA == "list__SF_CLP":
        SF_CLPList      = [ 
                            # 20.*SF_CLP,
                            # 15.*SF_CLP,
                            # 14.*SF_CLP,
                            # 13.*SF_CLP,
                            # 12.*SF_CLP,
                            # 11.*SF_CLP,
                            10.*SF_CLP,
                            9.0*SF_CLP,
                            8.0*SF_CLP,
                            7.0*SF_CLP,
                            6.0*SF_CLP,
                            5.0*SF_CLP,
                            4.0*SF_CLP,
                            3.0*SF_CLP,
                            2.0*SF_CLP,
                            1.5*SF_CLP,
                            1.0*SF_CLP,
                            # 0.5*SF_CLP,
                            # 0.2*SF_CLP,
                            0.1*SF_CLP, 
                            ]
        list_driftMax   = []
        list_SCTtest    = []
        for tag, SF_CLP in enumerate(SF_CLPList):
            S_CTtest = SF_CLP *S_MT
            list_SCTtest.append(S_CTtest)
            print(f"\n\n\n\n\n{'#'*65}")
            print(f"Running record {i_rec+1:02}/{numRecords:02}: {rec}")
            print(f"Sa = {S_MT:.3f}*g SF = {SF_CLP:.3f} ==> S_CTtest = {S_CTtest:.3f}")
            print(f"{'#'*65}\n\n")
            
            ops.wipe(); exec(open("MAIN.py").read())
            list_driftMax.append(100 *abs(driftMax))
            
            durIDA1     = time.time() - t_begIDA1; mins = int(durIDA1 /60)
            fp.plotIDA(list_driftMax, list_SCTtest, outputDirIDA, rec, mins)
            
        # Now Interpolate/Extrapolate to get S_CT
        durIDA1 = time.time() - t_begIDA1; mins = int(durIDA1 /60)
        S_CT    = fp.plotIDA(list_driftMax, list_SCTtest, outputDirIDA, rec, mins, True); print(f"S_CT = {S_CT}")
        list_S_CT.append([rec, S_CT])
        
    elif AlgorithmIDA == "manual__SF_CLP":
        if "list_driftMax" not in globals(): list_driftMax = []
        if "list_SCTtest"  not in globals(): list_SCTtest  = []
        S_CTtest = SF_CLP *S_MT
        list_SCTtest.append(S_CTtest)
        print(f"\n\n\n\n\n{'#'*65}")
        print(f"Running record {i_rec+1:02}/{numRecords:02}: {rec}")
        print(f"Sa = {S_MT:.3f}*g SF = {SF_CLP:.3f} ==> S_CTtest = {S_CTtest:.3f}")
        print(f"{'#'*65}\n\n")
        
        ops.wipe(); exec(open("MAIN.py").read())
        list_driftMax.append(100 *abs(driftMax))
        
        durIDA1     = time.time() - t_begIDA1; mins = int(durIDA1 /60)
        fp.plotIDA(list_driftMax, list_SCTtest, outputDirIDA, rec, mins)
        
        
print(f"list_S_CT = {list_S_CT}") 

finish_timeIDA  = time.time()
elapsedTime     = finish_timeIDA - start_timeIDA
mins = int(elapsedTime /60); secs = int(elapsedTime %60)
print("\n\n\n$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$")
print(f"IDA Finished in {mins}min+{secs}sec.")
print("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n\n\n")

if recordToLogIDA == True:
    sys.stdout.close()
    sys.stdout = sys.__stdout__




fa.replace_line('MAIN.py', 27, "recordToLog     = True                      # True, False")




























