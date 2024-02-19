import os, sys, time
import openseespy.opensees     as ops
import functions.FuncPlot      as fp
import functions.FuncAnalysis  as fa
# import numpy                   as np
# import matplotlib.pyplot       as plt


start_timeIDA = time.time()
sys.stdout = open('logIDA.txt', 'w') 
g = 9.80665
outputDirIDA = "Output/IDA"
os.makedirs(outputDirIDA, exist_ok=True)

recList     = fa.get_file_names("Input/GM")
recList     = recList[39:40]
# recList     = recList[5:10]
SaTarget    = 1.5 * 0.63 # MCE = 1.5*DBE
extraTime   = 0
numRecords  = len(recList)
list_S_CT   = []
for i_rec, rec in enumerate(recList):
    filePath    = f"Input/GM/{rec}" 
    T1          = 0.25
    # T1          = fa.analyzeEigen(1)
    dtGM        = float(rec[:7])
    NPTS        = int(rec[8:13])
    SaGM,PGV    = fa.get_spectral_acceleration(filePath, dtGM, T1, outputDirIDA)
    duration    = NPTS *dtGM +extraTime
    scaleFactor = SaTarget /SaGM # SaTarget is actually S_MT
    list_S_T    = []
    list_driftMax=[]
    tag         = 1
    t_begIDA1   = time.time()
    while True:
        print(f"\n{'#'*65}\nRunning record {i_rec+1:02}/{numRecords:02}: {rec} \nSa = {SaTarget:.3f}*g SF = {scaleFactor:.3f}\n{'#'*65}")
        ops.wipe()
        exec(open("MAIN.py").read())
        list_S_T.append(scaleFactor *SaGM)
        list_driftMax.append(abs(driftMax)*100)
        durIDA1     = time.time() - t_begIDA1; mins = int(durIDA1 /60)
        fp.plotIDA(list_driftMax, list_S_T, outputDirIDA, rec, mins)
        
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
            scaleFactor *= times
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
            scaleFactor *= times
            sign        = "b"
        tag += 1
    durIDA1     = time.time() - t_begIDA1; mins = int(durIDA1 /60)
    S_CT = fp.plotIDA(list_driftMax, list_S_T, outputDirIDA, rec, mins, True)
    list_S_CT.append([rec, S_CT])
    print(f"S_CT = {S_CT}")
print(f"list_S_CT = {list_S_CT}") 

finish_timeIDA  = time.time()
elapsedTime     = finish_timeIDA - start_timeIDA
mins = int(elapsedTime /60); secs = int(elapsedTime %60)
print("\n\n\n$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$")
print(f"IDA Finished in {mins}min+{secs}sec.")
print("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n\n\n")

sys.stdout.close()
sys.stdout = sys.__stdout__

    # scaleFactorList = [ 
    #                     20.*scaleFactor,
    #                     14.*scaleFactor,
    #                     13.*scaleFactor,
    #                     12.*scaleFactor,
    #                     11.*scaleFactor,
    #                     10.*scaleFactor,
    #                     # 9.0*scaleFactor,
    #                     # 8.0*scaleFactor,
    #                     # 7.0*scaleFactor,
    #                     # 6.0*scaleFactor,
    #                     5.0*scaleFactor,
    #                     # 4.0*scaleFactor,
    #                     # 3.0*scaleFactor,
    #                     # 2.0*scaleFactor,
    #                     # 1.5*scaleFactor,
    #                     # 1.0*scaleFactor,
    #                     # 0.5*scaleFactor,
    #                     # 0.2*scaleFactor,
    #                     0.1*scaleFactor, 
    #                     ]
    # for tag, scaleFactor in enumerate(scaleFactorList):
    #     print(f"\n{'#'*65}\nRunning record {i_rec+1:02}/{numRecords:02}: {rec} for Sa = {SaTarget}*g SF = {scaleFactor}\n{'#'*65}")
    #     ops.wipe()
    #     exec(open("MAIN.py").read())
    #     list_S_T.append(scaleFactor *SaTarget/g)
    #     list_driftMax.append(abs(driftMax)*100)
    #     fp.plotIDA(list_driftMax, list_S_T)
    # S_CT = fp.plotIDA(list_driftMax, list_S_T, True)































