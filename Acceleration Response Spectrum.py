import os
import functions.FuncAnalysis  as fa


g = 9.80665
outputDirIDA = "Output/IDA/02-story/Acceleration Response Spectrum"
os.makedirs(outputDirIDA, exist_ok=True); 

recList     = fa.get_file_names("Input/GM")
# recList     = recList[1:10]
list_S_CT   = []
for i_rec, rec in enumerate(recList):
    filePath    = f"Input/GM/{rec}" 
    T1          = 0.25
    dtGM        = float(rec[:7])
    NPTS        = int(rec[8:13])
    SaGM,PGV    = fa.get_spectral_acceleration(filePath, dtGM, T1, outputDirIDA)

    






    
    