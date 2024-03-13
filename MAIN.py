import time, os, sys, winsound
import openseespy.opensees     as ops
import opsvis                  as opv
import functions.FuncModel     as fm
import functions.FuncAnalysis  as fa
import functions.FuncRecorders as fr
import functions.FuncPlot      as fp
import numpy                   as np






#=============================================================================
#    Input File
#=============================================================================

exec(open("Input/unitsSI.py").read())       # This determines the OUTPUT units: unitsUS.py/unitsSI.py
exec(open('Input/inputData.py').read())
# exec(open("Input/materialParameters.py").read())
ops.logFile("logOpenSEES.txt")
#=============================================================================
#    Define Variables
#=============================================================================
# Modeling Options
recordToLog     = True                      # True, False
modelFoundation = True
rotSpring       = True
exertGravityLoad= True
linearity       = False
typeBuild       = 'CantileverColumn'            # 'CantileverColumn', 'coupledWalls', 'buildBeam', 'ShearCritBeam'
typeCB          = 'discritizedBothEnds'     # 'discretizedAllFiber', 'FSF', 'FSW', discritizedBothEnds (FSF = FlexureShearFlexure, FSW = FlexureShearWall)
typeAnalysis    = ['monotonic']             # 'monotonic', 'cyclic', 'NTHA'

Lw              = Section['wall']['propWeb'][1] + 2*Section['wall']['propFlange'][1]
PHL_wall        = 2/3 * Section['wall']['propWeb'][1]
PHL_beam        = 2/3 * Section['beam']['propWeb'][1]
numSegWall      = 7                         # If numSegWall=0, the model will be built only with one linear elastic element connecting the base node to top node
numSegBeam      = 15
SBL             = 0.3 *m                    # Length of Shear Link (Shear Beam)
# Monotonic Pushover Analysis
incrMono        = 0.5*((H_typical*n_story)/4000)
numIncrInit     = 5
drift           = 0.01
dispTarget      = drift*(H_typical*n_story)
# Cyclic Pushover Analysis
incrCycl        = incrMono
dY              = 9 *mm
CPD1            = 1                         # CPD = cyclesPerDisp; which should be an integer
CPD2            = 1


dispTarList     = [ 
                    *(CPD1*[dY/3]), 
                    *(CPD1*[2/3*dY]), 
                    *(CPD1*[dY]),   
                    *(CPD1*[1.5*dY]), 
                    *(CPD1*[2*dY]),
                    *(CPD1*[3*dY]), 
                    *(CPD1*[4*dY]),   
                    *(CPD1*[5*dY]), 
                    *(CPD2*[6*dY]),   
                    *(CPD2*[7*dY]),
                    *(CPD2*[8*dY]), 
                    *(CPD2*[9*dY]),   
                    *(CPD2*[10*dY])
                    ]


# Plotting Options:
buildingWidth1=20.; buildingHeight1=17.
plot_undefo     = True
plot_loaded     = True
plot_defo       = True
sfac            = 10

plot_MomCurv    = True
plot_Analysis   = True
plot_StressStrain=False
plot_section    = False
typeSpring      = "elastic"  # "elastic", "IMK_Pinching"
Pu_1wall        = -load['wallG']
# Pu_1wall        = -load['wallG']
#=============================================================================
#    MAIN
#=============================================================================
start_time = time.time()
if recordToLog == True:
    sys.stdout = open('logMAIN.txt', 'w')    

numFolder = 1
for types in typeAnalysis:
    outputDir       = f'Output/Pushover/{types}/{numFolder:03}'
    outputDirWalls  = f'Output/Pushover/{types}/{numFolder:03}/wall'
    outputDirBeams  = f'Output/Pushover/{types}/{numFolder:03}/beams'
    outputDirNTHA   = "Output/NTHA"
    
    os.makedirs(outputDir,      exist_ok=True)
    os.makedirs(outputDirWalls, exist_ok=True)
    os.makedirs(outputDirBeams, exist_ok=True)
    os.makedirs(outputDirNTHA,  exist_ok=True)
    
    # Build Model
    ops.wipe()
    ops.model('basic', '-ndm', 2, '-ndf', 3)
            
    if typeBuild == "CantileverColumn":
        Py = 1
        tagNodeControl, tagNodeBase, tagEleListToRecord_wall, wall = fm.buildCantileverN(L, Py, PHL_wall, numSegWall, "wall", modelFoundation, linearity, typeSpring)
        fa.analyzeEigen(1)
    elif typeBuild == "CantileverBeam":
        Py = 1
        tagNodeControl, tagNodeBase, tagEleListToRecord_wall, wall = fm.buildCantileverN(L, Py, PHL_beam, numSegBeam, "beam", modelFoundation, linearity, typeSpring)
        fa.analyzeEigen(1)
    elif typeBuild == 'buildBeam':
        tagNodeControl, tagNodeBase, tagEleListToRecord_wall, wall = fm.buildBeam(L, PHL_beam, numSegBeam, rotSpring, linearity, typeSpring)
    elif typeBuild == 'coupledWalls':
        P = n_story * load['wall']
        tagNodeControl, tagNodeBase, buildingWidth, buildingHeight, coords, wall, tagEleListToRecord_wall, beam, tagEleListToRecord_beam, tagNodeLoad = fm.coupledWalls(H_story_List, L_Bay_List, Lw, P, load, numSegBeam, numSegWall, PHL_wall, PHL_beam, SBL, typeCB, plot_section, modelFoundation, rotSpring, linearity, typeSpring)
        # fa.analyzeEigen(n_story, True)
    else:
        tagNodeControl, tagNodeBase  = fm.buildShearCritBeam(L)
        
    # Plot Model
    if plot_undefo == True:
        if "buildingWidth" not in globals(): buildingWidth=buildingHeight = 10
        opv.plot_model(node_labels=1, element_labels=1, fig_wi_he=(buildingWidth+buildingWidth1, buildingHeight+buildingHeight1),
                       fmt_model={'color': 'blue', 'linestyle': 'solid', 'linewidth': 0.6, 'marker': '.', 'markersize': 3})
    
    # Run Gravity Analysis
    if exertGravityLoad == True:
        if typeBuild == 'coupledWalls':
            fa.gravity(load, tagNodeLoad)
        elif typeBuild == 'CantileverColumn':
            # Axial Force Capacity of Walls (Pno)
            # Pno = wall.Pno
            # Pno = 0
            # fa.gravity(ALR*Pno, tagNodeControl)
            fa.gravity(Pu_1wall, tagNodeControl)
    
    # Record Lateral Loading Analysis Results
    fr.recordPushover(tagNodeControl, tagNodeBase, outputDir)
    Hw = wall.Hw; tf = wall.tf; Hc2 = wall.Hc2
    if typeBuild == "CantileverColumn" or typeBuild == "buildBeam":
        fr.recordStressStrain(outputDirWalls, tagEleListToRecord_wall,  wall)
    elif typeBuild == "coupledWalls":
        fr.recordStressStrain(outputDirWalls, tagEleListToRecord_wall,  wall)
        fr.recordStressStrain(outputDirBeams, tagEleListToRecord_beam,  beam)
        
    for tagEle in tagEleListToRecord_wall:
        tagNode = ops.eleNodes(tagEle)
        print(f"tagEleListToRecord_wall = {tagEleListToRecord_wall}")
        print(f"tagEle = {tagEle}")
        print(f"tagNode = {tagNode}")
        fr.recordMomCurv(tagNode, tagEle, wall, outputDirWalls)
    
    # Run Lateral Loading Analysis Results 
    if types == 'monotonic':
        # if True:
        # if False:
        if linearity == False:
            start_time_monotonic = time.time()
            print("\n\n\n$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$")
            print(f"Monotonic Pushover Analysis Initiated at {(start_time_monotonic - start_time):.0f}sec.")
            print("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n\n\n")
            if typeBuild == 'coupledWalls':
                fa.pushoverDCF(dispTarget, incrMono, numIncrInit, tagNodeControl, tagNodeLoad['wall'])
            else:
                fa.pushoverDCF(dispTarget, incrMono, numIncrInit, tagNodeControl, tagNodeControl, distributeOnWalls=False)
            
            finish_time_monotonic = time.time()
            mins = int((finish_time_monotonic - start_time_monotonic)/60)
            secs = int((finish_time_monotonic - start_time_monotonic)%60)
            print("\n\n\n$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$")
            print(f"Monotonic Pushover Analysis Finished in {mins}min+{secs}sec.")
            print("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n\n\n")
        
        else:
        # if linearity == True:
        # if False:
            # tagNodeLoadH = []
            # if type(tagNodeLoad) == dict:
            #     for tagNode in tagNodeLoad['wall']:
            #         tagCoordX = f"{tagNode}"[-3:-1]
            #         if tagCoordX == '00':
            #             tagNodeLoadH.append(tagNode)
            # elif type(tagNodeLoad) == int:
            #     tagNodeLoadH = tagNodeLoad
            plot_Analysis       = False
            plot_StressStrain   = False
            T1, driftMaximum, V_base, Vr_CB = fa.pushoverLCF(tagNodeControl, tagNodeBase, tagEleListToRecord_beam)
            print(f"driftMax = {driftMaximum*100:.5f}%")
        
        if plot_loaded == True:
            if "buildingWidth" not in globals(): buildingWidth=buildingHeight = 10
            opv.plot_loads_2d(nep=17, sfac=False, fig_wi_he=(buildingWidth+buildingWidth1, buildingHeight+buildingHeight1), fig_lbrt=False, fmt_model_loads={'color': 'black', 'linestyle': 'solid', 'linewidth': 1.2, 'marker': '', 'markersize': 1}, node_supports=True, truss_node_offset=0, ax=False)
        if plot_defo == True:
            if "buildingWidth" not in globals(): buildingWidth=buildingHeight = 10
            sfac = opv.plot_defo(fig_wi_he=(buildingWidth+buildingWidth1, buildingHeight+buildingHeight1),
                                 #fmt_defo={'color': 'blue', 'linestyle': 'solid', 'linewidth': 0.6, 'marker': '.', 'markersize': 3}
                                 )
            # opv.plot_defo(sfac)
    elif types == 'cyclic':
        start_time_cyclic = time.time()
        print("\n\n\n$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$")
        print(f"Cyclic Pushover Analysis Initiated at {(time.time() - start_time):.0f}sec.")
        print("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n\n\n")
        fa.cyclicAnalysis(dispTarList, incrCycl, tagNodeControl, numIncrInit)
        finish_time_cyclic = time.time()
        mins = int((finish_time_cyclic - start_time_cyclic)/60)
        secs = int((finish_time_cyclic - start_time_cyclic)%60)
        print("\n$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$")
        print(f"Cyclic Pushover Analysis Finished in {mins}min+{secs}sec.")
        print("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n\n\n")
        if plot_loaded == True:
            opv.plot_loads_2d(nep=17, sfac=False, fig_wi_he=(buildingWidth+buildingWidth1, buildingHeight+buildingHeight1), fig_lbrt=False, fmt_model_loads={'color': 'black', 'linestyle': 'solid', 'linewidth': 1.2, 'marker': '', 'markersize': 1}, node_supports=True, truss_node_offset=0, ax=False)
        if plot_defo == True:
            sfac = opv.plot_defo(fig_wi_he=(buildingWidth+buildingWidth1, buildingHeight+buildingHeight1),
                                 #fmt_defo={'color': 'blue', 'linestyle': 'solid', 'linewidth': 0.6, 'marker': '.', 'markersize': 3}
                                 )
    elif types == 'NTHA':
        plot_Analysis   = False
        plot_StressStrain=False
        start_time_NTHA = time.time()
        print("\n\n\n$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$")
        print(f"NTHA Initiated at {(time.time() - start_time):.0f}sec.")
        print("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n\n\n")
        if 'filePath' not in globals():
            filePath    = "Input/GM/0.02000_01500_RSN825_CAPEMEND_CPM000.txt"
            rec         = filePath[9:]
            SaGM        = 1.034481507651492 # for example
            scaleFactor = 8.96144583680962
            dtGM        = float(rec[:7])
            NPTS        = int(rec[8:13])
            # NPTS        = 20
            duration    = dtGM *NPTS
            tag         = 1
            outputDirIDA= "Output/IDA"
        a  = fa.read_ground_motion_record(filePath); t = np.array([i *dtGM for i in range(NPTS)])
        ta = np.column_stack((t, a[:NPTS]))
        fr.recordDataNTHA(tagNodeBase, tagNodeControl, outputDirNTHA, tag)
        fa.NTHA1(tagNodeControl, tagNodeBase, filePath, scaleFactor, dtGM, NPTS, duration, tag)
        if type(tagNodeControl) != list:
            tagNodeControl = [tagNodeControl]
        nFloors = len(tagNodeControl)
        driftMax = fp.plotNTHA(H_typical, H_first, nFloors, outputDirNTHA, ta, tag, scaleFactor, SaGM, rec[:-4])
        finish_time_NTHA = time.time()
        mins = int((finish_time_NTHA - start_time_NTHA)/60)
        secs = int((finish_time_NTHA - start_time_NTHA)%60)
        print("\n$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$")
        print(f"NTHA Finished in {mins}min+{secs}sec.")
        print("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n\n\n")
        if plot_loaded == True:
            opv.plot_loads_2d(nep=17, sfac=False, fig_wi_he=(buildingWidth+buildingWidth1, buildingHeight+buildingHeight1), fig_lbrt=False, fmt_model_loads={'color': 'black', 'linestyle': 'solid', 'linewidth': 1.2, 'marker': '', 'markersize': 1}, node_supports=True, truss_node_offset=0, ax=False)
        if plot_defo == True:
            sfac = opv.plot_defo(fig_wi_he=(buildingWidth+buildingWidth1, buildingHeight+buildingHeight1),
                                 #fmt_defo={'color': 'blue', 'linestyle': 'solid', 'linewidth': 0.6, 'marker': '.', 'markersize': 3}
                                 )
    else:
        print("UNKNOWN Pushover Analysis type!!!");sys.exit()
    ops.wipe()
#=============================================================================
#    Plot
#=============================================================================
    if plot_Analysis == True:
        fp.plotPushoverX(outputDir) 
    if plot_StressStrain == True:
        if typeBuild == "CantileverColumn" or typeBuild == "buildBeam":
            fp.plotStressStrain(outputDirWalls,tagEleListToRecord_wall)
        elif typeBuild == 'coupledWalls':
            fp.plotStressStrain(outputDirWalls,tagEleListToRecord_wall)
            fp.plotStressStrain(outputDirBeams,tagEleListToRecord_beam)
    
    if plot_MomCurv == True:
        EIeff_walls = []
        for tagEle in tagEleListToRecord_wall:
            EIeff_walls.append(fp.plotMomCurv(outputDirWalls, tagEle, wall, typeBuild))

end_time        = time.time()
elapsed_time    = end_time - start_time
mins            = int(elapsed_time/60)
secs            = int(elapsed_time%60)
print(f"\nElapsed time: {mins} min + {secs} sec")
print("\n$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$")
print("The analysis was run successfully.")
print("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$")

if recordToLog == True:
    sys.stdout.close()
    sys.stdout = sys.__stdout__

winsound.Beep(440, 300)  # generate a 440Hz sound that lasts 300 milliseconds
winsound.Beep(440, 300)







