import sys
import openseespy.opensees     as ops
import opsvis                  as opv

H_story_List    = [6, *(3*[6.])]
L_Bay_List      = [7., 7., 5.]
Lw              = 4.
PHL             = 2. 
numSegBeam      = 6
numSegWall      = 3


def coupledWalls(H_story_List, L_Bay_List, Lw, tagSec, numSegBeam, numSegWall, PHL):
    
    
    for L_Bay in L_Bay_List:
        if L_Bay <= Lw:
            print(f"L_Bay={L_Bay} <= Lw={Lw}")
            print('The program will exit now!'); sys.exit()
    
    # gridList        = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z']
    H_story_List    += [0]
    L_Bay_List      += [0]
    
    
    
    coords          = {}
    y               = 0.
    for storyNum, storyH in enumerate(H_story_List):
        # storyKey        = f"Story_{storyNum:02}"
        # print(storyKey)
        # print(f"Story Elevation = {y}")
        
        x           = 0.
        for gridIndex, L_Bay in enumerate(L_Bay_List):
            # gridKey         = f"{gridList[gridIndex]}"
            # print(gridKey)
            # print(f"Grid x = {x}")
            tagNode         = int(f"1{storyNum:02}{gridIndex:02}0")
            coords[tagNode] = [x, y]
            
            if gridIndex == 0 and y != 0:
                tagNode     = int(f"1{storyNum:02}{gridIndex:02}2")
                coords[tagNode] = [x+Lw/2, y]
            elif gridIndex > 0 and gridIndex < len(L_Bay_List)-1 and y != 0:
                tagNode     = int(f"1{storyNum:02}{gridIndex:02}1")
                coords[tagNode] = [x-Lw/2, y]
                tagNode     = int(f"1{storyNum:02}{gridIndex:02}2")
                coords[tagNode] = [x+Lw/2, y]
                
            x               += L_Bay
        y               += storyH
            
    #   Build Model
    ops.wipe()
    ops.model('basic', '-ndm', 2, '-ndf', 3)
    
    #   Define Nodes
    for tagNode, coord in coords.items():
        # tagCoordY   = f"{tagNode}"[1:-3]
        # tagSuffix   = f"{tagNode}"[-1]
        ops.node(tagNode, *coord)
        # if tagCoordY == '00' and tagSuffix != '0': #this is to skip the nodes at base whose suffices are not 0
        #     continue
        # else:
        #     # print(f"tagNode = {tagNode}\tcoord = {coord}")
        #     ops.node(tagNode, *coord)
    
    #   Put base node tags into a list
    tagNodeBaseList = []
    for tagNode, coord in coords.items():
        tagCoordY   = f"{tagNode}"[1:-3]
        if tagCoordY == '00':
            tagNodeBaseList.append(tagNode)
    
    #   Assign Constraints
    ops.fixY(0, *[1, 1, 1], '-tol', 1e-3)
    
    
    #   Define Geometric Transformation
    tagGTLinear = 1
    tagGTPDelta = 2
    ops.geomTransf('Linear', tagGTLinear)
    ops.geomTransf('PDelta', tagGTPDelta)
    
    #       Define beamIntegrator
    # tagInt      = 1
    # NIP         = 5
    # ops.beamIntegration('Legendre', tagInt, tagSec, NIP)  # 'Lobatto', 'Legendre' for the latter NIP should be odd integer.
    
    
    #   Define material and sections
    A, E, I = 1e5, 200e9, 1
    ops.uniaxialMaterial('Elastic', 1, E)
    
    #######################################################################################################
    # Define Element
    #######################################################################################################
    #   Walls:
    ##  Define tags of Walls and LeaningColumns
    
    def discretizeWall(tagNodeI, tagNodeJ, tagCoordXI, tagCoordYI, tagCoordYJ, Walls, coordsGlobal, PHL, numSegWall=1):
        
        xI  = coordsGlobal[tagNodeI][0];    yI  = coordsGlobal[tagNodeI][1]
        xJ  = coordsGlobal[tagNodeJ][0];    yJ  = coordsGlobal[tagNodeJ][1]
        
        Lx  = xJ - xI; Ly = yJ - yI
        L   = (Lx**2 + Ly**2)**0.5
        PHR = PHL/L
        lx  = PHR*Lx/numSegWall; ly = PHR*Ly/numSegWall
        
        coordsLocal = {}
        for i in range(0, numSegWall+1):
            tagNode = tagNodeI + i
            coordsLocal[tagNode] = [xI + i*lx, yI + i*ly]
            if i > 0:
                ops.node(tagNode, *coordsLocal[tagNode])
                tagElement = int(f"5{tagCoordXI}{tagCoordYI}{tagCoordYJ}{i}")
                Walls[tagElement]  = [tagNode-1, tagNode ]
                print(f"Wall{tagElement} = {Walls[tagElement]}")
                print(f"NodeI({tagNode-1}) = {coordsLocal[tagNode-1]}")
                print(f"NodeJ({tagNode}) = {coordsLocal[tagNode]}")
        tagElement = int(f"5{tagCoordXI}{tagCoordYI}{tagCoordYJ}{0}")
        Walls[tagElement] = [tagNode,   tagNodeJ]
        print(f"Wall{tagElement} = {Walls[tagElement]}")
        print(f"NodeI({tagNode}) = {coordsLocal[tagNode]}")
        print(f"NodeJ({tagNodeJ}) = {coordsGlobal[tagNodeJ]}")
        print("End")
        # return(0)
        
    gridLeaningColumn = f"{(len(L_Bay_List)-1):02}"
    
    Walls           = {}
    LeaningColumns  = {}
    for tagNode, coord in coords.items():
        print("LOOP1:")
        print(f"tagNode = {tagNode}\tcoord = {coord}")
        if f"{tagNode}"[-1] == '0':
            tagNodeI    = tagNode
            tagCoordXI  = f"{tagNodeI}"[3:-1]
            tagCoordYI  = f"{tagNodeI}"[1:-3]
            print(f"tagCoordXI = {tagCoordXI}\ttagCoordYI = {tagCoordYI}")
        for tagNode, coord in coords.items():
            print("LOOP2:")
            print(f"tagNode = {tagNode}\tcoord = {coord}")
            if f"{tagNode}"[-1] == '0':
                tagNodeJ    = tagNode
                tagCoordXJ  = f"{tagNodeJ}"[3:-1]
                tagCoordYJ  = f"{tagNodeJ}"[1:-3]
                print(f"tagCoordXJ = {tagCoordXJ}\ttagCoordYJ = {tagCoordYJ}")
                
                if tagCoordXI == tagCoordXJ: # this makes it a column
                    print("tagCoordXI == tagCoordXJ")
                    if int(tagCoordYJ) - int(tagCoordYI) == 1: 
                        print("int(tagCoordYJ) - int(tagCoordYI) == 1")
                        print(f"tagCoordXI={tagCoordXI}    gridLeaningColumn={gridLeaningColumn}")
                        if tagCoordXI != gridLeaningColumn:
                            print("tagCoordXI != gridLeaningColumn")
                            print(f"{tagNodeI} VS {tagNodeJ} ==> tagWall = 5{tagCoordXI}{tagCoordYI}{tagCoordYJ}")
                            if int(tagCoordYJ) == 1:
                                print("int(tagCoordYJ) == 1")
                                discretizeWall(tagNodeI, tagNodeJ, tagCoordXI, tagCoordYI, tagCoordYJ, Walls, coords, PHL, numSegWall)
                                # Walls[f"5{tagCoordXI}{tagCoordYI}{tagCoordYJ}"] = [tagNodeI, tagNodeJ]  #Prefix 5 is for Walls
                            else:
                                print("int(tagCoordYJ) != 1")
                                tagElement = int(f"5{tagCoordXI}{tagCoordYI}{tagCoordYJ}{0}")
                                Walls[tagElement] = [tagNodeI, tagNodeJ]  #Prefix 5 is for Walls
                        else:
                            tagElement = int(f"2{tagCoordXI}{tagCoordYI}{tagCoordYJ}")
                            LeaningColumns[tagElement] = [tagNodeI, tagNodeJ]  #Prefix 2 is for LeaningColumns
    
    ##  Define Walls
    for tagElement, tagNodes in Walls.items():
        # print(f"tagElement = {tagElement} & tanNodes = {tagNodes}")
        ops.element('elasticBeamColumn', tagElement, *tagNodes, A, E, I, tagGTPDelta)
        
    ##  Define LeaningColumns
    for tagElement, tagNodes in LeaningColumns.items():
        # print(f"tagElement = {tagElement} & tanNodes = {tagNodes}")
        ops.element('Truss', tagElement, *tagNodes, A, 1)
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    #   Rigid Beams:
    ##  Define tags of Rigid Beams
    RBeams = {}
    for tagNode, coord in coords.items():
        tagNodeI    = tagNode
        tagCoordXI  = f"{tagNodeI}"[3:-1]
        tagCoordYI  = f"{tagNodeI}"[1:-3]
        tagSuffixI  = f"{tagNodeI}"[-1]
        for tagNode, coord in coords.items():
            tagNodeJ    = tagNode
            tagCoordXJ  = f"{tagNodeJ}"[3:-1]
            tagCoordYJ  = f"{tagNodeJ}"[1:-3]
            tagSuffixJ  = f"{tagNodeJ}"[-1]
            
            if tagCoordXI == tagCoordXJ and tagCoordYI == tagCoordYJ and (tagCoordYI != '00' or tagCoordYJ != '00'):
                if tagSuffixJ == '0' or tagSuffixI == '0':
                    if int(tagSuffixJ)-int(tagSuffixI) == -1 or int(tagSuffixJ)-int(tagSuffixI) == 2:
                        # print(f"{tagNodeI} VS {tagNodeJ} ==> tagRBeam = 5{tagCoordYI}{tagCoordXI}{tagSuffixI}{tagSuffixJ}")
                        tagElement = int(f"6{tagCoordYI}{tagCoordXI}{tagSuffixI}{tagSuffixJ}")
                        RBeams[tagElement] = [tagNodeI, tagNodeJ]  #Prefix 6 is for RBeams
                    
    ##  Define Rigid Beams
    for tagElement, tagNodes in RBeams.items():
        # print(f"tagElement = {tagElement} & tanNodes = {tagNodes}")
        ops.element('elasticBeamColumn', tagElement, *tagNodes, 1e10*A, 1e10*E, 1e10*I, tagGTLinear)
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    #   Beams and Trusses:
    ##  Define tags of  Beams and Trusses
    def discretizeBeam(tagNodeI, tagNodeJ, tagCoordYI, tagCoordXI, tagCoordXJ, Beams, coordsGlobal, numSegBeam=1):
        
        xI  = coordsGlobal[tagNodeI][0];    yI  = coordsGlobal[tagNodeI][1]
        xJ  = coordsGlobal[tagNodeJ][0];    yJ  = coordsGlobal[tagNodeJ][1]
        
        Lx = xJ - xI; Ly = yJ - yI
        
        lx = Lx/numSegBeam; ly = Ly/numSegBeam
        
        coordsLocal = {}
        for i in range(0, numSegBeam):
            tagNode = tagNodeI + i
            coordsLocal[tagNode] = [xI + i*lx, yI + i*ly]
            if i > 0:
                ops.node(tagNode, *coordsLocal[tagNode])
                tagElement = int(f"4{tagCoordYI}{tagCoordXI}{tagCoordXJ}{i}")
                Beams[tagElement]  = [tagNode-1, tagNode ]
        tagElement = int(f"4{tagCoordYI}{tagCoordXI}{tagCoordXJ}{numSegBeam}")
        Beams[tagElement] = [tagNode,   tagNodeJ]
        
        # return(0)
        
    Beams   = {}
    Trusses = {}
    for tagNode, coord in coords.items():
        tagNodeI    = tagNode
        tagCoordXI  = f"{tagNodeI}"[3:-1]
        tagCoordYI  = f"{tagNodeI}"[1:-3]
        tagSuffixI  = f"{tagNodeI}"[-1]
        for tagNode, coord in coords.items():
            tagNodeJ    = tagNode
            tagCoordXJ  = f"{tagNodeJ}"[3:-1]
            tagCoordYJ  = f"{tagNodeJ}"[1:-3]
            tagSuffixJ  = f"{tagNodeJ}"[-1]
            
            if tagCoordXI != tagCoordXJ and tagCoordYI == tagCoordYJ and (tagCoordYI != '00' or tagCoordYJ != '00'):
                # build beam
                if tagCoordXJ != gridLeaningColumn:
                    if int(tagSuffixJ)-int(tagSuffixI) == -1 and int(tagCoordXJ)-int(tagCoordXI) == 1:
                        if tagSuffixI != '0' and tagSuffixJ != '0':
                            # print(f"{tagNodeI} VS {tagNodeJ} ==> tagBeam = 4{tagCoordYI}{tagCoordXI}{tagCoordXJ}")
                            discretizeBeam(tagNodeI, tagNodeJ, tagCoordYI, tagCoordXI, tagCoordXJ, Beams, coords, numSegBeam)
                            # Beams[f"4{tagCoordYI}{tagCoordXI}{tagCoordXJ}"] = [tagNodeI, tagNodeJ]  #Prefix 4 is for Beams
                # build truss
                elif tagCoordXJ == gridLeaningColumn and tagCoordXI == f"{(len(L_Bay_List)-2):02}":
                    if int(tagSuffixJ)-int(tagSuffixI) == -2:
                        # print(f"{tagNodeI} VS {tagNodeJ} ==> tagTruss = 2{tagCoordYI}{tagCoordXI}{tagCoordXJ}")
                        tagElement = int(f"3{tagCoordYI}{tagCoordXI}{tagCoordXJ}")
                        Trusses[tagElement] = [tagNodeI, tagNodeJ]  #Prefix 3 is for Trusses
    
    ##  Define Beams
    for tagElement, tagNodes in Beams.items():
        # print(f"tagElement = {tagElement} & tanNodes = {tagNodes}")
        ops.element('elasticBeamColumn', tagElement, *tagNodes, A, E, I, tagGTLinear)
        # ops.element('dispBeamColumn',    tagElement, *tagNodes, tagGTLinear, tagInt)
    ##  Define Trusses
    for tagElement, tagNodes in Trusses.items():
        # print(f"tagElement = {tagElement} & tanNodes = {tagNodes}")
        ops.element('Truss', tagElement, *tagNodes, A, 1)
        
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #   Define Top-Left corner node as Control Node
    for tagNode, coord in coords.items():
        tagCoordXI  = f"{tagNode}"[3:-1]
        tagCoordYI  = f"{tagNode}"[1:-3]
        tagSuffixI  = f"{tagNode}"[-1]
        if tagSuffixI == '0' and tagCoordXI == '00' and tagCoordYI == f"{storyNum:02}":
            tagNodeControl = tagNode
            print(f"tagNodeControl = {tagNodeControl}")

    return(tagNodeControl, tagNodeBaseList, x, y, coords)

tagSec = 1
tagNodeControl, tagNodeBaseList, buildingWidth, buildingHeight, coords = coupledWalls(H_story_List, L_Bay_List, Lw, tagSec, numSegBeam, numSegWall, PHL)

opv.plot_model(node_labels=0, element_labels=1, 
               fig_wi_he=(buildingWidth+10., buildingHeight+7.),
               fmt_model={'color': 'blue', 'linestyle': 'solid', 'linewidth': 1.2, 'marker': '.', 'markersize': 10})





















