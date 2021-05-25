# -*- coding: utf-8 -*-
"""
Created on Fri Mar 15 13:06:51 2019

@author: fli034
"""
import os
import random
import time
import pygimli as pg
from pygimli.physics import TravelTimeManager
import matplotlib.pyplot as plt
import numpy as np
from pygimli.physics.traveltime import load as loadGimliTT

def resamplePicks(pickFile,fract2Rmv,topoFile,minError,maxError):
    f = open(pickFile,'r')
    #Split the lines based on the newline charachter (nre row in list ever \n)
    lines = f.read().split('\n')
    #Close the file for good file keeping
    f.close()
    
    #Get total number of picks
    totalPicks = len(lines)
    #Calculate the total number to keep (1 - fract2Remove)
    picks2Keep = int(np.round((1-fract2Rmv)*totalPicks,0))
    
    
    topo = np.loadtxt(topoFile)    
    #geoLocs is the topographically adjusted x-location
    geoLocs = topo[:,0]
    #oldLocs is the x-location of the geophnes. This needs to match the pick file
    oldLocs = topo[:,2]   
    
    #Get a list of indicies to keep using the random package
    lines = random.sample(lines,picks2Keep)
    gimliPicks = np.zeros((1,3))   
    #Begin to loop over all ines in the file
    for i in range(0,len(lines)):
        #extract the first line and store it as a string (check with print(tmpString))
        tmpString = lines[i]        
        #If the first charachter is not a space then proceed
        if tmpString != '':
            #Split the string into a list based on spaces. Modify here if it's comma separated
            data = tmpString.split()
            #xPick = float(data[1])/gx + 1
            #This uses logic to find the index of the value closes to the geophone location
            #The index is the line number in the topography file
            xPick = np.argmin(np.sqrt((float(data[1])-oldLocs)**2))+1
            #print(xPick)
            #Extract the picked travel time
            tPick = float(data[2])
            #sPick = float(data[0])/gx + 1
            #Extract the index of the shot location
            ##The index is the line number in the topography file
            sPick = np.argmin(np.sqrt((float(data[0])-oldLocs)**2))+1       
        #print([xPick,sPick])
        #Do some error checking. Pygimli does not like repeat picks or the pick at the shot location
        #If the shot location does not equal (!=) the geophone location and the pick is greater than zero
        #Then add the pick to the 3 column matrix that is going to be written out.
        if tPick > 0 and xPick != sPick:
            tmpMat = [sPick,xPick,tPick]
            gimliPicks = np.row_stack((gimliPicks,tmpMat))
            #print([sPick,xPick,tPick])       
        #Go to the next line in the file
        i = i + 1   
    #Remove the first row the matrix because I initalized it with zeros.
    gimliPicks = np.delete(gimliPicks,0,axis=0)
   
    #Calculate Error
    minShot = geoLocs[int(np.min(gimliPicks[:,0])-1)]
    maxShot = geoLocs[int(np.max(gimliPicks[:,0])-1)]
    lineLength = maxShot - minShot
    print('Line Length = ' + str(lineLength) + ' m')
    
    
    """
    ****************** GIMLI INPUT FILE DESCRIPTION *************************
    numberOfRows(integer) #Comment line
    x1 e1
    x2 e2
    .
    .
    .
    xn en
    numberOfPicks (integer) #Comment Line
    # s g t e[optional] <-- REQUIRED LINE NOT A COMMENT!!
    s1 g1 t1
    s2 g2 t2
    .
    .
    .
    sn gn tn
    
    **Note sn, gn, tn are index values to the row number in the topography above
    
    The gimli input file has two pieces to it. The first piece is the topography
    the first comment line tells you that its the geophone locations and elevaitons
    these values will be referenced as indexs. Thus index = 1 will be the first
    geophone/elevation pair listed
    
    After the topogrpahy is written a line telling gimli how many picks there
    are is required. A comment line with teh flags s, g, t, e[optiona] is required
    s = source index, g = geophone index, t = travel-time pick (seconds),
    e = error in seconts (picking error for chi^2 calc)
    
    The bigest difference between the input file and gimli input is that the 
    gimli input references everythign as indicies not values. The indicies are
    from the topgraphy at the top of the file.
    """
    #Write out inputfile
    #Rename the input file. Remove the .txt by spliting the file using a .
    #Then add _gimliInput.txt to the file
    #Possible add path here? Bugs will occcur if a file name has a . in it...
    gimliInFile = pickFile.split('.')[0] + '_gimliInput.txt' 
    #Open/create this file for writing
    f = open(gimliInFile,'w')    
    
    f.write(str(topo.shape[0])+ ' # geophone locations\n')
    for i in range(0,topo.shape[0]):
        f.write(str(topo[i,0])+' ' +str(topo[i,1]) + '\n')
        
    f.write(str(gimliPicks.shape[0]) + ' #measurements\n')
    f.write('#s g t err\n')     
    for i in range(0,gimliPicks.shape[0]):
        offset = np.abs(geoLocs[int(gimliPicks[i,0]-1)]-geoLocs[int(gimliPicks[i,1]-1)])
        err2write = maxError*offset/lineLength + minError
        str2Write = '{:4d} {:4d} {:0.5f} {:0.5f} \n'.format(int(gimliPicks[i,0]), int(gimliPicks[i,1]),gimliPicks[i,2],err2write)
        f.write(str2Write)
    f.close()
    
    return gimliInFile
    
    
    
  

def calcUnitVector(SoL,EoL):
    de = EoL[0] - SoL[0]
    dn = EoL[1] - SoL[1]
    mag = np.sqrt(de**2 + dn**2)
    ux = de/mag
    uy = dn/mag

    return ux,uy

def spatiallyLocateVTKFile(vtkFile,SoL,EoL):
    
    
    file1 = open(vtkFile, 'r') 
    outFileName = vtkFile.split('.')[0] + '_UTM.vtk' 
    file2 = open(outFileName,'w')
    count = 0
     
    ux,uy = calcUnitVector(SoL,EoL)
    
    while True: 
        count += 1
      
        # Get next line from file 
        line = file1.readline() 
        if count == 5:
            nPoints = line.split(' ')
            nPoints = int(nPoints[1]) + count
            print(nPoints)
        
        if count > 5 and count <= nPoints:
            #print(line)
            tempLine = line.split('\t')
            xDist =float(tempLine[0])
            elev = float(tempLine[1])
            easting = xDist*ux + SoL[0]
            northing = xDist*uy + SoL[1]
            line2Write = '{0:10.2f}\t{1:10.2f}\t{2:10.2f}\n'.format(easting,northing,elev)
            file2.write(line2Write)
            #print(line2Write)
        else:
            file2.write(line)
            
        # if line is empty 
        # end of file is reached 
        if not line: 
            break
      
    file1.close() 
    file2.close()

def plotAndSaveResults(path,folderName,mesh,ra,start,end,finalModel):
    
    #Plotting paramters. Here instead of in function call for now.
    minVelocity = 300
    maxVelocity = 4500
    minVertGrad = 0
    maxVertGrad = 300
    minHist = 0
    maxHist = 25
    minRayCoverage = 0
    maxRayCoverage = 10
    
    #Create the folder
    os.makedirs(path+folderName)
    
  
    
    #Calculate Vertical Gradient
    if finalModel:
        vels = np.loadtxt(path + 'allVelocityModels.txt')
        vel = np.mean(vels,axis=1)
        vertGrad = pg.solver.grad(mesh,vel)
        std = np.std(vels,axis=1)
        
        response = ra.fw.fop.response(1./vel)
        chi2 = ra.fw.chi2(response)
        modelRMS = pg.utils.rms(ra.fw.dataVals-response)*1000
        
        #Print out important model Inversion results
        s2 = "Chi2 = " + str("{:0.4f}".format(chi2)) + "\n"
        s3 = "RMS = " + str("{:0.4f}".format(modelRMS)) + " ms\n"
        print(s2)
        print(s3)
        outFile = open(path+folderName+'/'+'InversionStats.txt','w')
        outFile.write(s2+s3)
        outFile.close()  
    else:
        vel = ra.paraModel()
        vertGrad = pg.solver.grad(mesh,vel)
        std = np.zeros(len(vel))
        
        #Print out important model Inversion results
        s1 = "Time Elapse (min): " + str("{:0.2f}".format((end - start)/60)) + " min\n"
        s2 = "Chi2 = " + str("{:0.4f}".format(ra.inv.chi2())) + "\n"
        s3 = "RMS = " + str("{:0.4f}".format(ra.inv._inv.absrms())) + " s\n"
        s4 = "Total Iterations = " + str("{:0.0f}".format(ra.inv.maxIter)) + "\n"
        print(s1)
        print(s2)
        print(s3)
        print(s4)
        outFile = open(path+folderName+'/'+'InversionStats.txt','w')
        outFile.write(s1+s2+s3+s4)
        outFile.close()    
    
    
    #Save Results as a .vtk file (Paraview)
    mesh.addData('Velocity',vel)
    mesh.addData('Vertical Velocity Gradient',-vertGrad[:,1])
    mesh.addData('Ray Coverage',ra.coverage())
    mesh.addData('Standard Ray Coverage',ra.standardizedCoverage())
    mesh.exportVTK(path+folderName+'/'+pickFile.split('.')[0] + '_gimliResults.vtk' )
    
    
    #Spatially Locate the Profile
    spatiallyLocateVTKFile(path+folderName+'/'+pickFile.split('.')[0] + '_gimliResults.vtk',SoL,EoL)
    
    
    
    
    
    #******************** PLOT AND SAVE RESULTS ******************************
    fig1 = plt.figure('Model Results')
    fig1.set_size_inches([8,8])
    ax = fig1.add_subplot(211)
    pg.show(mesh, vel, coverage=ra.standardizedCoverage(),cMin=minVelocity,cMax=maxVelocity,cMap='jet_r',ax=ax)
    pg.viewer.mpl.drawSensors(ax, data.sensorPositions(), diam=0.5, color="k")
    ax.set_xlabel('Distance (m)')
    ax.set_ylabel('Elevation (m)')
    ax.set_title('Velocity (m/s)')
    ax2 = fig1.add_subplot(212)
    pg.show(mesh, data=-vertGrad[:,1], coverage=ra.standardizedCoverage(),cMin=minVertGrad,cMax=maxVertGrad,ax=ax2,cMap='plasma')
    pg.viewer.mpl.drawSensors(ax2, data.sensorPositions(), diam=0.5, color="k")
    ax2.set_xlabel('Distance (m)')
    ax2.set_ylabel('Elevation (m)')
    ax2.set_title('Vertical Gradient (m/s/m)')
    fig1.savefig(path+folderName+'/'+'ModelResults.png',dpi=600)
    #ra.showResult(cMin=500,cMax=4500,cmap='jet_r')
    
        
    fig2 = plt.figure('Ray Coverage')
    fig2.set_size_inches([8,8])
    ax = fig2.add_subplot(211)
    pg.show(ra.paraDomain, vel, coverage=ra.standardizedCoverage(),cMin=minVelocity,cMax=maxVelocity,cMap='jet_r',ax=ax)
    pg.viewer.mpl.drawSensors(ax, data.sensorPositions(), diam=0.5, color="k")
    ra.drawRayPaths(ax,color="k")
    ax.set_xlabel('Distance (m)')
    ax.set_ylabel('Elevation (m)')
    ax.set_title('Ray Paths')
    ax2 = fig2.add_subplot(212)
    tmpCov = ra.coverage()
    tmpCov[np.isneginf(tmpCov)] = 0
    pg.show(mesh, data=tmpCov, coverage=ra.standardizedCoverage(),cMin=minRayCoverage,cMax=maxRayCoverage,ax=ax2,cMap='plasma')
    pg.viewer.mpl.drawSensors(ax2, data.sensorPositions(), diam=0.5, color="k")
    ax2.set_xlabel('Distance (m)')
    ax2.set_ylabel('Elevation (m)')
    ax2.set_title('Ray Coverage (# of rays in each cell)')
    fig2.savefig(path+folderName+'/'+'ModelRayCoverage.png',dpi=600)
    
    
    fig3 = plt.figure('Model Fits')
    fig3.set_size_inches([11,5])
    ax1 = fig3.add_subplot(121)
    pickTT = np.array(ra.data('t'))
    modTT = np.array(ra.inv.response)
    ax1.hist2d(pickTT*1000,(pickTT-modTT)*1000,bins=75,vmin=minHist,vmax=maxHist,cmap='plasma')
    ax1.plot([0,200],[0,0],'k--',linewidth=2)
    ax1.set_xlabel('Observed Travel-times (ms)')
    ax1.set_ylabel('Residual (Obs - Mod) (ms)')
    
    ax2 = fig3.add_subplot(122)
    ax2.hist2d(pickTT*1000,modTT*1000,bins=100,vmin=minHist,vmax=maxHist,cmap='plasma')
    ax2.plot([0,200],[0,200],'k',linewidth=2)
    ax2.plot([0,200],[-5,195],'k--',linewidth=2)
    ax2.plot([0,200],[5,205],'k--',linewidth=2)
    ax2.set_xlabel('Observed Travel-times (ms)')
    ax2.set_ylabel('Modelted Travel-times (ms)')
    fig3.savefig(path+folderName+'/'+'ModelFits.png',dpi=600)
    
    shotLocs = np.array(ra.data('s'))
    pickTT = np.array(ra.data('t'))
    modTT = np.array(ra.inv.response)
    uniqShots_pwave = np.unique(shotLocs) 
    rmsUniqShots_pwave = np.zeros(len(uniqShots_pwave))
    shotNum = 0
    for shot in uniqShots_pwave:
        ind = np.where(shot==shotLocs)[0]
        rmsUniqShots_pwave[shotNum] = np.sqrt(np.sum((pickTT[ind]-modTT[ind])**2)/len(ind))*1000
        shotNum = shotNum + 1
        
    fig4 = plt.figure('Missfit by Shot Location')
    fig4.set_size_inches([11.63,  4.52])
    
    ax2 = fig4.add_subplot(111)
    ax2.bar(uniqShots_pwave,rmsUniqShots_pwave,width=5,color='tab:blue')
    ax2.set_xlabel('Shot Location')
    ax2.set_ylabel('RMS Missfit (ms)')
    ax2.set_title('P-wave Inversion')
    fig4.tight_layout()
    fig4.savefig(path+folderName+'/'+'MissfitByShot.png',dpi=600)
    
    modTT = np.array(ra.inv.response)
    sensorPostions_x  =np.array(ra.data.sensorPositions())[:,0]
    #test = np.column_stack((ra.data('s'),ra.data('g'),modTT,(pickTT-modTT)*1000))
    test = np.column_stack((ra.data('s'),ra.data('g'),modTT))
    for i in range(0,len(test)):
        test[i,0] = sensorPostions_x[int(test[i,0])]
        test[i,1] = sensorPostions_x[int(test[i,1])]
    np.savetxt(path+folderName+'/'+pickFile.split('.')[0]+'_Pwave_modeledTravelTimes.txt',test,fmt='%10.4f')
    
    
    
if __name__ == "__main__":
    
    #************************** INPUTS ***************************************
    
    """
    Space seperated file with shot location as the first column, geophone 
    location as the second column, and travel-time (in seconds) in the third
    column 
    """
    pickFile = 'TLV_2016_0816_picks.txt' #Pick file from picker.py
    
    """
    The topo file is a three column file sepearted by spaces. The first column
    is the x-locaiton that has been adjusted to map space. This number will be
    a few centimeters smaller than the field value because it is the land
    distance projected into map distance. The second column is the elevation
    in meters and the third column is the geophone location at its orignal
    spacinge. I call it it the old locaitons. This distance will be too long
    on lines with a lot of topography. The transformation is important for 
    projecting the VTK files to the topography. The old locs column will 
    match the locations in the pick file!!
    """
    topoFile = 'TLV_2016_0816_topoLiDARm.txt'
    
    
    #Line Location for spatially located VTK and GMT Files
    #MUST MATCH THE LOCATION USED TO EXTRACT THE TOPOGRAPHY!!!
    SoL = [464085.77,4563854.83]
    EoL = [464000.61,4563481.57]
    
    """    
    Set the errors for the inversion weighting. This uses a simple linear
    function to generate errors in seconds:
        error = offset*maxError/lineLength + minError
    It is advised that you don't go less than 0.001 (1 ms) for the minErro.
    I usally fit about 5-6 ms on the farthest offset picks. Thus default values
    are maxError = 0.005 and minError = 0.001
    """
    maxError = 0.005
    minError = 0.001
    
    """
    MESH PARAMTERS
    maxDepth = maximum depth of the model sets (paraDepth in function)
    meshDx = 0.33; a decimal number defining how many nodes between geophones. 
        0.33 = 3. I you wanted 10 it would be 0.1. If you wanted two it would 
        be 0.5 (sets paraDX in function)
    maxTriArea = maximum area of each triangle. The smaller you make this
        number the more cells you will have. This means the more accurate your
        forward model is but at the cost of cmputaiton time. (sets 
        paraMaxCellSize in functino call)
    """
    maxDepth = 60
    meshDx = 0.3
    maxTriArea = 7
    
    """
    INVERSION PARAMTERS
    smoothing = sets lam in the model. Make this number larger for smoothing. 
        This is the scalar value on teh weighting matrix
    vertWeight = sets the zWeight in the function. This value controls how
        much vertical versus horizontal smoothing there is. This value needs
        to stay between 0 and 1. A value of 1 is a very laterally smooth model
        a value close to zero lets a lot of lateral hetergeneity in.
    minVel = this is the minimum velocity at the top of the model
    maxVel = this is the maximum velocity at teh bottom of the model    
    """
    smoothing = 100
    vertWeight = 0.1
    minVel = 1000
    maxVel = 4500
    maxIterations = 10
    
    """
    Bootstrapping parameters. This algorithm will run the inversion over and 
    over again witht the same inversion parameters except that it removes
    a percentage (or fraction) of the data prior to each run. This works well
    if you have a couple thousand picks.
    fract2Remove = the fraction of picks to remove each bootstrap iteration
    nRuns = total number of bootstrap iterations
    path2Runs = path to the runs folder it must have the / at the end since it's a directory
    """
    fract2Remove = 0.4
    nRuns = 2
    path2Runs = "/Users/bflinch/Dropbox/Clemson/Research/GithubDevel/PyGimliRefractionInversionScripts/BootstrapInversion/runs/"
    #*************************************************************************
    
    
    #********** DO INVERSION***************************************************
    if nRuns > 1:
        for i in range(0,nRuns):
            folderName = "run"+str(i+1)
            #Convert the file
            inFile = resamplePicks(pickFile,fract2Remove,topoFile, minError, maxError)
            #Read in the file with gimli method
            data = loadGimliTT(inFile)
            
            #Create the structure required for inversion using TravelTimeManager and data container
            ra = TravelTimeManager(data)
            #Generate the mesh
            mesh = ra.createMesh(data=ra.data,paraDepth=maxDepth,paraDX=meshDx,paraMaxCellSize=maxTriArea)
                    
            
            #Do the Inversion (time the inversion)
            start = time.time()
            ra.inv.maxIter = maxIterations
            ra.invert(data=ra.data,mesh=mesh,lam=smoothing,zWeight=vertWeight,useGradient=True,vTop=minVel,vBottom=maxVel)
            end = time.time()
            
            if i == 0:
                velModels = np.zeros((len(mesh.cellCenters()),nRuns))
                velModels[:,i] = ra.paraModel()
            else:
                velModels[:,i] = ra.paraModel()
                np.savetxt(path2Runs+'allVelocityModels.txt',velModels) #fail safe write the data each iteration in case of crash
            
            print('\n\nCompleted Run ' + str(i+1) + ' of ' + str(nRuns) + ' Inversion...\n\n')
            plotAndSaveResults(path2Runs,folderName,mesh,ra,start,end,finalModel=False)
            data.save(path2Runs+folderName+'/'+pickFile.split('.')[0]+'_run'+str(i+1)+'.txt')
        
        #AVERAGE THE VELOCITIES AND RERUN THE MODEL
        print('\n\n******Starting Final Averaged Run*****\n\n')
        inFile = resamplePicks(pickFile,0,topoFile,minError,maxError)
        #Read in the file with gimli method
        data = loadGimliTT(inFile)
        #Create the structure required for inversion using TravelTimeManager and data container
        ra = TravelTimeManager(data)
        #Generate the mesh
        mesh = ra.createMesh(data=ra.data,paraDepth=maxDepth,paraDX=meshDx,paraMaxCellSize=maxTriArea)
        avgVel = np.mean(velModels,axis=1)
        start = time.time()
        ra.inv.maxIter = 0
        ra.invert(data=ra.data,mesh=mesh,lam=smoothing,zWeight=vertWeight,startModel=1./avgVel,verbose=True)
        end = time.time()   
        

        print('\n\n******Completed Final Averaged Run*****\n\n')
        folderName = 'FinalAveragedRun'
        finalModel = True
        plotAndSaveResults(path2Runs,folderName,mesh,ra,start,end,finalModel)
        np.savetxt(path2Runs+folderName+'/'+'allVelocityModels.txt',velModels)
    else:
            #Convert the file
            folderName = "singleInversionResults"
            inFile = resamplePicks(pickFile,0,topoFile, minError, maxError)
            #Read in the file with gimli method
            data = loadGimliTT(inFile)
            
            #Create the structure required for inversion using TravelTimeManager and data container
            ra = TravelTimeManager(data)
            
            #Generate the mesh
            mesh = ra.createMesh(data=ra.data,paraDepth=maxDepth,paraDX=meshDx,paraMaxCellSize=maxTriArea)
            #Plot the mesh
            #pg.show(mesh=mesh)
            
            #Do the Inversion (time the inversion)
            start = time.time()
            ra.inv.maxIter = maxIterations
            ra.invert(data=ra.data,mesh=mesh,lam=smoothing,zWeight=vertWeight,useGradient=True,vTop=minVel,vBottom=maxVel)
            end = time.time()
            
            print("\n\nCompleted Single Inversion...\n\n")
            plotAndSaveResults(path2Runs,folderName,mesh,ra,start,end)
















