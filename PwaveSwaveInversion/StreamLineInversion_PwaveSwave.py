# -*- coding: utf-8 -*-
"""
Created on Fri Mar 15 13:06:51 2019

@author: fli034
"""
import time
import pygimli as pg
from pygimli.physics import TravelTimeManager
import pygimli.meshtools as mt
import matplotlib.pyplot as plt
import numpy as np
from pygimli.physics.traveltime import load as loadGimliTT

def convertPicks2Gimli(pickFile,topoFile, minError, maxError):
    
    #Read in the topography File
    topo = np.loadtxt(topoFile)    
    #geoLocs is the topographically adjusted x-location
    geoLocs = topo[:,0]
    #oldLocs is the x-location of the geophnes. This needs to match the pick file
    oldLocs = topo[:,2]   
    #Rread in the file -- oldschool into a list.
    f = open(pickFile,'r')
    #Split the lines based on the newline charachter (nre row in list ever \n)
    lines = f.read().split('\n')
    #Close the file for good file keeping
    f.close()
    #Initialize gimliInput
    gimliPicks = np.zeros((1,3))   
    #Begin to loop over all ines in the file
    for i in range(0,len(lines)):
        #extract the first line and store it as a string (check with print(tmpString))
        tmpString = tmpString = lines[i]        
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
    
if __name__ == "__main__":
    
    #************************** INPUTS ***************************************
    
    """
    Space seperated file with shot location as the first column, geophone 
    location as the second column, and travel-time (in seconds) in the third
    column 
    """
    pickFile = 'CEF_2021_0315_Pwave_samePicks_after.txt' #Pick file from picker.py
    pickFile2 = 'CEF_2021_0315_Swave_samePicks_after.txt'
    
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
    topoFile = 'CEF_2021_0315_LiDAR.txt'
    
    
    #Line Location for spatially located VTK and GMT Files
    #MUST MATCH THE LOCATION USED TO EXTRACT THE TOPOGRAPHY!!!
    SoL = [327658.40,3847006.43]
    EoL = [327839.56,3847031.66]
    
    """    
    Set the errors for the inversion weighting. This uses a simple linear
    function to generate errors in seconds:
        error = offset*maxError/lineLength + minError
    It is advised that you don't go less than 0.001 (1 ms) for the minErro.
    I usally fit about 5-6 ms on the farthest offset picks. Thus default values
    are maxError = 0.005 and minError = 0.001
    """
    maxError = 0.005
    minError = 0.0015
    
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
    maxDepth = 70
    meshDx = 0.2
    maxTriArea = 5
    
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
    smoothing = 80
    vertWeight = 0.1
    minVel = 300
    maxVel = 4500
    maxIterations = 15
    #*************************************************************************
    
    
    #********** DO INVERSION***************************************************
    #Convert the file
    inFile = convertPicks2Gimli(pickFile,topoFile, minError, maxError)
    inFile2 = convertPicks2Gimli(pickFile2,topoFile, minError, maxError)
    #Read in the file with gimli method
    data = loadGimliTT(inFile)
    data2 = loadGimliTT(inFile2)
    
    #Create the structure required for inversion using TravelTimeManager and data container
    ra = TravelTimeManager(data)
    ra2 = TravelTimeManager(data2)
    

    #Generate the mesh
    mesh = ra.createMesh(data=ra.data,paraDepth=maxDepth,paraDX=meshDx,paraMaxCellSize=maxTriArea)
    #Plot the mesh
    #pg.show(mesh=mesh)
    
    #Do the Inversion (time the inversion)
    start1 = time.time()
    ra.inv.maxIter = maxIterations
    #ra.invert(data=ra.data,mesh=mesh,lam=smoothing,zWeight=vertWeight,useGradient=True,vTop=minVel,vBottom=maxVel,maxIter=maxIterations,limits=[100,5000])
    #ra2.invert(data=ra2.data,mesh=mesh,lam=smoothing,zWeight=vertWeight,useGradient=True,vTop=minVel-minVel/1.6,vBottom=maxVel-maxVel/1.6,maxIter=maxIterations,limits=[25,3000])
    ra.invert(data=ra.data,mesh=mesh,lam=smoothing,zWeight=vertWeight,useGradient=True,vTop=minVel,vBottom=maxVel,verbose=True)
    end1 = time.time()
    
    vel = ra.paraModel()
    vel2Start = vel/2
    ra2.inv.startModel = vel2Start
    ra2.inv.maxIter = maxIterations
    start2 = time.time()
    ra2.invert(data=ra2.data,mesh=mesh,lam=smoothing,zWeight=vertWeight,verbose=True,startModel=1/vel2Start)
    #ra2.invert(data=ra2.data,mesh=mesh,lam=smoothing,zWeight=vertWeight,useGradient=True,vTop=minVel/2,vBottom=maxVel/2,verbose=True)
    end2 = time.time()
    
    
    #Print out important model Inversion results
    print('***********P-WAVE RESULTS***********')
    s1 = "Time Elapse (min): " + str("{:0.2f}".format((end1 - start1)/60)) + " min\n"
    s2 = "Chi2 = " + str("{:0.4f}".format(ra.inv.chi2())) + "\n"
    s3 = "RMS = " + str("{:0.4f}".format(1000*ra.inv._inv.absrms())) + " ms\n"
    s4 = "Total Iterations = " + str("{:0.0f}".format(ra.inv.maxIter)) + "\n"
    print(s1)
    print(s2)
    print(s3)
    print(s4)
    outFile = open('InversionStats.txt','w')
    outFile.write(s1+s2+s3+s4)
    outFile.close()
    print('\n***********S-WAVE RESULTS***********')
    s1 = "Time Elapse (min): " + str("{:0.2f}".format((end2 - start2)/60)) + " min\n"
    s2 = "Chi2 = " + str("{:0.4f}".format(ra2.inv.chi2())) + "\n"
    s3 = "RMS = " + str("{:0.4f}".format(1000*ra2.inv._inv.absrms())) + " ms\n"
    s4 = "Total Iterations = " + str("{:0.0f}".format(ra2.inv.maxIter)) + "\n"
    print(s1)
    print(s2)
    print(s3)
    print(s4)
    outFile = open('InversionStats2.txt','w')
    outFile.write(s1+s2+s3+s4)
    outFile.close()
    #*************************************************************************
    
    #Calculate Vertical Gradient
    vel = ra.paraModel()
    sVel = ra2.paraModel()
    
    #vertGrad = pg.solver.grad(mesh,vel)
    vpvs = vel/sVel
    pr = (vel**2 - 2*sVel**2)/(2*(vel**2-sVel**2))
    #pr = np.array(pr)
    #ind = np.where(pr<0)
    #pr[ind] = 0
    #ind = np.where(pr>0.5)
    #pr[ind] = 0.5
    
    #Save Results as a .vtk file (Paraview)
    mesh.addData('Pwave Velocity',vel)
    mesh.addData('Swave Velocity',sVel)
    mesh.addData('Vp/Vs',vpvs)
    mesh.addData('cov',ra.coverage())
    mesh.addData('stdCov',ra.standardizedCoverage())
    mesh.exportVTK(pickFile.split('.')[0] + '_gimliResults.vtk' )
    
    
    #Spatially Locate the Profile
    spatiallyLocateVTKFile(pickFile.split('.')[0] + '_gimliResults.vtk',SoL,EoL)
    
    
    
    
    
    #******************** PLOT AND SAVE RESULTS ******************************
    fig1 = plt.figure('Model Results')
    fig1.set_size_inches([12,8])
    ax = fig1.add_subplot(321)
    vel = ra.paraModel()
    pg.show(mesh, vel, coverage=ra.standardizedCoverage(),cMin=100,cMax=4000,cmap='jet_r',ax=ax)
    pg.viewer.mpl.drawSensors(ax, data.sensorPositions(), diam=0.5, color="k")
    ax.set_xlabel('Distance (m)')
    ax.set_ylabel('Elevation (m)')
    ax.set_title('P-wave Velocity (m/s)')
    
    ax = fig1.add_subplot(322)
    vel = ra.paraModel()
    vertGrad = pg.solver.grad(mesh,vel)
    pg.show(mesh, -vertGrad[:,1], coverage=ra.standardizedCoverage(),cMin=50,cMax=400,cmap='viridis',ax=ax)
    pg.viewer.mpl.drawSensors(ax, data.sensorPositions(), diam=0.5, color="k")
    ax.set_xlabel('Distance (m)')
    ax.set_ylabel('Elevation (m)')
    ax.set_title('P-wave Vertical Velocity Gradient (m/s/m)')
    
    
    ax2 = fig1.add_subplot(323)
    pg.show(mesh, data=sVel, coverage=ra2.standardizedCoverage(),cMin=100,cMax=2500,cmap='jet_r',ax=ax2)
    pg.viewer.mpl.drawSensors(ax2, data.sensorPositions(), diam=0.5, color="k")
    ax2.set_xlabel('Distance (m)')
    ax2.set_ylabel('Elevation (m)')
    ax2.set_title('S-wave Velocity (m/s)')
    
    ax = fig1.add_subplot(324)
    vel = ra2.paraModel()
    vertGrad = pg.solver.grad(mesh,vel)
    pg.show(mesh, -vertGrad[:,1], coverage=ra2.standardizedCoverage(),cMin=50,cMax=200,cmap='viridis',ax=ax)
    pg.viewer.mpl.drawSensors(ax, data.sensorPositions(), diam=0.5, color="k")
    ax.set_xlabel('Distance (m)')
    ax.set_ylabel('Elevation (m)')
    ax.set_title('S-wave Vertical Velocity Gradient (m/s/m)')
    
    
    ax3 = fig1.add_subplot(325)
    pg.show(mesh, data=vpvs, coverage=ra.standardizedCoverage(),cMin=1,cMax=5,cmap='plasma',ax=ax3)
    pg.viewer.mpl.drawSensors(ax2, data.sensorPositions(), diam=0.5, color="k")
    ax3.set_xlabel('Distance (m)')
    ax3.set_ylabel('Elevation (m)')
    ax3.set_title('Vp/Vs')
    
    
    ax4 = fig1.add_subplot(326)
    pg.show(mesh, data=pr, coverage=ra.standardizedCoverage(),cMin=0,cMax=0.5,cmap='plasma',ax=ax4)
    pg.viewer.mpl.drawSensors(ax2, data.sensorPositions(), diam=0.5, color="k")
    ax4.set_xlabel('Distance (m)')
    ax4.set_ylabel('Elevation (m)')
    ax4.set_title('Poisson Ratio')
    fig1.tight_layout()
    fig1.savefig('ModelResults.png',dpi=600)
    
    
    
    
    
    
    
    
    fig5 = plt.figure('Poisson Ratio')
    ax4 = fig5.add_subplot(111)
    pg.show(mesh, data=pr, coverage=ra.standardizedCoverage(),cMin=0,cMax=0.5,cmap='plasma',ax=ax4)
    pg.viewer.mpl.drawSensors(ax2, data.sensorPositions(), diam=0.5, color="k")
    ax4.set_xlabel('Distance (m)')
    ax4.set_ylabel('Elevation (m)')
    ax4.set_title('Poisson Ratio')
    #ra.drawRayPaths(ax4,color="k")
    #ra2.drawRayPaths(ax4,color="b")
    fig5.tight_layout()
    fig5.savefig('PoissonsRatio.png',dpi=600)
    
    #ra.showResult(cMin=500,cMax=4500,cmap='jet_r')
    
        
    fig2 = plt.figure('Ray Coverage P-wave')
    vel = ra.paraModel()
    fig2.set_size_inches([8,8])
    ax = fig2.add_subplot(211)
    pg.show(ra.paraDomain, vel, coverage=ra.standardizedCoverage(),cMin=100,cMax=4000,cmap='jet_r',ax=ax)
    pg.viewer.mpl.drawSensors(ax, data.sensorPositions(), diam=0.5, color="k")
    #ra.drawRayPaths(ax,color="k")
    ax.set_xlabel('Distance (m)')
    ax.set_ylabel('Elevation (m)')
    ax.set_title('Ray Paths')
    ax2 = fig2.add_subplot(212)
    tmpCov = ra.coverage()
    tmpCov[np.isneginf(tmpCov)] = 0
    pg.show(mesh, data=tmpCov, coverage=ra.standardizedCoverage(),cMin=0,cMax=10,ax=ax2,cmap='plasma')
    pg.viewer.mpl.drawSensors(ax2, data.sensorPositions(), diam=0.5, color="k")
    ax2.set_xlabel('Distance (m)')
    ax2.set_ylabel('Elevation (m)')
    ax2.set_title('Ray Coverage (# of rays in each cell)')
    fig2.savefig('ModelRayCoverage_pwave.png',dpi=600)
    
            
    fig2 = plt.figure('Ray Coverage S-wave')
    fig2.set_size_inches([8,8])
    ax = fig2.add_subplot(211)
    pg.show(ra2.paraDomain, sVel, coverage=ra2.standardizedCoverage(),cMin=100,cMax=4000,cmap='jet_r',ax=ax)
    pg.viewer.mpl.drawSensors(ax, data2.sensorPositions(), diam=0.5, color="k")
    #ra2.drawRayPaths(ax,color="k")
    ax.set_xlabel('Distance (m)')
    ax.set_ylabel('Elevation (m)')
    ax.set_title('Ray Paths')
    ax2 = fig2.add_subplot(212)
    tmpCov = ra2.coverage()
    tmpCov[np.isneginf(tmpCov)] = 0
    pg.show(mesh, data=tmpCov, coverage=ra2.standardizedCoverage(),cMin=0,cMax=10,ax=ax2,cmap='plasma')
    pg.viewer.mpl.drawSensors(ax2, data.sensorPositions(), diam=0.5, color="k")
    ax2.set_xlabel('Distance (m)')
    ax2.set_ylabel('Elevation (m)')
    ax2.set_title('Ray Coverage (# of rays in each cell)')
    fig2.savefig('ModelRayCoverage_swave.png',dpi=600)
    
    
    
    fig3 = plt.figure('S-wave Model Fits')
    fig3.set_size_inches([11,5])
    ax1 = fig3.add_subplot(121)
    pickTT = np.array(ra2.data('t'))
    modTT = np.array(ra2.inv.response)
    ax1.hist2d(pickTT*1000,(pickTT-modTT)*1000,bins=75,vmin=0,vmax=10,cmap='plasma')
    ax1.plot([0,200],[0,0],'k--',linewidth=2)
    ax1.set_xlabel('Observed Travel-times (ms)')
    ax1.set_ylabel('Residual (Obs - Mod) (ms)')
    
    ax2 = fig3.add_subplot(122)
    ax2.hist2d(pickTT*1000,modTT*1000,bins=100,vmin=0,vmax=10,cmap='plasma')
    ax2.plot([0,200],[0,200],'k',linewidth=2)
    ax2.plot([0,200],[-5,195],'k--',linewidth=2)
    ax2.plot([0,200],[5,205],'k--',linewidth=2)
    ax2.set_xlabel('Observed Travel-times (ms)')
    ax2.set_ylabel('Modelted Travel-times (ms)')
    fig3.savefig('Swave_ModelFits.png',dpi=600)
    
    
    fig3 = plt.figure('P-wave Model Fits')
    fig3.set_size_inches([11,5])
    ax1 = fig3.add_subplot(121)
    pickTT = np.array(ra.data('t'))
    modTT = np.array(ra.inv.response)
    ax1.hist2d(pickTT*1000,(pickTT-modTT)*1000,bins=75,vmin=0,vmax=10,cmap='plasma')
    ax1.plot([0,200],[0,0],'k--',linewidth=2)
    ax1.set_xlabel('Observed Travel-times (ms)')
    ax1.set_ylabel('Residual (Obs - Mod) (ms)')
    
    ax2 = fig3.add_subplot(122)
    ax2.hist2d(pickTT*1000,modTT*1000,bins=100,vmin=0,vmax=10,cmap='plasma')
    ax2.plot([0,200],[0,200],'k',linewidth=2)
    ax2.plot([0,200],[-5,195],'k--',linewidth=2)
    ax2.plot([0,200],[5,205],'k--',linewidth=2)
    ax2.set_xlabel('Observed Travel-times (ms)')
    ax2.set_ylabel('Modelted Travel-times (ms)')
    fig3.savefig('Pwave_ModelFits.png',dpi=600)
    
    
    #S-wave Misfit by Shot
    shotLocs = np.array(ra2.data('s'))
    pickTT = np.array(ra2.data('t'))
    modTT = np.array(ra2.inv.response)
    uniqShots_swave = np.unique(shotLocs) 
    rmsUniqShots_swave = np.zeros(len(uniqShots_swave))
    shotNum = 0
    for shot in uniqShots_swave:
        ind = np.where(shot==shotLocs)[0]
        rmsUniqShots_swave[shotNum] = np.sqrt(np.sum((pickTT[ind]-modTT[ind])**2)/len(ind))*1000
        shotNum = shotNum + 1
    
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
    ax1 = fig4.add_subplot(122)
    ax1.bar(uniqShots_swave,rmsUniqShots_swave,width=5,color='tab:blue')
    ax1.set_xlabel('Shot Location')
    ax1.set_ylabel('RMS Missfit (ms)')
    ax1.set_title('S-wave Inversion')
    
    ax2 = fig4.add_subplot(121)
    ax2.bar(uniqShots_pwave,rmsUniqShots_pwave,width=5,color='tab:blue')
    ax2.set_xlabel('Shot Location')
    ax2.set_ylabel('RMS Missfit (ms)')
    ax2.set_title('P-wave Inversion')
    fig4.tight_layout()
    fig4.savefig('MissfitByShot.png',dpi=600)
    
    # ##OUTPUT
    modTT = np.array(ra2.inv.response)
    sensorPostions_x  =np.array(ra2.data.sensorPositions())[:,0]
    #test = np.column_stack((ra.data('s'),ra.data('g'),modTT,(pickTT-modTT)*1000))
    test = np.column_stack((ra2.data('s'),ra2.data('g'),modTT))
    for i in range(0,len(test)):
        test[i,0] = sensorPostions_x[int(test[i,0])]
        test[i,1] = sensorPostions_x[int(test[i,1])]
    np.savetxt(pickFile2.split('.')[0]+'_Swave_modeledTravelTimes.txt',test,fmt='%10.4f')


    modTT = np.array(ra.inv.response)
    sensorPostions_x  =np.array(ra.data.sensorPositions())[:,0]
    #test = np.column_stack((ra.data('s'),ra.data('g'),modTT,(pickTT-modTT)*1000))
    test = np.column_stack((ra.data('s'),ra.data('g'),modTT))
    for i in range(0,len(test)):
        test[i,0] = sensorPostions_x[int(test[i,0])]
        test[i,1] = sensorPostions_x[int(test[i,1])]
    np.savetxt(pickFile.split('.')[0]+'_Pwave_modeledTravelTimes.txt',test,fmt='%10.4f')