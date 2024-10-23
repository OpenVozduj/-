import numpy as np
import std_atmosphere as sa
import os
import shutil
import scriptsEjekatl as se
from subprocess import call

#%% flow conditions
Re = 1.5e6
Ma = 0.15
c = 1
rho, T, p, cp, mu, Pr, U = sa.properties_M_Re(Re, Ma, c)

#%% Directories
meshDir = 'Meshes'
caseName = 'M015_Re15'
externalDir = '/media/labgpu104/Sukhoi/OpenVozduj/incompressibleFlow/M015_Re15_2/'
rootDir = os.getcwd()

#%% Simulation
nameAirfoils = np.load('names_airfoils.npy').tolist()
alphas = np.load('alphas.npy')
maxIter = 5000
cores = 8
exFoam = '/opt/openfoam11/bin/foamExec'
pRun =' mpirun.openmpi --hostfile machines -np '+str(cores)+' '+exFoam+' foamRun -parallel'

for i in range(3):
    airfoilDir = caseName+'/'+nameAirfoils[i]
    os.mkdir(airfoilDir)
    alphas_i = alphas[i]
    cl = np.array([])
    cd = np.array([])
    pe = np.array([])
    alpha = np.array([])
    for j in range(alphas_i.size):
        alphaDir = airfoilDir+'/alpha_'+str(alphas_i[j])
        os.mkdir(alphaDir)
        meshFile = meshDir+'/'+nameAirfoils[i]+'.msh2'
        shutil.copy(meshFile, alphaDir+'/'+nameAirfoils[i]+'.msh2')
        dir0 = alphaDir+'/0'
        os.mkdir(dir0)
        os.chdir(dir0)
        se.dir_0(rho, T, p, cp, mu, Pr, U, alphas_i[j])
        os.chdir(rootDir)
        dirSystem = alphaDir+'/system'
        os.mkdir(dirSystem)
        os.chdir(dirSystem)
        se.dir_system(cores, maxIter)
        os.chdir(rootDir)
        shutil.copy('machines', alphaDir+'/machines')
        os.chdir(alphaDir)
        call([exFoam, 'gmshToFoam', nameAirfoils[i]+'.msh2'])
        se.boundary()
        se.physical_propertiers()
        se.momentum_transport()
        shutil.move('physicalProperties', 'constant/physicalProperties')
        shutil.move('momentumTransport', 'constant/momentumTransport')
        call([exFoam, 'decomposePar'])
        os.system(pRun)
        call([exFoam, 'reconstructPar'])
        for c in range(cores):
            delProc = 'processor'+str(c)
            shutil.rmtree(delProc)
        cd_a, cl_a = se.read_coeffs()
        os.chdir(rootDir)
        alpha = np.append(alpha, alphas_i[j])
        cl = np.append(cl, cl_a)
        cd = np.append(cd, cd_a)
        pe = np.append(pe, cl_a**1.5/cd_a)
        if pe.size > 1:
            if pe[-1] < pe[-2]:
                break
    coeffs = np.row_stack((alpha, cl))
    coeffs = np.row_stack((coeffs, cd))
    coeffs = np.row_stack((coeffs, pe))
    np.save(airfoilDir+'/coeffs_'+nameAirfoils[i]+'.npy', coeffs)
