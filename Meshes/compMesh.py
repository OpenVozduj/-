import numpy as np
import cst
from scipy.interpolate import InterpolatedUnivariateSpline
import gmsh

#%% real subsonics airfoils 
airfoils = np.load('airfoilsCST_sub.npy')
nameAirfoils = np.load('airfoilsNames_sub.npy')
nameAirfoils = nameAirfoils.tolist()
for i in range(len(nameAirfoils)):
    nameAirfoils[i] = nameAirfoils[i].replace('airfoilsCST_sub/', '')
    nameAirfoils[i] = nameAirfoils[i].replace('_CST6.json', '')
#%% real supercritics airfoils
# airfoils = np.load('airfoilsCST_trans.npy')
# nameAirfoils = np.load('airfoilsNames_trans.npy')
# nameAirfoils = nameAirfoils.tolist()
# for i in range(len(nameAirfoils)):
#     nameAirfoils[i] = nameAirfoils[i].replace('airfoilsCST_trans/', '')
#     nameAirfoils[i] = nameAirfoils[i].replace('_CST6.json', '')
#%% fake subsonics airfoils
# airfoils = np.load('airfoilsCSTgan_sub2.npy')
# nameAirfoils = np.load('AirfoilsNamesGan_sub2.npy')
# nameAirfoils = nameAirfoils.tolist()
# for i in range(len(nameAirfoils)):
#     nameAirfoils[i] = nameAirfoils[i].replace('airfoil_gant_', 'gans_')
#%%
# airfoils = np.load('airfoilsCSTgan_trans2.npy')
# nameAirfoils = np.load('AirfoilsNamesGan_trans2.npy')
# nameAirfoils = nameAirfoils.tolist()
# for i in range(len(nameAirfoils)):
#     nameAirfoils[i] = nameAirfoils[i].replace('airfoil_gant_', 'gant_')

#%% Optimal airfoils
# airfoils = np.load('nsga2foils.npy')
# nameAirfoils = ['kfoil2_1', 'kfoil2_2', 'kfoil2_3', 'kfoil2_4', 'kfoil2_5']

#%% Mesher
for airfoil in range(len(nameAirfoils)):
    #%% preprocessing of geometry
    x, yu, yl = cst.cstN6(airfoils[airfoil])
    cc = 0.15
    Pu = InterpolatedUnivariateSpline(x, yu, k=1)
    pcu = Pu(cc) 
    Pl = InterpolatedUnivariateSpline(x, yl, k=1)
    pcl = Pl(cc)
    for i in range(x.size):
        if x[i]>cc:
            I = i-1
            break
        else:
            pass
    x_te = x[I+1:-1]
    yu_te = yu[I+1:-1]
    yl_te = yl[I+1:-1]
    x_le = np.concatenate((x[1:I][::-1], x[:I]), axis=0)
    y_le = np.concatenate((yu[1:I][::-1], yl[:I]), axis=0)
    nte = x_te.size
    nle = x_le.size    
    #%% control volume
    R = 20.0
    L = 45.0
    pay = R*np.sin(60*np.pi/180)
    pax = R*np.cos(60*np.pi/180)
    numid = nameAirfoils[airfoil]
    gmsh.initialize()
    gmsh.model.add(numid)
    # Geometry
    lc = 1
    #Points
    gmsh.model.geo.addPoint(1.0, 0.0, 0.0, lc, 1)
    gmsh.model.geo.addPoint(cc, pcu, 0.0, lc, 2)
    gmsh.model.geo.addPoint(cc, pcl, 0.0, lc, 3)
    for i in range(nle):
        gmsh.model.geo.addPoint(x_le[i], y_le[i], 0.0, lc, i+4)
    for i in range(nte):
        gmsh.model.geo.addPoint(x_te[i], yu_te[i], 0.0, lc, i+nle+4)
    for i in range(nte):
        gmsh.model.geo.addPoint(x_te[i], yl_te[i], 0.0, lc, i+nle+nte+4)
    gmsh.model.geo.addPoint(1.0, R, 0.0, lc, 500)
    gmsh.model.geo.addPoint(1-pax, pay, 0.0, lc, 501)
    gmsh.model.geo.addPoint(1-pax, -pay, 0.0, lc, 502)
    gmsh.model.geo.addPoint(1.0, -R, 0.0, lc, 503)
    gmsh.model.geo.addPoint(L, -R, 0.0, lc, 504)
    gmsh.model.geo.addPoint(L, 0.0, 0.0, lc, 505)
    gmsh.model.geo.addPoint(L, R, 0.0, lc, 506)
    #Curves
    tagsSpl1 = [2]
    for i in range(4, nle+4):
        tagsSpl1.append(i)
    tagsSpl1.append(3)
    gmsh.model.geo.addSpline(tagsSpl1, 1)
    tagsSpl2 = [2]
    for i in range(4+nle, nte+nle+4):
        tagsSpl2.append(i)
    tagsSpl2.append(1)
    gmsh.model.geo.addSpline(tagsSpl2, 2)
    tagsSpl3 = [3]
    for i in range(4+nle+nte, 2*nte+nle+4):
        tagsSpl3.append(i)
    tagsSpl3.append(1)
    gmsh.model.geo.addSpline(tagsSpl3, 3)
    gmsh.model.geo.addCircleArc(501, 1, 500, 4)
    gmsh.model.geo.addCircleArc(501, 1, 502, 5)
    gmsh.model.geo.addCircleArc(502, 1, 503, 6)
    gmsh.model.geo.addLine(500, 506, 7)
    gmsh.model.geo.addLine(1, 505, 8)
    gmsh.model.geo.addLine(503, 504, 9)
    gmsh.model.geo.addLine(2, 501, 10)
    gmsh.model.geo.addLine(1, 500, 11)
    gmsh.model.geo.addLine(505, 506, 12)
    gmsh.model.geo.addLine(3, 502, 13)
    gmsh.model.geo.addLine(1, 503, 14)
    gmsh.model.geo.addLine(505, 504, 15)
    #Loops
    gmsh.model.geo.addCurveLoop([1, 13, -5, -10], 1)
    gmsh.model.geo.addCurveLoop([2, 11, -4, -10], 2)
    gmsh.model.geo.addCurveLoop([3, 14, -6, -13], 3)
    gmsh.model.geo.addCurveLoop([8, 12, -7, -11], 4)
    gmsh.model.geo.addCurveLoop([8, 15, -9, -14], 5)
    #Surfaces
    for i in range(5):
        gmsh.model.geo.addPlaneSurface([i+1], i+1)    
    #Mesh
    elem_a = 50 #50
    grad_a = 5 #5
    elem_w = 80 #81
    grad_w = 1.15 #1.11
    elem_l = 75 #81
    grad_l = 1.0872 #1.06
    elem_s = 30 #30
    gmsh.model.geo.mesh.setTransfiniteCurve(1, elem_a, meshType='Bump', coef=grad_a)
    gmsh.model.geo.mesh.setTransfiniteCurve(2, elem_s)
    gmsh.model.geo.mesh.setTransfiniteCurve(3, elem_s)
    gmsh.model.geo.mesh.setTransfiniteCurve(4, elem_s)
    gmsh.model.geo.mesh.setTransfiniteCurve(5, elem_a, meshType='Bump', coef=grad_a)
    gmsh.model.geo.mesh.setTransfiniteCurve(6, elem_s)
    gmsh.model.geo.mesh.setTransfiniteCurve(7, elem_l, meshType='Progression', coef=grad_l)
    gmsh.model.geo.mesh.setTransfiniteCurve(8, elem_l, meshType='Progression', coef=grad_l)
    gmsh.model.geo.mesh.setTransfiniteCurve(9, elem_l, meshType='Progression', coef=grad_l)
    gmsh.model.geo.mesh.setTransfiniteCurve(10, elem_w, meshType='Progression', coef=grad_w)
    gmsh.model.geo.mesh.setTransfiniteCurve(11, elem_w, meshType='Progression', coef=grad_w)
    gmsh.model.geo.mesh.setTransfiniteCurve(12, elem_w, meshType='Progression', coef=grad_w)
    gmsh.model.geo.mesh.setTransfiniteCurve(13, elem_w, meshType='Progression', coef=grad_w)
    gmsh.model.geo.mesh.setTransfiniteCurve(14, elem_w, meshType='Progression', coef=grad_w)
    gmsh.model.geo.mesh.setTransfiniteCurve(15, elem_w, meshType='Progression', coef=grad_w)
    # Surfaces
    for i in range(5):
        gmsh.model.geo.mesh.setTransfiniteSurface(i+1)
    for i in range(5):
        gmsh.model.geo.mesh.setRecombine(2, i+1)
    #Volumes
    lay = [(2,1),
            (2,2),
            (2,3),
            (2,4),
            (2,5)]
    gmsh.model.geo.extrude(lay, 0, 0, 1, numElements=[1], recombine=True)
    gmsh.model.geo.synchronize()    
    #Physical groups
    #CV
    vc_phys_tag = 1001
    vc_num_tags = [1, 2, 3, 4, 5]
    gmsh.model.addPhysicalGroup(3, vc_num_tags, tag=vc_phys_tag)
    gmsh.model.setPhysicalName(3, vc_phys_tag, "CV")
    #Airfoil
    wall_tag = 2001
    wall_tags = [24, 46, 68]
    gmsh.model.addPhysicalGroup(2, wall_tags, tag=wall_tag)
    gmsh.model.setPhysicalName(2, wall_tag, "airfoil")
    #FrontBack
    front_tag = 3001
    front_tags = [1, 2, 3, 4, 5, 37, 59, 81, 103, 125]
    gmsh.model.addPhysicalGroup(2, front_tags, tag=front_tag)
    gmsh.model.setPhysicalName(2, front_tag, "frontBack")
    #Inlet
    inlet_tag = 4001
    inlet_tags = [32, 54, 76, 98, 120]
    gmsh.model.addPhysicalGroup(2, inlet_tags, tag=inlet_tag)
    gmsh.model.setPhysicalName(2, inlet_tag, "inlet")
    #Outlet
    outlet_tag = 5001
    outlet_tags = [94, 116]
    gmsh.model.addPhysicalGroup(2, outlet_tags, tag=outlet_tag)
    gmsh.model.setPhysicalName(2, outlet_tag, "outlet")
    #Mesh creation
    gmsh.model.mesh.generate(3)
    gmsh.write(numid+".msh2")
    #Sys
    # if '-nopopup' not in sys.argv:
    #     gmsh.fltk.run()
    gmsh.finalize()
