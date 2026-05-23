using GLMakie, HOHQMesh

JetRIS_verysimp_ap_new = newProject("JetRIS_verysimp_ap_new", joinpath(@__DIR__, "meshes"))
setPolynomialOrder!(JetRIS_verysimp_ap_new, 1)
setMeshFileFormat!(JetRIS_verysimp_ap_new, "ABAQUS")
#setMeshFileFormat!(JetRIS_verysimp_ap_new, "ISM-V2")


#define dimensions
lefft = -0.040 *100
righht = 0.600 *100
lowergascell = -0.125 *100
lowerjetchamber = -0.125 *100
rightgascell = 0.0 *100
leftjetchamber = 0.010 *100
ap = 0.0 *100
sym = 0.0005 *100

h = 0.015 *100
grid_x0 = lefft - 0.02 *100
grid_y0 = min(lowergascell, lowerjetchamber) - 0.02 *100
lower_left = [grid_x0, grid_y0, 0.0]
spacing = [h, h, 0.0]
num_intervals = [80, 30, 0] 
addBackgroundGrid!(JetRIS_verysimp_ap_new, lower_left, spacing, num_intervals)
plotProject!(JetRIS_verysimp_ap_new, GRID)


#define outer boundary curves (all in meters)  
inflow = newEndPointsLineCurve("inflow",[lefft,sym,0.000],[lefft,lowergascell,0.000])
downgascell = newEndPointsLineCurve("downgascell",[lefft,lowergascell,0.000],[rightgascell,lowergascell,0.000])
wallgascell = newEndPointsLineCurve("wallgascell",[rightgascell,lowergascell,0.000],[rightgascell,ap,0.000])
aperture = newEndPointsLineCurve("aperture",[rightgascell,ap,0.000],[leftjetchamber,ap,0.000])
walljetchamber = newEndPointsLineCurve("walljetchamber",[leftjetchamber,ap,0.000],[leftjetchamber,lowerjetchamber,0.000])
downjetchamber = newEndPointsLineCurve("downjetchamber",[leftjetchamber,lowerjetchamber,0.000],[righht,lowerjetchamber,0.000])
outflow = newEndPointsLineCurve("outflow",[righht,lowerjetchamber,0.000],[righht,sym,0.000])

#symmetry line
sym=newEndPointsLineCurve(":symmetry",[righht,sym,0.000],[lefft,sym,0.000])

# Add outer boundary curves
addCurveToOuterBoundary!(JetRIS_verysimp_ap_new,inflow)
addCurveToOuterBoundary!(JetRIS_verysimp_ap_new,downgascell)
addCurveToOuterBoundary!(JetRIS_verysimp_ap_new,wallgascell)
addCurveToOuterBoundary!(JetRIS_verysimp_ap_new,aperture)
addCurveToOuterBoundary!(JetRIS_verysimp_ap_new,walljetchamber)
addCurveToOuterBoundary!(JetRIS_verysimp_ap_new,downjetchamber)
addCurveToOuterBoundary!(JetRIS_verysimp_ap_new,outflow)
addCurveToOuterBoundary!(JetRIS_verysimp_ap_new,sym)

#Refinement
ref_ap = newRefinementLine("ref_ap", "smooth", [rightgascell-0.005 *100, ap, 0.0], [leftjetchamber+0.005 *100, ap, 0.0], 0.001 *100, 4.0/1000 *100) 
ref_wallgescell = newRefinementLine("ref_wallges", "smooth", [rightgascell, 0.0, 0.0], [rightgascell, lowergascell, 0.0], 0.002 *100, 2.0/1000 *100) 
ref_walljetchamber = newRefinementLine("ref_walljet", "smooth", [leftjetchamber, 0.002 *100, 0.0], [leftjetchamber, lowerjetchamber, 0.0], 0.0005 *100, 2.0/1000 *100)
ref_out = newRefinementCenter("center", "smooth", [ 0.010 *100,0.0, 0.000],0.001 *100, 8.0/1000 *100) 


#add the Refinement
addRefinementRegion!(JetRIS_verysimp_ap_new,ref_ap)
addRefinementRegion!(JetRIS_verysimp_ap_new,ref_wallgescell)
addRefinementRegion!(JetRIS_verysimp_ap_new,ref_walljetchamber)
addRefinementRegion!(JetRIS_verysimp_ap_new,ref_out)

generate_mesh(JetRIS_verysimp_ap_new)
