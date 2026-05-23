"
The simulation code in Ap_NvSt_newIC.jl works with the mesh that Newaperture_mesh.jl creates. 
If you want to use the mesh that is created with this code here you must put JetRIS_verysimp_ap5.inp in mesh_file = joinpath(@__DIR__, .... instead of JetRIS_verysimp_ap_new.inp in the file Ap_NvSt_newIC.jl
You must also completely change the boundary dictionaries (yes, both the hyperbolic and the parabolic) in the file Ap_NvSt_newIC.jl. 
Every boundary that is defined in JetRIS_verysimp_ap5.inp must get assigned a hyperbolic and parabolic meaning in Ap_NvSt_newIC.jl.

Warning: this mesh code has problems regarding elongating the domain.
"

using GLMakie, HOHQMesh

JetRIS_verysimp_ap5 = newProject("JetRIS_verysimp_ap5", joinpath(@__DIR__, "meshes"))
setPolynomialOrder!(JetRIS_verysimp_ap5, 1)
setMeshFileFormat!(JetRIS_verysimp_ap5, "ABAQUS")
#setMeshFileFormat!(JetRIS_verysimp_ap5, "ISM-V2")

f=0.37#Change to control the x-diameter of the right chamber
g=0.1#Change to control the x-diameter of the left chamber
d=0.37#Make bigger to set the lower left corner lower (y-direction)
lower_left = [-0.215, -0.132-d, 0.0]
#Comment the 18 mm option and uncomment an other option for a finer mesh
#spacing = [0.018, 0.018, 0.0] 
#num_intervals = [77, 14, 0] 
spacing = [0.017, 0.017, 0.0]
num_intervals = [90, 45, 0] 
#spacing = [0.003, 0.003, 0.0] 
#num_intervals = [462, 80, 0] 
#spacing = [0.006, 0.006, 0.0] 
#num_intervals = [231, 40, 0] 
#spacing = [0.012, 0.012, 0.0] 
#num_intervals = [154, 22, 0] 
#spacing = [0.006, 0.006, 0.0] 
#num_intervals = [308, 44, 0] 
addBackgroundGrid!(JetRIS_verysimp_ap5, lower_left, spacing, num_intervals)
plotProject!(JetRIS_verysimp_ap5, GRID)

#define outer boundary curves (all in meters)

#capping circle:
Ax_lw = 0.088829 - 0.005 * sind(270.000-315.251)
Ay_lw = (-0.001-0.005*cosd(270.000-315.251))
Bx_lw = 0.0934+0.0006
By_lw = (-0.001-0.005*cosd(270.000-315.251))
radius= (Bx_lw - Ax_lw)/(cosd(30))
Mx_lw = Bx_lw
My_lw = By_lw-radius/2

#capping circle:
Bx_u = 0.088229+0.0006+0.005*sind(90-44.749)
By_u = 0.001+0.005*cosd(90-44.749)
Ax_u = 0.0934+0.0006
Ay_u = 0.001+0.005*cosd(90-44.749)
radius= (Bx_u - Ax_u)/(cosd(180+30))
Mx_u = Ax_u
My_u = Ay_u+radius/2

b=My_lw+radius*sind(90)+0.0005#Change the +0.0005 to a larger number to make the aperture diameter larger

c=0.27#Change to make the left chamber shorter 0.114 works
a=0.57#Change to make the right chamber longer 0.114 works
#a=c must be fulfilled in orfer to generate a useful mesh.

l_ulw_vwall_inflow = newEndPointsLineCurve("l_ulw_vwall_inflow",[-0.215+c,0.0,0.000],[-0.215+c,-0.126+g,0.000])
l_lw_hwall_lower = newEndPointsLineCurve("l_lw_hwall_lower",[-0.215+c,-0.126+g,0.000],[0.086,-0.126+g,0.000])
r_lw_vwall_leftorange = newEndPointsLineCurve("r_lw_vwall_leftorange",[0.086,-0.126+g,0.000],[0.086,-0.006-b-0.001,0.000])
r_lw_circwall_funnel_third = newCircularArcCurve("r_lw_circwall_funnel_third",[0.086+0.001,-0.006-b-0.001,0.000],0.001,180.000,90.000,"degrees")
r_lw_hwall_funnel_sixth = newEndPointsLineCurve("r_lw_hwall_funnel_sixth",[0.086+0.001,-0.006-b,0.000],[0.088229+0.0006,-0.006-b,0.000]) 
r_lw_circwall_funnel_seventh = newCircularArcCurve("r_lw_circwall_funnel_seventh",[0.088229+0.0006,-0.001-b,0.000],0.005,270.000,315.251,"degrees")  
r_lw_diagwall_funnel_eighth = newCircularArcCurve("r_lw_diagwall_funnel_eighth", [Mx_lw,My_lw-b,0.000],radius,150.000,90.000,"degrees")
r_lw_vwall_funnel_nineth = newEndPointsLineCurve("r_lw_vwall_funnel_nineth",[Mx_lw,My_lw+radius*sind(90)-b,0.000],[0.0934+0.0006,-0.126-d,0.000]) #  
r_lw_hwall_lower = newEndPointsLineCurve("r_lw_hwall_lower",[0.0934+0.0006,-0.126-d,0.000],[0.414+a,-0.126-d,0.000]) #  
r_lwu_vwall_outflow = newEndPointsLineCurve("r_lwu_vwall_outflow",[0.414+a,-0.126-d,0.000],[0.414+a,0.0,0.000])

#symmetry line
sym=newEndPointsLineCurve(":symmetry",[0.414+a,0.0,0.0],[-0.215+c,0.0,0.0])

# Add outer boundary curves
addCurveToOuterBoundary!(JetRIS_verysimp_ap5,l_ulw_vwall_inflow)
addCurveToOuterBoundary!(JetRIS_verysimp_ap5,l_lw_hwall_lower)
addCurveToOuterBoundary!(JetRIS_verysimp_ap5,r_lw_vwall_leftorange)
addCurveToOuterBoundary!(JetRIS_verysimp_ap5,r_lw_circwall_funnel_third)
addCurveToOuterBoundary!(JetRIS_verysimp_ap5,r_lw_hwall_funnel_sixth)
addCurveToOuterBoundary!(JetRIS_verysimp_ap5,r_lw_circwall_funnel_seventh)
addCurveToOuterBoundary!(JetRIS_verysimp_ap5,r_lw_diagwall_funnel_eighth)
addCurveToOuterBoundary!(JetRIS_verysimp_ap5,r_lw_vwall_funnel_nineth)
addCurveToOuterBoundary!(JetRIS_verysimp_ap5,r_lw_hwall_lower)
addCurveToOuterBoundary!(JetRIS_verysimp_ap5,r_lwu_vwall_outflow)
addCurveToOuterBoundary!(JetRIS_verysimp_ap5,sym)

#Refinement
#ref_lw_funnelhupperside=newRefinementLine("line", "smooth",  [ 0.086,-0.006-b, 0.000],[ 0.088229,-0.006-b, 0.000],0.002, 2.0/1000)  
#ref_lw_funnelcircupperside_first=newRefinementLine("line", "smooth",  [ 0.088229,-0.006-b, 0.000],[ 0.089210,-0.005903-b, 0.000],0.001, 2.0/1000)  
#ref_lw_funnelcircupperside_second=newRefinementLine("line", "smooth",  [ 0.089210,-0.005903-b, 0.000],[ 0.090153,-0.005615-b, 0.000],0.001, 2.0/1000)  
#ref_lw_funnelcircupperside_third=newRefinementLine("line", "smooth",  [ 0.090153,-0.005615-b, 0.000],[ 0.091020,-0.005148-b, 0.000],0.001, 2.0/1000)  
#ref_lw_funnelcircupperside_fourth=newRefinementLine("line", "smooth",  [ 0.091020,-0.005148-b, 0.000],[ 0.091780,-0.004520-b, 0.000],0.001, 2.0/1000)  
#ref_lw_funneldiagtoaperture=newRefinementLine("line", "smooth",  [ 0.091780,-0.004520-b, 0.000],[ 0.0934,-0.003-b, 0.000],0.002, 2.0/1000)  
#ref_r_lw_diagwall_funnel_eighth = newRefinementLine("line", "smooth",[0.088229+0.0006 - 0.005 * sind(270.000-315.251),(-0.001-0.005*cosd(270.000-315.251))-b,0.000],[0.0934+0.0006,(-0.001-0.005*cosd(270.000-315.251))-b,0.000],0.0005, 8.0/1000)  
#ref_r_u_diagwall_funnel_second = newRefinementLine("line", "smooth",[0.0934+0.0006,0.0045-b,0.000],[0.088229+0.0006+0.005*sind(90-44.749),0.001+0.005*cosd(90-44.749)-b,0.000],0.0005, 8.0/1000)  
ref_out = newRefinementLine("line", "smooth",[0.094,-0.001,0.000],[0.094,-0.002,0.000],0.6*0.0005, 3.0/1000)
ref_middleline = newRefinementLine("line","smooth",[0.092,0.0,0.0],[0.095,0.0,0.0],0.5*0.0005, 4.0/1000)

#ref_r_lw_circwall_funnel_third = newRefinementCenter("center","smooth",[0.086+0.001,-0.006-b-0.001,0.000],0.5*0.0005, 4.0/1000)

#add the Refinement
#addRefinementRegion!(JetRIS_verysimp_ap5,ref_lw_funnelhupperside)
#addRefinementRegion!(JetRIS_verysimp_ap5,ref_lw_funnelcircupperside_first)
#addRefinementRegion!(JetRIS_verysimp_ap5,ref_lw_funnelcircupperside_second)
#addRefinementRegion!(JetRIS_verysimp_ap5,ref_lw_funnelcircupperside_third)
#addRefinementRegion!(JetRIS_verysimp_ap5,ref_lw_funnelcircupperside_fourth)
#addRefinementRegion!(JetRIS_verysimp_ap5,ref_lw_funneldiagtoaperture)
addRefinementRegion!(JetRIS_verysimp_ap5,ref_out)
addRefinementRegion!(JetRIS_verysimp_ap5,ref_middleline)

#addRefinementRegion!(JetRIS_verysimp_ap5,ref_r_lw_circwall_funnel_third)


generate_mesh(JetRIS_verysimp_ap5)
