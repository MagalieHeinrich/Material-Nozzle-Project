Material Nozzle Project: JetRIS Simulation

This repository contains a numerical framework for simulating compressible fluid dynamics through the JetRIS apparatus. 
It is designed to model the expansion of Argon gas from a high-pressure stagnation chamber into a low-pressure ambient jet chamber.

The simulation leverages Trixi.jl for the numerical solvers and HOHQMesh.jl for unstructured mesh generation. 

---

Physics & Mathematical Model

The code models monatomic Argon gas gamma = 5/3 using a coupled hyperbolic-parabolic system:
Convective Fluxes: 2D Compressible Euler Equations.
Diffusive Fluxes: 2D Compressible Navier-Stokes Diffusion.

Core Parameters & Constants
Viscosity (mu): 22.7241 e-6 Pa*s (based on Lemmon & Jacobsen, 2003).
Prandtl Number (Pr): 2/3 (standard for monoatomic gases).
Molar Mass (M_m): 0.039948 kg/mol.
Stagnation State (Gas cell): Pressure P_{stagn} = 8000 Pa, Temperature T_stagn} = 300 K.
Ambient State (Jet chamber): Pressure $P_{amb} = 5 Pa.

---

Solver Configuration (Trixi.jl)

The numerical scheme uses a Discontinuous Galerkin Spectral Element Method (DGSEM) with a default polynomial degree of `my_poldeg = 2`. 

To handle strong shocks and expansions through the aperture without crashing, the solver utilizes:
Volume Integral: `VolumeIntegralShockCapturingHG` (Hennemann-Gassner indicator) tracking `density_pressure` variations with a max smoothing blending factor of alpha_{max} = 0.8.
Flux Solvers: `flux_ranocha_turbo` for volume fluxes combined with a robust Lax-Friedrichs / HLLC surface flux mixture.
Limiter: `PositivityPreservingLimiterZhangShu` targeting negative densities and pressures at a strict threshold of 1.0 e-12.
Time Integration: Strong Stability Preserving Runge-Kutta `SSPRK43` via `OrdinaryDiffEqLowStorageRK.jl` or newly `OrdinaryDiffEqSSPRK`.

---

Mesh Generation Architecture (HOHQMesh.jl)
The project includes two distinct mesh configurations inside the source layout. 
Both scripts generate 2D planar meshes in an ABAQUS (.inp) format compatible with Trixi's P4estMesh reader.

1. Default Mesh: JetRIS_verysimp_ap_new.inp (from Newaperture_mesh.jl)
A simplified, sharp-cornered block domain scaling structures via centimeter adjustments. It features targeted local path refinements near the aperture boundary (ref_ap) and wall boundaries (ref_walljetchamber).

2. Alternative Complex Mesh: JetRIS_verysimp_ap5.inp (from JetRIS_verysimp_ap5.jl)
A complex, curved-funnel mesh tracking customized circular arcs to smooth out the entry throat geometry. This mesh has problems with elongation of the domain.

Mesh Configuration Warning:

If you switch the simulation script to use JetRIS_verysimp_ap5.inp, you must completely map the updated boundary dictionaries (my_boundary_conditions and my_bcs_parabolic). 
Every named geometric edge segment created by HOHQMesh must match a valid boundary routine assignment in Trixi, or the initialization loop will fail.

---

How to Run the Simulations:

Ensure your Julia environment has the required packages installed. Open your Julia REPL and run:
using Pkg
Pkg.add(["Trixi", "OrdinaryDiffEq", "StaticArrays", "OrdinaryDiffEqLowStorageRK", "OrdinaryDiffEqSSPRK", "Trixi2Vtk", "LinearAlgebra", "GLMakie", "HOHQMesh"])

Generate the Mesh (run the mesh code) and then run the simulation

Post-Processing & Visualization:
The script terminates by running a trixi2vtk command, automatically processing all raw data files located inside h5_raw_data/ and storing ready-to-view .vtu spatial meshes.
You can open these .vtu files in paraview.
