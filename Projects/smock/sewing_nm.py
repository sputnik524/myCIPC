import sys
sys.path.insert(0, "../../Python")
import Drivers
from JGSL import *

if __name__ == "__main__":
    sim = Drivers.FEMDiscreteShellBase("double", 3)

    algI = 1
    if len(sys.argv) > 1:
        algI = int(sys.argv[1])

    clothI = 0
    if len(sys.argv) > 2:
        clothI = int(sys.argv[2])

    size = '85K'
    if len(sys.argv) > 3:
        size = sys.argv[3]
    
    membEMult = 1.0
    if len(sys.argv) > 4:
        membEMult = float(sys.argv[4])
    
    bendEMult = 100
    if len(sys.argv) > 5:
        bendEMult = float(sys.argv[5])

    sim.mu = 0.0

    sewing_both_sides = True
    
    # 1. add the smocking garment
    sim.add_shell_with_scale_3D_smock("input/S4.obj","input/S1.obj", Vector3d(0, 0, 0), Vector3d(0.05, 0.05, 0.05),\
        Vector3d(0, 0, 0), Vector3d(0, 0, 1), -90) # make sure the planar is orthogonal to the ground
    
    # set the upper edge of smocking mesh to DBC with range input
    if(not sewing_both_sides):
        sim.set_DBC_with_range(Vector3d(-0.1, -0.1, -0.1), Vector3d(1.1, 1.1, 1.1), 
            Vector3d(0, 0.0, 0), Vector3d(0, 0, 0), Vector3d(0, 1, 0), 0, Vector4i(0, 0, 10 * sim.fine_mesh_res, -1))

    # 2. add the garment with identical size that to be sewed 
    sim.add_shell_with_scale_3D("input/S3.obj", Vector3d(0, -1.1, 0), \
        Vector3d(0.05, 0.05, 0.05),Vector3d(0, 0, 0), Vector3d(0, 0, 1), -90)
    
    if sewing_both_sides:
        sim.add_shell_with_scale_3D("input/S3.obj", Vector3d(0, 1.1, 0), \
            Vector3d(0.05, 0.05, 0.05),Vector3d(0, 0, 0), Vector3d(0, 0, 1), -90)
        
        sim.set_DBC_with_range(Vector3d(-0.1, -0.1, -0.1), Vector3d(1.1, 1.1, 1.1), 
            Vector3d(0, 0.0, 0), Vector3d(0, 0, 0), Vector3d(0, 1, 0), 0, Vector4i(2 * sim.fine_mesh_res * sim.fine_mesh_res - sim.smock_size, 0, 2 * sim.fine_mesh_res * sim.fine_mesh_res + sim.fine_mesh_res - sim.smock_size, -1))
    
    # 3. add sewing info
    sim.add_stitching(both_sides = sewing_both_sides)

    sim.muComp = StdVectorXd([0, 0, sim.mu,  0, 0, sim.mu,  sim.mu, sim.mu, 0.1])

    sim.dt = 0.04
    sim.frame_dt = 0.04
    sim.frame_num = 100
    sim.withCollision = True
    # sim.smock = True
    # sim.smock_cons = 0.2 
    sim.uniform_stitching_ratio = 0.25
    sim.if_contact = False
    # sim.gravity = Vector3d(0, 0, 0) # cancel this term after setting DBC
    # sim.staticSolve = True

    # density, E, nu, thickness, initial displacement case
    if algI == 0:
        # iso
        sim.initialize(sim.cloth_density_iso[clothI], sim.cloth_Ebase_iso[clothI] * membEMult,
            sim.cloth_nubase_iso[clothI], sim.cloth_thickness_iso[clothI], 0)
        sim.bendingStiffMult = bendEMult / membEMult
        sim.kappa_s = Vector2d(1e3, 0)
        sim.s = Vector2d(sim.cloth_SL_iso[clothI], 0)
    elif algI == 1:
        # iso, no strain limit
        sim.initialize(sim.cloth_density_iso[clothI], sim.cloth_Ebase_iso[clothI] * membEMult,
            sim.cloth_nubase_iso[clothI], sim.cloth_thickness_iso[clothI], 0, smock = sim.smock, filepath_smock = "input/S2.obj", filepath_smock_pattern = "input/S1.obj")
        sim.bendingStiffMult = bendEMult / membEMult
        sim.kappa_s = Vector2d(0, 0)

    sim.initialize_OIPC(1e-3, 0)

    sim.run()
