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
    sim.PNTol = 1e-3

    sim.add_shell_with_scale_3D("input/dress/high_cloth_wraping.obj", Vector3d(0.0, 0.0, -0.0), Vector3d(0.01, 0.01, 0.01),\
        Vector3d(0, 0, 0), Vector3d(1, 0, 0), 0)
    
    solve_static = True
    sim.muComp = StdVectorXd([0, 0, sim.mu,  0, 0, sim.mu,  sim.mu, sim.mu, 0.1])

    sim.dt = 0.002
    sim.frame_dt = 0.002
    sim.frame_num = 0
    sim.use_s2 = True
    sim.use_dist = True
    sim.use_populate = True
    sim.gravity = Vector3d(0, 0, 0) 
    sim.staticSolve = solve_static
    sim.smock = True
    sim.smock_cons = 1.0
    sim.uniform_stitching_ratio_smock = 1.0
    sim.withCollision = True

    # wb front data
    stage_size = 12911
    stage_start = 135415
    start = 145825
    end = 139400
    scale = 1.9 * 0.01

    sim.populate(start, end, stage_start, stage_size, scale)

    # wb back
    stage_size = 12964
    stage_start = 15715
    start = 17334
    end = 24288

    sim.populate(start, end, stage_start, stage_size, scale)


    if(solve_static):
        sim.dt = 0.01
        sim.smock_cons = 1.0
        sim.k_stitch = 1e7
        sim.frame_num = 1
        sim.staticSolve = True
        sim.PNTol = 1e-4
        sim.withCollision = False
    
    # populate before sim init
    sim.initialize(sim.cloth_density_iso[clothI], sim.cloth_Ebase_iso[clothI] * membEMult,
    sim.cloth_nubase_iso[clothI], sim.cloth_thickness_iso[clothI], 0, sim.smock)
    sim.bendingStiffMult = bendEMult / membEMult
    sim.kappa_s = Vector2d(0, 0)
    
    offset_vec = -0.008
    sim.Offset_stitching(offset_vec, 72 * 2)

    sim.initialize_OIPC(1e-3, 0)

    sim.run()
