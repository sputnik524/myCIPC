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
    sim.PNTol = 2e-3

    # add stage garment
    sim.add_shell_with_scale_3D("input/dress/high_cloth_wraping_2.obj", Vector3d(0.0, 0.0, -0.0), Vector3d(1.0, 1.0, 1.0),\
        Vector3d(0, 0, 0), Vector3d(1, 0, 0), 0, True)

    # add mannequin
    meshCounter = sim.add_shell_with_scale_3D("input/f_average_A40.obj", Vector3d(0.0, -0.8, 0.0), \
        Vector3d(0.9, 0.95, 0.85),Vector3d(0, 0, 0), Vector3d(0, 0, 1), 0) 
    
    sim.set_DBC_with_range(Vector3d(-0.1, -0.1, -0.1), Vector3d(1.1, 1.1, 1.1), 
            Vector3d(0, 0.0, 0), Vector3d(0, 0, 0), Vector3d(0, 1, 0), 0, meshCounter)
    
    solve_static = False
    sim.muComp = StdVectorXd([0, 0, sim.mu,  0, 0, sim.mu,  sim.mu, sim.mu, 0.1])

    sim.dt = 0.003
    sim.frame_dt = 0.003
    sim.frame_num = 20
    sim.k_stitch = 1e3
    sim.use_s2 = True
    sim.use_dist = True
    sim.use_populate = True
    sim.gravity = Vector3d(0, -9.8, 0) 
    sim.staticSolve = solve_static
    sim.smock = True
    sim.smock_cons = 5.0
    sim.uniform_stitching_ratio_smock = 1.0
    sim.withCollision = True

    sim.load_frame("input/dress/shell_start_4.obj", True) 
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

    # left sleeve front 
    start = 95561
    end = 79153
    stage_start = 78174
    stage_size = 106821 - 78175 + 1
    x_axis = Vector3d(-56.063 + 59.3825, 5.46448 - 2.67906, 0.0) 
    y_axis = Vector3d(-98.9106 + 97.8388, 2.59655 - 1.31922, 0.0)

    sim.populate(start, end, stage_start, stage_size, scale, x_axis, y_axis)

    # left sleeve back 
    start = 243107
    end = 252690
    stage_start = 243072
    stage_size = 271279 - 243073 + 1
    # x_axis = Vector3d(-56.063 + 59.3825, 5.46448 - 2.67906, 0.0) 
    # y_axis = Vector3d(-98.9106 + 97.8388, 2.59655 - 1.31922, 0.0)
    sim.populate(start, end, stage_start, stage_size, scale, x_axis, y_axis)

    # right sleeve front 
    start = 106860
    end = 115120
    stage_start = 106821
    stage_size = 135415 - 106822 + 1
    x_axis = Vector3d(101.669 - 103.97, 18.4399 - 16.5086, 0.0) 
    y_axis = Vector3d(102.786 - 101.071 , 5.74841 - 3.70469, 0.0)

    sim.populate(start, end, stage_start, stage_size, scale, x_axis, y_axis)

    # right sleeve back 
    start = 223656
    end = 240784
    stage_start = 214793
    stage_size = 243072 - 214794 + 1

    sim.populate(start, end, stage_start, stage_size, scale, x_axis, y_axis)

    if(solve_static):
        sim.dt = 0.01
        sim.smock_cons = 1.0
        sim.k_stitch = 1e7
        sim.frame_num = 1
        sim.staticSolve = True
        sim.PNTol = 2.5e-4
        sim.withCollision = False
    # breakpoint()
    # populate before sim init
    sim.initialize(sim.cloth_density_iso[clothI], sim.cloth_Ebase_iso[clothI] * membEMult,
    sim.cloth_nubase_iso[clothI], sim.cloth_thickness_iso[clothI], 0, sim.smock)
    sim.bendingStiffMult = bendEMult / membEMult
    sim.kappa_s = Vector2d(0, 0)
    
    # offset_vec = -0.01
    # sim.Offset_stitching(offset_vec, 72 * 2 + 126 * 4)

    sim.load_frame("input/dress/shell_start_4.obj") 

    sim.initialize_OIPC(1e-3, 0)

    sim.run()
