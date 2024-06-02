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

    sim.add_shell_with_scale_3D("input/dress/high_cloth_wraping.obj", Vector3d(0.0, 0.0, -0.0), Vector3d(1.0, 1.0, 1.0),\
        Vector3d(0, 0, 0), Vector3d(1, 0, 0), 0)
    

    sim.dt = 0.002
    sim.frame_dt = 0.002
    sim.frame_num = 0

    sim.use_populate = True
    
    # wb front data
    stage_size = 12911
    stage_start = 135415
    start = 145238
    end = 139400
    scale = 1.0

    sim.populate(start, end, stage_start, stage_size, scale)
    
    sim.initialize(sim.cloth_density_iso[clothI], sim.cloth_Ebase_iso[clothI] * membEMult,
    sim.cloth_nubase_iso[clothI], sim.cloth_thickness_iso[clothI], 0)
    sim.bendingStiffMult = bendEMult / membEMult
    sim.kappa_s = Vector2d(0, 0)
    
    offset_vec = -0.15
    sim.Offset_stitching(offset_vec)

    sim.initialize_OIPC(1e-3, 0)

    sim.run()
