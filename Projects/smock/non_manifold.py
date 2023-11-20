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
    
    # sim.add_shell_with_scale_3D_smock("input/smock_box.obj","input/smock_box_pattern.obj", Vector3d(0, 6.0, 8.0), Vector3d(1.0, 1.0, 1.0),\
    #     Vector3d(0, 0, 0), Vector3d(1, 0, 0), -90)

    sim.add_shell("input/smock_box_rest.obj", Vector3d(0, 0.0, 0.0), \
        Vector3d(0, 0, 0), Vector3d(1, 0, 0), 0)

    sim.muComp = StdVectorXd([0, 0, sim.mu,  0, 0, sim.mu,  sim.mu, sim.mu, 0.1])

    sim.dt = 0.04
    sim.frame_dt = 0.04
    sim.frame_num = 100
    sim.withCollision = True

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
            sim.cloth_nubase_iso[clothI], sim.cloth_thickness_iso[clothI], 0)
        sim.bendingStiffMult = bendEMult / membEMult
        sim.kappa_s = Vector2d(0, 0)
    elif algI == 2:
        # aniso
        sim.initialize(sim.cloth_density[clothI], sim.cloth_Ebase[clothI], # actually only affect bending
            0, sim.cloth_thickness[clothI], 0)
        sim.bendingStiffMult = bendEMult
        sim.fiberStiffMult = sim.cloth_weftWarpMult[clothI] * membEMult
        sim.inextLimit = sim.cloth_inextLimit[clothI]
        sim.kappa_s = Vector2d(0, 0)
    # elif algI == 6:
        # coupled strain limiting
    # elif algI == 7:
        # split strain limiting

    sim.initialize_OIPC(1e-3, 0)

    sim.load_frame_nm("input/smock_box.obj","input/smock_box_pattern.obj")

    sim.run()
