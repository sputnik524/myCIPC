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

    # determine the smock pattern type
    smock_sizes = [130,48,64,32]
    fine_mesh_res_ = [106,130,130,190]
    smock_names = ['','braid','twist', 'arrow']
    smock_pattern_type = 1 # 0 for box, 1 for braid, 2 for twist, 3 for arrow
    sim.smock_size = smock_sizes[smock_pattern_type]
    sim.fine_mesh_res = fine_mesh_res_[smock_pattern_type]
    smock_name = smock_names[smock_pattern_type]

    # make the intial smocking
    # sim.add_shell_with_scale_3D_smock("input/"+smock_name +"/S3.obj","input/"+smock_name +"/S1.obj", Vector3d(0, 0, -0.3), Vector3d(0.045, 0.045, 0.045),\
    #     Vector3d(0, 0, 0), Vector3d(0, 0, 1), -90)
    sim.add_shell_with_scale_3D("input/"+smock_name +"/S3.obj", Vector3d(0, 0.3, 0.4), Vector3d(0.02, 0.02, 0.02),\
        Vector3d(0, 0, 0), Vector3d(1, 0, 0), -90)

    # Rescale S2!
    # Offset S3! 

    sim.muComp = StdVectorXd([0, 0, sim.mu,  0, 0, sim.mu,  sim.mu, sim.mu, 0.1])

    sim.dt = 0.04
    sim.frame_dt = 0.04
    sim.frame_num = 100
    sim.withCollision = True
    sim.smock = True
    sim.smock_cons = 0.1
    sim.uniform_stitching_ratio = 0.5
    sim.if_contact = False
    sim.gravity = Vector3d(0, 0, 0) 
    # sim.staticSolve = True
    sim.use_s2 = True

    sim.initialize(sim.cloth_density_iso[clothI], sim.cloth_Ebase_iso[clothI] * membEMult,
        sim.cloth_nubase_iso[clothI], sim.cloth_thickness_iso[clothI], 0, smock = sim.smock, filepath_smock = "input/"+smock_name +"/S2_rescaled.obj", filepath_smock_pattern = "input/"+smock_name +"/S1.obj")
    sim.bendingStiffMult = bendEMult / membEMult
    sim.kappa_s = Vector2d(0, 0)

    # sim.load_frame("input/offset.obj") 
    
    sim.initialize_OIPC(1e-3, 0)

    # sim.run()