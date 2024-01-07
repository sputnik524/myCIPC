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
    smock_names = ['box','braid','twist','arrow','leaf','braid_2','twist_2']
    smock_pattern_type = 1
    sim.smock_size = 0 # remain manifold
    sim.fine_mesh_res = 130
    smock_name = smock_names[smock_pattern_type]

    ################################ in-situ Smocking!! ####################  
    
    sim.add_shell_with_scale_3D("input/"+smock_name +"/S3.obj", Vector3d(0.3, 0.3, -0.05), Vector3d(0.04, 0.04, 0.04),\
        Vector3d(0, 0, 0), Vector3d(1, 0, 0), 0)

    # 2. add the garment with identical size that to be sewed 
    sim.add_shell_with_scale_3D("input/"+smock_name +"/S3.obj", Vector3d(0.3, 0.1, 0.45), #Vector3d(0, 0.2, -0.85)
        Vector3d(0.04, 0.07, 0.07),Vector3d(0, 0, 0), Vector3d(1, 0, 0), 0) 
    
    sim.add_shell_with_scale_3D("input/"+smock_name +"/S3.obj", Vector3d(1.0, 0.22, -0.05), \
        Vector3d(0.05, 0.05, 0.05),Vector3d(0, 0, 0), Vector3d(1, 0, 0), 0)
    
    sim.add_shell_with_scale_3D("input/"+smock_name +"/S3.obj", Vector3d(1.0, 0.1, 0.45), \
        Vector3d(0.05, 0.07, 0.07),Vector3d(0, 0, 0), Vector3d(1, 0, 0), 0)

    # 3. add mannequin
    meshCounter = sim.add_shell_with_scale_3D("input/woman_rot.obj", Vector3d(0.75, 0.57, 0.12), \
        Vector3d(2.3, 2.3, 2.3),Vector3d(0, 0, 0), Vector3d(1, 0, 0), 180) 
    
    # print("The body DBC range:")
    sim.set_DBC_with_range(Vector3d(-0.1, -0.1, -0.1), Vector3d(1.1, 1.1, 1.1), 
            Vector3d(0, 0.0, 0), Vector3d(0, 0, 0), Vector3d(0, 1, 0), 0, meshCounter)

    sim.muComp = StdVectorXd([0, 0, sim.mu,  0, 0, sim.mu,  sim.mu, sim.mu, 0.1])

    sim.dt = 0.01
    sim.frame_dt = 0.01
    sim.frame_num = 100
    sim.k_stitch = 1000
    sim.withCollision = True
    sim.smock = True
    sim.smock_cons = 1.0
    sim.uniform_stitching_ratio = 0.5
    sim.uniform_stitching_ratio_smock = 1.0
    sim.if_contact = False
    sim.use_s2 = True
    sim.gravity = Vector3d(9.81, 0, 0)

    sim.initialize(sim.cloth_density_iso[clothI], sim.cloth_Ebase_iso[clothI] * membEMult,
        sim.cloth_nubase_iso[clothI], sim.cloth_thickness_iso[clothI], 0, smock = sim.smock, filepath_smock = "input/"+smock_name +"/S2_rescaled_04.obj", filepath_smock_pattern = "input/"+smock_name +"/S1.obj")
    sim.bendingStiffMult = bendEMult / membEMult
    sim.kappa_s = Vector2d(0, 0)

    offset_vec = Vector3d(0, 0, 0.02)
    sim.Offset_smocking(offset_vec) 
    sim.add_stitching_withbody()
    
    sim.initialize_OIPC(1e-3, 0)
    # sim.load_frame("input/"+smock_name +"/shell0.obj")

    # sim.run()
