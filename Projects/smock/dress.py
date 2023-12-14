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
    smock_sizes = [130,48,64,72]
    fine_mesh_res_ = [106,130,130,130]
    smock_names = ['','braid','twist', 'arrow']
    smock_pattern_type = 3 # 0 for box, 1 for braid, 2 for twist, 3 for arrow
    sim.smock_size = smock_sizes[smock_pattern_type]
    sim.fine_mesh_res = fine_mesh_res_[smock_pattern_type]
    smock_name = smock_names[smock_pattern_type]

    ################################ Not using NH for not PSD in this case!!####################  
    
    # for web meshes, rot x-axis
    sim.add_shell_with_scale_3D_smock("input/"+smock_name +"/S4.obj","input/"+smock_name +"/S1.obj", Vector3d(0.5, 0.70, -0.1), Vector3d(0.065, 0.065, 0.065),\
        Vector3d(0, 0, 0), Vector3d(1, 0, 0), -90) 

    # 2. add the garment with identical size that to be sewed 
    sim.add_shell_with_scale_3D("input/"+smock_name +"/S3.obj", Vector3d(0.1, -0.2, 0.5), #Vector3d(0, 0.2, -0.85)
        Vector3d(0.055, 0.09, 0.055),Vector3d(0, 0, 0), Vector3d(1, 0, 0), 0) 
    
    sim.add_shell_with_scale_3D("input/"+smock_name +"/S3.obj", Vector3d(1.0, 0.2, -0.1), \
        Vector3d(0.07, 0.07, 0.07),Vector3d(0, 0, 0), Vector3d(1, 0, 0), 0)
    
    sim.add_shell_with_scale_3D("input/"+smock_name +"/S3.obj", Vector3d(1.0, -0.2, 0.5), \
        Vector3d(0.07, 0.09, 0.07),Vector3d(0, 0, 0), Vector3d(1, 0, 0), 0)

    # 3. add mannequin
    meshCounter = sim.add_shell_with_scale_3D("input/woman_rot.obj", Vector3d(0.75, 0.5, 0.12), \
        Vector3d(2.7, 2.7, 2.7),Vector3d(0, 0, 0), Vector3d(1, 0, 0), 180) 
    
    # print("The body DBC range:")
    sim.set_DBC_with_range(Vector3d(-0.1, -0.1, -0.1), Vector3d(1.1, 1.1, 1.1), 
            Vector3d(0, 0.0, 0), Vector3d(0, 0, 0), Vector3d(0, 1, 0), 0, meshCounter)

    
    # 4. add sewing info
    sim.add_stitching_withbody()

    sim.muComp = StdVectorXd([0, 0, sim.mu,  0, 0, sim.mu,  sim.mu, sim.mu, 0.1])

    sim.dt = 0.005
    sim.frame_dt = 0.005
    sim.frame_num = 100
    sim.withCollision = True
    sim.smock = True
    sim.smock_cons = 0.1
    sim.uniform_stitching_ratio = 0.1
    sim.if_contact = False
    sim.gravity = Vector3d(9.81, 0, 0) # cancel this term after setting DBC
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
            sim.cloth_nubase_iso[clothI], sim.cloth_thickness_iso[clothI], 0, smock = sim.smock, filepath_smock = "input/"+smock_name +"/S2.obj", filepath_smock_pattern = "input/"+smock_name +"/S1.obj")
        sim.bendingStiffMult = bendEMult / membEMult
        sim.kappa_s = Vector2d(0, 0)
    
    sim.initialize_OIPC(1e-3, 0)

    sim.run()
