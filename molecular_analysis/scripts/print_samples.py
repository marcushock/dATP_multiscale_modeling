import pyemma 
import os
import datetime 
import mdtraj as md
import pickle 

adp_samples_file = "/home/marcus/Documents/myosin_dynamics/1qvi/adp_sampled_frame_data.pkl"

f = open(adp_samples_file, 'rb')
adp_samples = pickle.load(f)
f.close()
print("# Loading samples from {}".format(adp_samples_file))
print("# Modified on:")
modified = os.path.getctime(adp_samples_file)
print("# {}".format(datetime.datetime.fromtimestamp(modified)))

traj_list = []

# for k in range(3):
#     print(k)
#     for traj_n, frame_n in adp_samples[k]:
#         #f traj_n == k:
#         print("{}-{}".format(traj_n,frame_n), end = ' ')
#         traj_list.append(md.load_frame('/media/marcus/OS/Users/marcus/Documents/myosin_dynamics/trajectories/1qvi_adp_{}_2us_10ps.nc'.format(traj_n+1),
#              top = '/media/marcus/OS/Users/marcus/Documents/myosin_dynamics/trajectories/1qvi_adp_strip.prmtop',
#              index = frame_n))
#     print()

file_path = '/media/marcus/OS/Users/marcus/Documents/myosin_dynamics/trajectories'
top_file = '1qvi_adp_strip.prmtop'
for k in range(3):
    print("# Metatsable State {}".format(k))
    print('parm {}/{}'.format(file_path, top_file))
    for traj_n, frame_n in adp_samples[k]:
        #f traj_n == k:
        print("trajin {}/1qvi_adp_{}_2us_10ps.nc {} {}".format(file_path, traj_n+1,frame_n,frame_n))
    
    print('align :1-760@CA')
    print('strip :811-813 parmout 1qvi_protein_test.parm7 ')
    print()
    #print('trajout test_metastatble_{}.pdb\n'.format(k))
    # trajout stripped_cluster.pdb pdb vdw nobox noter
    print('for i=1;i<=15;i++\n\ttrajout metastable_{}_$i.pdb pdb vdw nobox noter start $i stop $i \ndone\n'.format(k))
    print('run')
    print()
    