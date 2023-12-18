xyz_template="""1
%s
%s 0.0 0.0 0.0"""

pack_template="""tolerance 15.0
discale 1.5
filetype xyz
output {}

"""

repeat_template="""structure {}
  number {}
  inside box {} {} {} {} {} {}
  radius {}
end structure

"""


def gen_lattice_packmol(box_length,particle_types, particle_type_list, diameter_list, packmol_path="/scratch/projects/hockygroup/gmh4/codes/packmol/packmol", outprefix=None, box_fudge=0.99):
    import numpy as np
    import tempfile
    import os
    import shutil
    box_size = np.array([box_length,box_length,box_length])
    box_min = -box_size/2.*box_fudge
    box_max = box_size/2.*box_fudge
    unique_type_ids = list(np.unique(particle_type_list))
    unique_types = np.array(particle_types)[unique_type_ids]
    radius_list = np.array(diameter_list)[unique_type_ids]/2.

    run_dir = tempfile.mkdtemp()
    print("running in",run_dir)

    pack_input = os.path.join(run_dir,f"mixture.inp")
    pack_output = os.path.join(run_dir,f"mixture.xyz")

    inp_fh = open(pack_input,'w')
    inp_fh.write(pack_template.format(pack_output))

    for typeidx, type_id in enumerate(unique_type_ids):
        n_object = int(np.sum(particle_type_list==type_id))
        type_name = unique_types[typeidx]
        coordinates = xyz_template%(type_name, type_name)
        struct_file = os.path.join(run_dir,f"{type_name}.xyz")

        fh = open(struct_file,'w')
        fh.write(coordinates)
        fh.close()

        inp_fh.write(repeat_template.format(struct_file,n_object, \
                     box_min[0], box_min[1], box_min[2], \
                     box_max[0], box_max[1], box_max[2], \
                     radius_list[typeidx], \
                    ))
    inp_fh.close()
    os.system(f"{packmol_path} < {pack_input}")
    if outprefix is not None:
        shutil.copyfile(pack_input,outprefix+'.pack.inp')
        shutil.copyfile(pack_output,outprefix+'.pack.xyz')

    initial_cfg = np.loadtxt(pack_output,skiprows=2,usecols=(1,2,3)) 
    shutil.rmtree(run_dir)
    return initial_cfg



