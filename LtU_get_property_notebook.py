import h5py
import numpy as np
import arepo_package_notebook as arepo_package
import os


def get_particle_property_LTU(basePath,desired_property,p_type, desired_redshift):

    '''
    LtU version of getting particle properties
    '''
    output_redshift,output_snapshot=arepo_package.desired_redshift_to_output_redshift(basePath,desired_redshift,list_all=False)

    if (output_snapshot < 10):
        snap_shot = '00' + str(output_snapshot)
    else:
        snap_shot = '0' + str(output_snapshot)
        
    file_list = os.listdir(basePath+f'snapdir_{snap_shot}/')
    # print(file_list) # Trying to see the order in which data is loaded
    i=0
    # print('Current snapshot: ',output_snapshot)
    for ii in range(len(file_list)):
        # print(file)
        try:
            h=h5py.File(basePath+f'snapdir_{snap_shot}/'+f'snap_{snap_shot}.{ii}.hdf5')
            #print(file,list(h.keys()))
            PartType5=h.get('PartType%d'%p_type)
            property_list_file=PartType5.get(desired_property)[:]
            if (i==0):
                property_list = property_list_file
            else:
                property_list=np.append(property_list,property_list_file,axis=0)
            i+=1
        except Exception as e:
            continue
            # print(e)
            
    try:
        property_list_file
        
        return property_list,output_redshift
        
    except NameError:
        
        print(f"Failed to load {desired_property} at redshift {desired_redshift}!")

        return np.array(np.nan), output_redshift
                
