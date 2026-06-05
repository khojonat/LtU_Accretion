from LtU_get_property import get_particle_property_LTU
import numpy as np
import illustris_python.illustris_python as il
import arepo_package_notebook as arepo_package

def load_fields(path,parttype,fields,redshifts):
    ''' Load in all specified fields for a given parttype over a specified redshift range in code units '''

    Data = {'redshifts':redshifts}
    snap_list,redshift_list=arepo_package.get_snapshot_redshift_correspondence(path)
    
    # header=il.groupcat.loadHeader(path,snap_list[np.where(redshift_list==redshifts[0])])
    # h = header.get('HubbleParam')
    
    for i,z in enumerate(redshifts):
        for field in fields:
            if field not in Data:
                Data[field] = [get_particle_property_LTU(path,field,p_type=parttype, desired_redshift = z)[0]]
            else:
                Data[field].append(get_particle_property_LTU(path,field,p_type=parttype, desired_redshift = z)[0])
    
    return Data
    