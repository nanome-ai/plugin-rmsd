import numpy as np
import nanome
from nanome.util import Logs

def same_order(atoms1, atoms2):
    for index, _ in enumerate(atoms1):
        if atoms1[index].symbol != atoms2[index].symbol:
            return False
    return True

def get_positions(atom_list):
    return list(map(lambda atom: atom.position, atom_list))

def positions_to_array(pos_list):
    return np.asarray(list(map(lambda pos: position_to_array(pos), pos_list)))

def position_to_array(position):
    return np.asarray([position.x, position.y, position.z], dtype = float)

def array_to_position(array):
    return nanome.util.Vector3(array[0],array[1],array[2])

def strip_hydrogens(atoms):
    return list(filter(lambda a: a.symbol != "H", atoms))

def strip_non_selected(atoms):
    return list(filter(lambda a: a.selected, atoms))

def strip_non_backbone(atoms):
    return list(filter(lambda a: IsBackbone(a), atoms))

def IsBackbone(atom):
    atomName = atom.name
    return atomName == "N" or atomName == "CA" or atomName == "C" or atomName == "O" or atomName == "OXT" or atomName == "OC"

def strip_alternatives(atoms):
    rt = list(filter(lambda a: a._occupancy >= 0.5, atoms))
    oc_test = list(filter(lambda a: a._occupancy <1, rt))     
 
    pairs=set()
    remove_idx=[]
    for i,x in enumerate(rt):
        if x ._occupancy == 0.5:
            if x.name not in pairs:
                pairs.add(x.name)
            else:
                pairs.remove(x.name)
                remove_idx.append(i)
                #rt.remove(x)
    rt = [i for j, i in enumerate(rt) if j not in remove_idx]

    Logs.debug(pairs)
    oc_test = list(filter(lambda a: a._occupancy < 1, rt))     
    Logs.debug(list(map(lambda a: a.name, oc_test)))
    return rt

    
