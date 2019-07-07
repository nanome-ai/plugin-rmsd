def same_order(atoms1, atoms2):
    for index, _ in enumerate(atoms1):
        if atoms1[index].symbol != atoms2[index].symbol:
            return False
    return True

def strip_hydrogens(atoms):
    return list(filter(lambda a: a.symbol != "H", atoms))

def strip_nonselected(atoms):
    return list(filter(lambda a: a.selected, atoms))

def strip_non_backbone(atoms):
    return list(filter(lambda a: IsBackbone(a), atoms))

def get_coordinates(atoms):
    coords = list()

    for atom in atoms:
        coords.append(atom.position.get_copy())

    return coords

def IsBackbone(atom):
    atomName = atom.name
    return atomName == "N" or atomName == "CA" or atomName == "C" or atomName == "O" or atomName == "OXT" or atomName == "OC"