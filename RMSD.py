import nanome
import sys
import time
from rmsd_calculation import *
# from rmsd_menu import RMSDMenu
from  rmsd_new_menu import RMSDMenu
import rmsd_helpers as help
from nanome.util import Logs

class RMSD(nanome.PluginInstance):
    def start(self):
        Logs.debug("Start RMSD Plugin")
        self.args = RMSD.Args()
        self._menu = RMSDMenu(self)
        self._menu.build_menu()

    def on_run(self):
        menu = self.menu
        menu.enabled = True
        self._menu._request_refresh()

    def on_complex_added(self):
        nanome.util.Logs.debug("Complex added: refreshing")
        self.request_refresh()

    def on_complex_removed(self):
        nanome.util.Logs.debug("Complex removed: refreshing")
        # self.request_refresh()
        self.request_refresh2()

    def request_refresh(self):
        self.request_complex_list(self.on_complex_list_received)
        nanome.util.Logs.debug("Complex list requested")

    # 
    def request_refresh2(self):
        self.request_complex_list(self.on_complex_list_received2)
        nanome.util.Logs.debug("Complex list requested")


    def update_button(self, button):
        self.update_content(button)

    def make_plugin_usable(self):
        self._menu.make_plugin_usable()

    def on_complex_list_received(self, complexes):
        Logs.debug("complex received: ", complexes)
        self._menu.change_complex_list(complexes)

    def run_rmsd(self, mobile, target):
        self._mobile = mobile
        self._target = target
        self.request_workspace(self.on_workspace_received)

    def on_workspace_received(self, workspace):
        complexes = workspace.complexes
        for complex in complexes:
            if complex.index == self._mobile.index:
                mobile_complex = complex
            if complex.index == self._target.index:
                target_complex = complex

        result = self.align(target_complex, mobile_complex)
        if result :
            self.update_workspace(workspace)
        Logs.debug("RMSD done")
        self.make_plugin_usable()
        self.request_refresh()
    
    def update_args(self, arg, option):
        setattr(self.args, arg, option)

    class Args(object):
        def __init__(self):
            self.rotation = "kabsch" #alt: "quaternion", "none"
            self.reorder = False
            self.reorder_method = "hungarian" #alt "brute", "distance"
            self.use_reflections = False # scan through reflections in planes (eg Y transformed to -Y -> X, -Y, Z) and axis changes, (eg X and Z coords exchanged -> Z, Y, X). This will affect stereo-chemistry.
            self.use_reflections_keep_stereo = False # scan through reflections in planes (eg Y transformed to -Y -> X, -Y, Z) and axis changes, (eg X and Z coords exchanged -> Z, Y, X). Stereo-chemistry will be kept.
            #exclusion options
            self.no_hydrogen = False
            self.selected_only = False
            self.backbone_only = False
            self.align = True

        @property
        def update(self):
            if self.align and (self.no_hydrogen or self.selected_only):
                Logs.debug('Invalid options: cannot use align with "no hydrogen" or "selected only" options')
                return False
            return self.align

    def align(self, complex0, complex1):
        args = self.args
        p_atoms = list(complex0.atoms)
        q_atoms = list(complex1.atoms)

        p_size = len(p_atoms)
        q_size = len(q_atoms)

        if not p_size == q_size:
            Logs.debug("error: Structures not same size")
            return False

        if args.selected_only:
            p_atoms = help.strip_nonselected(p_atoms)
            q_atoms = help.strip_nonselected(q_atoms)

        if args.no_hydrogen:
            p_atoms = help.strip_hydrogens(p_atoms)
            q_atoms = help.strip_hydrogens(q_atoms)

        if args.backbone_only:
            p_atoms = help.strip_non_backbone(p_atoms)
            q_atoms = help.strip_non_backbone(q_atoms)

        p_atom_names, p_coord_orig = get_coordinates(p_atoms)
        q_atom_names, q_coord_orig = get_coordinates(q_atoms)
        q_atoms = np.asarray(q_atoms)
        if np.count_nonzero(p_atom_names != q_atom_names) and not args.reorder:
            #message should be sent to nanome as notification?
            msg = "\nerror: Atoms are not in the same order. \n reorder to align the atoms (can be expensive for large structures)."
            Logs.debug(msg)
            return False

        p_coord = copy.deepcopy(p_coord_orig)
        q_coord = copy.deepcopy(q_coord_orig)

        # Create the centroid of P and Q which is the geometric center of a
        # N-dimensional region and translate P and Q onto that center. 
        # http://en.wikipedia.org/wiki/Centroid
        p_cent = centroid(p_coord)
        q_cent = centroid(q_coord)
        p_coord -= p_cent
        q_coord -= q_cent

        # set rotation method
        if args.rotation.lower() == "kabsch":
            rotation_method = kabsch_rmsd
        elif args.rotation.lower() == "quaternion":
            rotation_method = quaternion_rmsd
        elif args.rotation.lower() == "none":
            rotation_method = None
        else:
            Logs.debug("error: Unknown rotation method:", args.rotation)
            return False

        # set reorder method
        # when reorder==False, set reorder_method to "None"
        if not args.reorder:
            reorder_method = None
        elif args.reorder_method.lower() == "hungarian":
            reorder_method = reorder_hungarian
        elif args.reorder_method.lower() == "brute":
            reorder_method = reorder_brute
        elif args.reorder_method.lower() == "distance":
            reorder_method = reorder_distance
        else:
            Logs.debug("error: Unknown reorder method:", args.reorder_method)
            Logs.debug("The value of reorder is: ",args.reorder)
            return False


        # Save the resulting RMSD
        result_rmsd = None

        if args.use_reflections or args.use_reflections_keep_stereo:
            result_rmsd, q_swap, q_reflection, q_review = check_reflections(
                p_atom_names,
                q_atom_names,
                p_coord,
                q_coord,
                reorder_method=reorder_method,
                rotation_method=rotation_method,
                keep_stereo=args.use_reflections_keep_stereo)
        elif args.reorder:
            q_review = reorder_method(p_atom_names, q_atom_names, p_coord, q_coord)
            q_coord = q_coord[q_review]
            q_atom_names = q_atom_names[q_review]
            q_atoms = q_atoms[q_review]
            if not all(p_atom_names == q_atom_names):
                Logs.debug("error: Structure not aligned")
                return False

        #calculate RMSD
        if result_rmsd:
            pass
        elif rotation_method is None:
            result_rmsd = rmsd(p_coord, q_coord)
        else:
            result_rmsd = rotation_method(p_coord, q_coord)
        Logs.debug("result: {0}".format(result_rmsd))
        self._menu.update_score(result_rmsd)

        # Logs.debug result
        if args.update:
            if args.reorder:
                if q_review.shape[0] != len(q_coord_orig):
                    Logs.debug("error: Reorder length error. Full atom list needed for --Logs.debug")
                    return False
                q_coord_orig = q_coord_orig[q_review]
                q_atoms = q_atoms[q_review]
            
            q_complex_position = complex1.position

            # Get rotation matrix
            U = kabsch(q_coord_orig, p_coord)

            # recenter all atoms and rotate all atoms
            q_coord_orig -= q_cent
            q_coord_orig = np.dot(q_coord_orig, U)

            # center q on p's original coordinates
            q_coord_orig += p_cent

            print("num q_atoms", (q_atoms.shape[0]))
            print("q_atoms.type", type(q_atoms[0]))
            # done and done
            print(q_atoms[0].position.x, q_coord_orig[0][0])
            for coord, atom in zip(q_coord_orig, q_atoms):
                atom.position.x = coord[0]
                atom.position.y = coord[1]
                atom.position.z = coord[2]
            complex1.name = complex1.name[::-1]
            #complex1.position = complex0.position
            #complex1.rotation = complex0.rotation
            Logs.debug("Finished update")
        return result_rmsd

if __name__ == "__main__":
    # Creates the server, register SimpleHBond as the class to instantiate, and start listening
    plugin = nanome.Plugin("RMSD", "A simple plugin that aligns complexes through RMSD calculation", "Test", False)
    plugin.set_plugin_class(RMSD)
    plugin.run('127.0.0.1', 8888)