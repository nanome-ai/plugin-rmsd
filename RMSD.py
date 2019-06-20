import nanome
import sys
import time
from rmsd_calculation import *
from rmsd_menu import RMSDMenu
import rmsd_helpers as help
from nanome.util import Logs

class RMSD(nanome.PluginInstance):
    def start(self):
        Logs.debug("Start RMSD Plugin")
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
        self.request_refresh()

    def request_refresh(self):
        self.request_complex_list(self.on_complex_list_received)
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

            self.align = True

        @property
        def update(self):
            if self.align and (self.no_hydrogen or self.selected_only):
                Logs.debug('Invalid options: cannot use align with "no hydrogen" or "selected only" options')
                return False
            return self.align

    def align(self, complex0, complex1):
        args = RMSD.Args()

        p_atoms = list(complex0.atoms)
        q_atoms = list(complex1.atoms)

        p_size = len(p_atoms)
        q_size = len(q_atoms)

        if not p_size == q_size:
            Logs.debug("error: Structures not same size")
            return

        if not help.same_order(p_atoms, q_atoms) and not args.reorder:
            #message should be sent to nanome as notification?
            msg = "\nerror: Atoms are not in the same order. \n reorder to align the atoms (can be expensive for large structures)."
            Logs.debug(msg)
            return

        if args.selected_only:
            p_atoms = help.strip_nonselected(p_atoms)
            q_atoms = help.strip_nonselected(q_atoms)

        if args.no_hydrogen:
            p_atoms = help.strip_hydrogens(p_atoms)
            q_atoms = help.strip_hydrogens(q_atoms)

        p_atoms, p_coord = get_coordinates(p_atoms)
        q_atoms, q_coord = get_coordinates(q_atoms)

        p_atom_names = list(map(lambda a: a.name, p_atoms))
        q_atom_names = list(map(lambda a: a.name, q_atoms))

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
            return


        # set reorder method
        if not args.reorder:
            reorder_method = None

        if args.reorder_method == "hungarian":
            reorder_method = reorder_hungarian

        elif args.reorder_method == "brute":
            reorder_method = reorder_brute

        elif args.reorder_method == "distance":
            reorder_method = reorder_distance

        else:
            Logs.debug("error: Unknown reorder method:", args.reorder_method)
            return


        # Save the resulting RMSD
        result_rmsd = None


        if args.use_reflections:

            result_rmsd, q_swap, q_reflection, q_review = check_reflections(
                p_atom_names,
                q_atom_names,
                p_coord,
                q_coord,
                reorder_method=reorder_method,
                rotation_method=rotation_method)

        elif args.use_reflections_keep_stereo:

            result_rmsd, q_swap, q_reflection, q_review = check_reflections(
                p_atom_names,
                q_atom_names,
                p_coord,
                q_coord,
                reorder_method=reorder_method,
                rotation_method=rotation_method,
                keep_stereo=True)

        elif args.reorder:

            q_review = reorder_method(p_atom_names, q_atom_names, p_coord, q_coord)
            q_coord = q_coord[q_review]
            q_atom_names = q_atom_names[q_review]

            if not all(p_atom_names == q_atom_names):
                Logs.debug("error: Structure not aligned")
                return


        #calculare RMSD
        if result_rmsd:
            pass

        elif rotation_method is None:
            result_rmsd = rmsd(p_coord, q_coord)

        else:
            result_rmsd = rotation_method(p_coord, q_coord)
        print(type(result_rmsd))
        Logs.debug("{0}".format(result_rmsd))

        # Logs.debug result
        if args.update:

            if args.reorder:

                if q_review.shape[0] != len(q_coord):
                    Logs.debug("error: Reorder length error. Full atom list needed for --Logs.debug")
                    return

                q_coord = q_coord[q_review]
                q_atoms = q_atoms[q_review]

            # Get rotation matrix
            U = kabsch(q_coord, p_coord)

            # recenter all atoms and rotate all atoms
            q_coord -= q_cent
            q_coord = np.dot(q_coord, U)

            # center q on p's original coordinates
            q_coord += p_cent

            # done and done
            for coord, atom in zip(q_coord, q_atoms):
                atom.position.x = coord[0]
                atom.position.y = coord[1]
                atom.position.z = coord[2]
            Logs.debug("Finished update")
        return

        
if __name__ == "__main__":
    # Creates the server, register SimpleHBond as the class to instantiate, and start listening
    plugin = nanome.Plugin("RMSD", "A simple plugin that aligns complexes through RMSD calculation", "Test", False)
    plugin.set_plugin_class(RMSD)
    plugin.run('127.0.0.1', 8888)
