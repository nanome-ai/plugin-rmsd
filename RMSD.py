import nanome
import sys
import time
from rmsd_calculation import *
from rmsd_menu import RMSDMenu

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
        print(self._mobile.index)
        print(self._target.index)
        for complex in complexes:
            if complex.index == self._mobile.index:
                mobile_complex = complex
            if complex.index == self._target.index:
                target_complex = complex

        mobile_complex = self.align(target_complex, mobile_complex)
        self.update_workspace(workspace)
        Logs.debug("RMSD done")
        self.make_plugin_usable()
 
    def align(self, complex0, complex1):
        # flags
        rotation = "kabsch"
        # mutually exclusive flags
        no_hydrogen = False
        remove_idx = 0
        add_idx = 0
        p_all_atoms, p_all = get_coordinates(complex0)
        q_all_atoms, q_all = get_coordinates(complex1)
        p_common_atoms, p_common, q_common_atoms, q_common = get_common_coordinates(complex0, complex1)

        p_size = p_common.shape[0]
        q_size = q_common.shape[0]

        if not p_size == q_size:
            print("error: Structures not same size")
            print("complex0 size = " + str(p_size))
            print("complex1 size = " + str(q_size))
            quit()
        if p_size == 0:
            print("error: no common residue found")
            quit()

        p_common_symbols = list()
        q_common_symbols = list()
        for atom in p_common_atoms:
            p_common_symbols.append(atom.symbol)
        for atom in q_common_atoms:
            q_common_symbols.append(atom.symbol)

        if np.count_nonzero(p_common_symbols != q_common_symbols) and not args.reorder:
            msg = """
error: Atoms are not in the same order.

Use --reorder to align the atoms (can be expensive for large structures).

Please see --help or documentation for more information or
https://github.com/charnley/rmsd for further examples.
"""
            print(msg)
            exit()

        # Set local view
        p_view = None
        q_view = None

        if no_hydrogen:
            p_view = np.where(p_common_atoms != 'H')
            q_view = np.where(q_common_atoms != 'H')

        elif remove_idx:
            index = range(p_size)
            index = set(index) - set(remove_idx)
            index = list(index)
            p_view = index
            q_view = index

        elif add_idx:
            p_view = add_idx
            q_view = add_idx

        # Set local view
        if p_view is None:
            p_coord = copy.deepcopy(p_common)
            q_coord = copy.deepcopy(q_common)

        else:
            p_coord = copy.deepcopy(p_common[p_view])
            q_coord = copy.deepcopy(q_common[q_view])

        # Create the centroid of P and Q which is the geometric center of a
        # N-dimensional region and translate P and Q onto that center.
        # http://en.wikipedia.org/wiki/Centroid
        p_cent = centroid(p_coord)
        q_cent = centroid(q_coord)
        p_coord -= p_cent
        q_coord -= q_cent

        # set rotation method
        if rotation.lower() == "kabsch":
            rotation_method = kabsch_rmsd

        elif rotation.lower() == "quaternion":
            rotation_method = quaternion_rmsd

        elif rotation.lower() == "none":
            rotation_method = None

        else:
            print("error: Unknown rotation method:", rotation)
            quit()

        # Save the resulting RMSD
        result_rmsd = None

        # Get rotation matrix
        U = kabsch(q_coord, p_coord)

        # recenter all atoms and rotate all atoms
        q_all -= q_cent
        q_all = np.dot(q_all, U)

        # center q on p's original coordinates
        q_all += p_cent

        # done and done
        for i in range(0, int(q_all.shape[0])):
            q_all_atoms[i].position = nanome.util.Vector3(q_all[i][0], q_all[i][1], q_all[i][2])

        # align the complex itself
        Logs.debug(complex1.position)
        Logs.debug(complex0.position)
        complex1.position = complex0.position
        complex1.rotation = complex0.rotation

        return complex1

        
if __name__ == "__main__":
    # Creates the server, register SimpleHBond as the class to instantiate, and start listening
    plugin = nanome.Plugin("RMSD", "A simple plugin that aligns complexes through RMSD calculation", "Test", False)
    plugin.set_plugin_class(RMSD)
    plugin.run('127.0.0.1', 8888)
