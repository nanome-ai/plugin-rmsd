import nanome
import sys
import time
import math
from .rmsd_calculation import *
# from rmsd_menu import RMSDMenu
from .rmsd_menu import RMSDMenu
from . import rmsd_helpers as help
from nanome.util import Logs
# from .quaternion import Quaternion
import numpy as np
from . import rmsd_selection as selection
from itertools import combinations

class RMSD(nanome.PluginInstance):
    def start(self):
        Logs.debug("Start RMSD Plugin")
        self.args = RMSD.Args()
        self._menu = RMSDMenu(self)
        self._menu.build_menu()
        self.selected_before = []
        self._mobile = []
        self._target = None
        # passed from rmsd_menu.py to compare the index 
        # and autoselect entry menu complex
        self.compare_index = None


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
        self._menu._selected_mobile = []
        self._menu._selected_target = None
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
        self.request_complexes([self._target.index] + [x.index for x in self._mobile], self.on_complexes_received)

    def update_mobile(self, mobile):
        self._mobile = mobile

    def update_target(self, target):
        self._target = target

    def select_all_atoms(self, complex_index):
        complex_list = list(self._mobile)
        complex_list.append(self._target)
        complex_index_list = [x.index for x in complex_list if x != None]
        self.compare_index = complex_index
        self.request_complexes(complex_index_list,self.on_select_atoms_received)

    def on_select_atoms_received(self, complexes):
        for x in complexes:
            if x != None and x.index == self.compare_index:
                for y in x.atoms:
                    y.selected = True
                self.update_structures_deep([x])


    # def on_workspace_received(self, workspace):
    def on_complexes_received(self,complexes):
        target_complex = complexes[0]
        mobile_complex = complexes[1:]
        result = 0

        total_percentage = len(self._mobile)
        percentage_count = 0
        for x in mobile_complex:
            result += self.align(target_complex, x)
            percentage_count += 1
            self._menu.change_loading_percentage(percentage_count/(total_percentage*2))
        if result :
            self._menu.update_score(result)
        self.update_mobile(mobile_complex)
        self.update_target(target_complex)
        Logs.debug("RMSD done")
        self._menu.lock_image()
        self._menu.hide_loading_bar()
        self.update_structures_deep([target_complex])
        self.update_structures_deep(mobile_complex)
        self._menu.loadingBar.percentage = 0
       
    def update_args(self, arg, option):
        setattr(self.args, arg, option)

    class Args(object):
        def __init__(self):
            self.rotation = "kabsch" #alt: "quaternion", "none"
            self.reorder = False
            self.reorder_method = "hungarian" #alt "brute", "distance"
            self.select = "global"
            self.use_reflections = False # scan through reflections in planes (eg Y transformed to -Y -> X, -Y, Z) and axis changes, (eg X and Z coords exchanged -> Z, Y, X). This will affect stereo-chemistry.
            self.use_reflections_keep_stereo = False # scan through reflections in planes (eg Y transformed to -Y -> X, -Y, Z) and axis changes, (eg X and Z coords exchanged -> Z, Y, X). Stereo-chemistry will be kept.
            #exclusion options
            self.no_hydrogen = True
            self.selected_only = True
            self.backbone_only = True
            self.align = True
            self.align_box = False
            self.align_sequence = True

        @property
        def update(self):
            return self.align

        def __str__(self):
            ln = "\n"
            tab = "\t"
            output  = "args:" + ln
            output += tab + "rotation:" + str(self.rotation) + ln
            output += tab + "reorder:" + str(self.reorder) + ln
            output += tab + "reorder_method:" + str(self.reorder_method) + ln
            output += tab + "use_reflections:" + str(self.use_reflections) + ln
            output += tab + "use_reflections_keep_stereo:" + str(self.use_reflections_keep_stereo) + ln
            output += tab + "no_hydrogen:" + str(self.no_hydrogen) + ln
            output += tab + "selected_only:" + str(self.selected_only) + ln
            output += tab + "backbone_only:" + str(self.backbone_only) + ln
            output += tab + "align:" + str(self.align) + ln
            output += tab + "align box:" +str(self.align_box) + ln
            return output

    def align(self, p_complex, q_complex):

        #p is fixed q is mobile
        args = self.args
        p_atoms = list(p_complex.atoms)
        q_atoms = list(q_complex.atoms)

        if args.selected_only:
            p_atoms = help.strip_non_selected(p_atoms)
            q_atoms = help.strip_non_selected(q_atoms)

        if args.no_hydrogen:
            p_atoms = help.strip_hydrogens(p_atoms)
            q_atoms = help.strip_hydrogens(q_atoms)

        if args.backbone_only:
            p_atoms = help.strip_non_backbone(p_atoms)
            q_atoms = help.strip_non_backbone(q_atoms)

        p_size = len(p_atoms)
        q_size = len(q_atoms)

        p_atom_names = get_atom_types(p_atoms)
        q_atom_names = get_atom_types(q_atoms)
        p_pos_orig = help.get_positions(p_atoms)
        q_pos_orig = help.get_positions(q_atoms)
        q_atoms = np.asarray(q_atoms)
        if p_size == 0 or q_size == 0:
            Logs.debug("error: sizes of selected complexes are 0")
            self._menu.change_error("zero_size")
            return False
        if not p_size == q_size:
            Logs.debug("error: Structures not same size receptor size:", q_size, "target size:", p_size)
            self._menu.change_error("different_size")
            return False
        if np.count_nonzero(p_atom_names != q_atom_names) and not args.reorder:
            #message should be sent to nanome as notification?
            msg = "\nerror: Atoms are not in the same order. \n reorder to align the atoms (can be expensive for large structures)."
            Logs.debug(msg)
            self._menu.change_error("different_order")
            return False
        else:
            if(self._menu.error_message.text_value!="Loading..."):
                self._menu.change_error("clear")

        p_coords = help.positions_to_array(p_pos_orig)
        q_coords = help.positions_to_array(q_pos_orig)

        # Create the centroid of P and Q which is the geometric center of a
        # N-dimensional region and translate P and Q onto that center. 
        # http://en.wikipedia.org/wiki/Centroid
        p_cent = centroid(p_coords)
        q_cent = centroid(q_coords)
      
        p_coords -= p_cent
        q_coords -= q_cent

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
                p_coords,
                q_coords,
                reorder_method=reorder_method,
                rotation_method=rotation_method,
                keep_stereo=args.use_reflections_keep_stereo)

        elif args.reorder:
            q_review = reorder_method(p_atom_names, q_atom_names, p_coords, q_coords)
            q_coords = q_coords[q_review]
            q_atom_names = q_atom_names[q_review]
            q_atoms = q_atoms[q_review]
            if not all(p_atom_names == q_atom_names):
                Logs.debug("error: Structure not aligned")
                return False

        #calculate RMSD
        if result_rmsd:
            pass
        elif rotation_method is None:
            result_rmsd = rmsd(p_coords, q_coords)
        else:
            result_rmsd = rotation_method(p_coords, q_coords)
        Logs.debug("result: {0}".format(result_rmsd))

        # Logs.debug result
        if args.update:
            #resetting coords
            p_coords = help.positions_to_array(p_pos_orig)
            q_coords = help.positions_to_array(q_pos_orig)

            p_coords -= p_cent
            q_coords -= q_cent

            #reordering coords  ?
            if args.reorder:
                if q_review.shape[0] != len(q_coords):
                    Logs.debug("error: Reorder length error. Full atom list needed for --Logs.debug")
                    return False
                q_coords = q_coords[q_review]
                q_atoms = q_atoms[q_review]

            # Get rotation matrix
            U = kabsch(p_coords, q_coords)

            #update rotation
            U_matrix = nanome.util.Matrix(4,4)
            for i in range(3):
                for k in range(3):
                    U_matrix[i][k] = U[i][k]
            U_matrix[3][3] = 1
            rot_quat = p_complex.rotation
            rot_matrix = nanome.util.Matrix.from_quaternion(rot_quat)

            result_matrix = rot_matrix * U_matrix
            result_quat = nanome.util.Quaternion.from_matrix(result_matrix)
            q_complex.rotation = result_quat
            Logs.debug("Finished update")
 
            #align centroids
            p_cent = p_complex.rotation.rotate_vector(help.array_to_position(p_cent))
            q_cent = q_complex.rotation.rotate_vector(help.array_to_position(q_cent))

            q_complex.position = p_complex.position + p_cent - q_cent

            if self.args.align_box:
                # maybe it will work here?

                # save the aligned global coord
                matrix1 = q_complex.get_complex_to_workspace_matrix()
                global_pos = map(lambda atom: matrix1 * atom.position, q_complex.atoms)

                # set q rotation to p rotation
                q_complex.rotation = p_complex.rotation

                # restore aligned atom positions from global
                matrix2 = q_complex.get_workspace_to_complex_matrix()
                for (atom, gPos) in zip(q_complex.atoms, global_pos):
                    atom.position = matrix2 * gPos 

            if(self._menu.error_message.text_value=="Loading..."):
                self._menu.change_error("clear")
            
            p_complex.locked = True 
            q_complex.locked = True
        
        return result_rmsd

    # auto select with global/local alignment
    def select(self,mobile,target):
        self._mobile = mobile
        self._target = target
        self.request_workspace(self.on_select_received) 

    def on_select_received(self, workspace):
        complexes = workspace.complexes
        mobile_index_list = list(map(lambda a: a.index, self._mobile))
        self._mobile = []
        for complex in complexes:
            if complex.index in mobile_index_list:
                self._mobile.append(complex)
            if complex.index == self._target.index:
                self._target = complex
       

        if (self.args.select.lower() == "global" and self.args.align_sequence):
            # self.selected_before = [[list(map(lambda a:a.selected,x.atoms)) for x in self._mobile],
            #                         list(map(lambda a:a.selected,self._target.atoms))]
            # 1. DUMMY METHOD
            total_percentage = len(self._mobile) * 2 - 1
            percentage_count = 0
            for x in self._mobile:
                selection.global_align(x , self._target)  
                percentage_count += 1/(total_percentage*2)
                self._menu.change_loading_percentage(percentage_count)

            for x in self._mobile:
                selection.global_align(x , self._target) 
                percentage_count += 1/(total_percentage*2)
                self._menu.change_loading_percentage(percentage_count )

        elif (self.args.select.lower() == "local" and self.args.align_sequence):
            total_percentage = len(self._mobile) * 2 - 1
            percentage_count = 0
            for x in self._mobile:
                selection.local_align(x , self._target)  
                percentage_count += 1/(total_percentage*2)
                self._menu.change_loading_percentage(percentage_count)

            for x in self._mobile:
                selection.local_align(x , self._target) 
                percentage_count += 1/(total_percentage*2)
                self._menu.change_loading_percentage(percentage_count )

        self.workspace = workspace
        self.update_workspace(workspace)
        if self._menu.check_resolve_error():
            self._menu._run_rmsd()

        self._menu.hide_loading_bar() 
        self._menu.loadingBar.percentage = 0
        # self._menu.change_error("clear")
        
    
    # find the two sequences whose distance is the smallest
    # called in on_select_received, used in the clustalW part
    def min_dist(self, matrix):
        if len(matrix) < 2 or len(matrix[0]) < 2:
            Logs.debug("distance matrix size is too small")
            return -1,-1
        else:
            min_val = math.inf
            for x in range(len(matrix)):
                for y in range(len(matrix[0])):
                    if x !=y and matrix[x][y] < min_val:
                        min_val = matrix[x][y]
                        seq1 = x
                        seq2 = y

            return seq1, seq2
 


def main():
    plugin = nanome.Plugin("RMSD", "Aligns complexes using RMSD calculations.", "Alignment", False)
    plugin.set_plugin_class(RMSD)
    plugin.run('127.0.0.1', 8888)

if __name__ == "__main__":
    main()
