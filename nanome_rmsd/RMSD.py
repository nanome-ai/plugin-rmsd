import nanome
import sys
import time
from .rmsd_calculation import *
# from rmsd_menu import RMSDMenu
from .rmsd_menu import RMSDMenu
from . import rmsd_helpers as help
from nanome.util import Logs
# from .quaternion import Quaternion
import numpy as np
from . import rmsd_selection as selection

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
        self._menu.change_error("loading")
        self._mobile = mobile
        self._target = target
        self.request_workspace(self.on_workspace_received) 


    def on_workspace_received(self, workspace):
        complexes = workspace.complexes
        mobile_index_list = list(map(lambda a: a.index, self._mobile))
        mobile_complex = []
        for complex in complexes:
            if complex.index in mobile_index_list:
                mobile_complex.append(complex)
            if complex.index == self._target.index:
                target_complex = complex
        self.workspace = workspace
        result = 0
        for x in mobile_complex:
            result += self.align(target_complex, x)
        if result :
            self._menu.update_score(result)
            
        
        Logs.debug("RMSD done")
        self.make_plugin_usable()
        self.update_workspace(workspace)
        #self.lock_result()
        
        # self.request_refresh()
    
    def lock_result(self):
        for x in self._mobile:
                x.locked = True
        self._target.locked = True
        # self.update_workspace(workspace)

        self.update_structures_shallow([self._target])
        self.update_structures_shallow(self._mobile)

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

        # p_atoms = help.strip_alternatives(p_atoms)
        # q_atoms = help.strip_alternatives(q_atoms)

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

            #reordering coords  ???
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
           
            if(self._menu.error_message.text_value=="Loading..."):
                self._menu.change_error("clear")
            
            p_complex.locked = True
            q_complex.locked = True
        
        return result_rmsd

    # auto select with global/local alignment
    def select(self,mobile,target):
        # create seq1 and seq2
        #run the alignment
        self._menu.change_error("loading")
        #self._mobile = mobile
        self._mobile = mobile
        # Logs.debug("mobile list before ",self._mobile)
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
        
        if (self.args.select.lower() == "global"):
            self.selected_before = [[list(map(lambda a:a.selected,x.atoms)) for x in self._mobile],
                                    list(map(lambda a:a.selected,self._target.atoms))]
            for x in self._mobile:
                selection.global_align(x , self._target)  
            for x in self._mobile[:-1]:
                selection.global_align(x , self._target)  

        if (self.args.select.lower() == "local"):
            selection.local_align(self._mobile,self._target)
        if (self.args.select.lower() == "none"):
            if self.selected_before:
                self.change_selected(self._mobile,self._target,self.selected_before[0],self.selected_before[1])
        self.workspace = workspace
        self.make_plugin_usable()
        self.update_workspace(workspace)
        self._menu.change_error("clear")
        

    def change_selected(self,mobile,target,mobile_selected,target_selected):
        if [len(list(map(lambda a:a,x.atoms))) for x in mobile]==[len(x) for x in mobile_selected]\
            and len(list(map(lambda a:a,target.atoms))) == len(target_selected):
            for j,y in enumerate(mobile):    
                for i,x in enumerate(y.atoms):
                    x.selected = mobile_selected[j][i]
            for i,x in enumerate(target.atoms):
                x.selected = target_selected[i]
        else:
            Logs.debug("selected complexes changed")
            self._menu.change_error("selected_changed")


def main():
    plugin = nanome.Plugin("RMSD", "A simple plugin that aligns complexes through RMSD calculation", "Test", False)
    plugin.set_plugin_class(RMSD)
    plugin.run('127.0.0.1', 8888)

if __name__ == "__main__":
    main()
