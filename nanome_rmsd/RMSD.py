import nanome
import sys
import time
from .rmsd_calculation import *
# from rmsd_menu import RMSDMenu
from .rmsd_menu import RMSDMenu
from . import rmsd_helpers as help
from nanome.util import Logs
# from .quaternion import Quaternion

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
        self._menu._selected_mobile = None
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
        self.request_workspace(self.on_workspace_received)

    def on_workspace_received(self, workspace):
        complexes = workspace.complexes
        for complex in complexes:
            if complex.index == self._mobile.index:
                mobile_complex = complex
            if complex.index == self._target.index:
                target_complex = complex
        self.workspace = workspace
        result = self.align(target_complex, mobile_complex)
        if result :
            self.update_workspace(workspace)
        Logs.debug("RMSD done")
        self.make_plugin_usable()
        # self.request_refresh()
    
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
            self.no_hydrogen = False
            self.selected_only = False
            self.backbone_only = False
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
        self._menu.update_score(result_rmsd)

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
        return result_rmsd

    # auto select with global/local alignment
    def select(self,mobile,target):
        # create seq1 and seq2
        #run the alignment
        self._mobile = mobile
        self._target = target
        if (self.args.select.lower() == "global"):
            Logs.debug("type of mobile is ",type(mobile))
            self.global_align(self._mobile , self._target)
        if (self.args.select.lower() == "local"):
            self.local_align(self._mobile,self._target)

    # needleman wunsch algorithm
    def global_align(self, res1, res2, gap_penalty = -1, mismatch_penalty = -1, match_reward = 2):
        Logs.debug("type in global align function is",type(res1))
        seq1 = res1.residues
        seq2 = res2.residues
        #create the table
        m, n = len(seq1), len(seq2)
        score = zeros((m+1, n+1))      
    
        # file the first column and row
        for i in range(0, m + 1):
            score[i][0] = gap_penalty * i
        for j in range(0, n + 1):
            score[0][j] = gap_penalty * j
        # fill the table wtih scores
        for i in range(1, m + 1):
            for j in range(1, n + 1):
                if seq1[i-1] == seq2[j-1]:
                    match = score[i - 1][j - 1] + match_reward
                else:
                    match = score[i - 1][j - 1] + mismatch_penalty
                delete = score[i - 1][j] + gap_penalty
                insert = score[i][j - 1] + gap_penalty
                score[i][j] = max(match, delete, insert)

        # Traceback and compute the alignment 
        # align1, align2 = '', ''
        # start from the bottom right cell
        i,j = m,n
        while i > 0 and j > 0: 
            score_current = score[i][j]
            score_diagonal = score[i-1][j-1]
            score_up = score[i][j-1]
            score_left = score[i-1][j]
            # two residuses match, only deselect when the selected atoms are not matched
            if score_current == score_diagonal + mismatch_penalty and seq1[i-1] == seq2[j-1]:
                # align1 += seq1[i-1]
                # align2 += seq2[j-1]
                
                i -= 1
                j -= 1
            # two of the residues do not match, the deselect both
            elif score_current == score_diagonal + match_reward and seq1[i-1] != seq2[j-1]:
                for x in res1.atoms:
                    x.select = False
                for y in res2.atoms:
                    y.select = False
                i -= 1
                j -= 1
            # seq1 has an extra residue, deselect it
            elif score_current == score_left + gap_penalty:
                # align1 += seq1[i-1]
                # align2 += '-'
                for x in res1.atoms:
                    x.select = False
                i -= 1
            # seq2 has an extra residue, deselect it
            elif score_current == score_up + gap_penalty:
                # align1 += '-'
                # align2 += seq2[j-1]
                for x in res2.atoms:
                    x.select = False
                j -= 1

        # Finish tracing up to the top left cell
        # while i > 0:
        #     align1 += seq1[i-1]
        #     align2 += '-'
        #     i -= 1
        # while j > 0:
        #     align1 += '-'
        #     align2 += seq2[j-1]
        #     j -= 1

        # finalize(align1, align2)

    def local_align(res1,res2):
        pass
                

def main():
    plugin = nanome.Plugin("RMSD", "A simple plugin that aligns complexes through RMSD calculation", "Test", False)
    plugin.set_plugin_class(RMSD)
    plugin.run('127.0.0.1', 8888)

if __name__ == "__main__":
    main()
