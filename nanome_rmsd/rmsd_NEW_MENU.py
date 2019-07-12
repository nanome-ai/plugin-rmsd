import nanome
from nanome.util import Logs

import os

class RMSDMenu():
    def __init__(self, rmsd_plugin):
        self._menu = rmsd_plugin.menu
        self._plugin = rmsd_plugin
        self._selected_mobile = None # button
        self._selected_target = None # button
        self._run_button = None
        self._current_tab = "receptor" #receptor = 0, target = 1
        self._drop_down_dict={"rotation":["None", "Kabsch","Quaternion"],"reorder_method":["None","Hungarian","Brute", "Distance"]}
        self._current_reorder = "None"
        self._current_rotation = "None"

    def _request_refresh(self):
        self._plugin.request_refresh()

    # run the rmsd algorithm
    def _run_rmsd(self):
        if self._selected_mobile != None or self._selected_target != None:
            self._plugin.run_rmsd(self._selected_mobile.complex, self._selected_target.complex)
        else:
            # TODO, need to let the RUN button able to be selected again
            pass

    # change the args in the plugin
    def update_args(self,arg,option):
        self._plugin.update_args(arg,option)

    def update_score(self,value):
        Logs.debug("update score called: ",value)
        self.rmsd_score_label.text_value = str("%.3g"%value)
        self._plugin.update_menu(self._menu)

    def make_plugin_usable(self, state = True):
        self._run_button.unusable = not state
        self._plugin.update_button(self._run_button)

    # change the complex list
    def change_complex_list(self, complex_list):

        def mobile_pressed(button):
            if self._selected_mobile != None:
                self._selected_mobile.selected = False
            button.selected = True
            self._selected_mobile = button
            self.receptor_text.text_value ="Receptor: "+ button.complex.name
            self._plugin.update_menu(self._menu)

        def target_pressed(button):
            if self._selected_target != None:
                self._selected_target.selected = False
            button.selected = True
            self._selected_target = button
            self.target_text.text_value = "Target: "+button.complex.name
            self._plugin.update_menu(self._menu)

        self._mobile_list = []
        self._target_list = []

        for complex in complex_list:
            clone = self._complex_item_prefab.clone()
            ln_btn = clone.get_children()[0]
            btn = ln_btn.get_content()
            btn.set_all_text(complex.name)
            btn.complex = complex
            btn.register_pressed_callback(mobile_pressed)
            self._mobile_list.append(clone)
            
            #clone1 = clone.clone()
            clone1 = self._complex_item_prefab.clone()
            ln_btn = clone1.get_children()[0]
            btn = ln_btn.get_content()
            btn.set_all_text(complex.name)
            btn.complex = complex
            btn.register_pressed_callback(target_pressed)
            self._target_list.append(clone1)
        if self._selected_mobile == None:
            self.receptor_text.text_value ="Receptor: Unselected"
        if self._selected_target == None:
            self.target_text.text_value ="Target: Unselected "
 
        if self._current_tab == "receptor":
            self._show_list.items=self._mobile_list
        else:
            self._show_list.items=self._target_list

        self._plugin.update_menu(self._menu)

    # build the menu
    def build_menu(self):
        # refresh the lists
        def refresh_button_pressed_callback(button):
            self._request_refresh()

        # press the run button and run the algorithm
        def run_button_pressed_callback(button):
            self.make_plugin_usable(False)
            self._run_rmsd()

        # show the target list when the receptor tab is pressed
        def receptor_tab_pressed_callback(button):
            self._current_tab="receptor"
            receptor_tab.selected = True
            target_tab.selected = False
            self._show_list.items = self._mobile_list
            self._plugin.update_menu(self._menu)

        # show the target list when the target tab is pressed
        def target_tab_pressed_callback(button):
            self._current_tab="target"
            target_tab.selected = True
            receptor_tab.selected = False
            self._show_list.items = self._target_list
            self._plugin.update_menu(self._menu)
            
        
        # no hydrogen = ! no hydrogen
        def no_hydrogen_button_pressed_callback(button):
            self.update_args("no_hydrogen", False)
            no_hydrogen_button.selected = not no_hydrogen_button.selected
            self._plugin.update_menu(self._menu)

        # use reflections = ! use reflections
        # def use_reflections_button_pressed_callback(button):
        #     self.update_args("use_reflections", False)
        #     use_reflections_button.selected = not use_reflections_button.selected
        #     self._plugin.update_menu(self._menu)

        # backbone only = ! backbone only
        def backbone_only_button_pressed_callback(button):
            self.update_args("backbone_only", False)
            backbone_only_button.selected = not backbone_only_button.selected
            self._plugin.update_menu(self._menu)
        
        # selected only = ! selected only
        def selected_only_button_pressed_callback(button):
            self.update_args("selected_only", False)
            selected_only_button.selected = not selected_only_button.selected            
            self._plugin.update_menu(self._menu)

        # change Reorder to the next option
        def reorder_button_pressed_callback(button):
            drop_down  = self._drop_down_dict["reorder_method"]
            temp_length=len(drop_down)
            
            pre_index = drop_down.index(self._current_reorder)
            post_index = (pre_index + 1) % temp_length

            post_option = drop_down[post_index]

            reorder_button.selected = post_option == "None"
            reorder_button.set_all_text(post_option)
            
            # tell the plugin and update the menu
            self._current_reorder = post_option
            self.update_args("reorder_method", post_option)
            self.update_args("reorder", post_option != "None")
            self._plugin.update_menu(self._menu)

        # change Rotation to the next option
        def rotation_button_pressed_callback(button):
            drop_down  = self._drop_down_dict["rotation_method"]
            temp_length=len(drop_down)
            
            pre_index = drop_down.index(self._current_rotation)
            post_index = (pre_index + 1) % temp_length

            post_option = drop_down[post_index]

            rotation_button.selected = post_option == "None"
            rotation_button.set_all_text(post_option)
            
            # tell the plugin and update the menu
            self._current_rotation = post_option
            self.update_args("rotation_method", post_option)
            self._plugin.update_menu(self._menu)
        
        # Create a prefab that will be used to populate the lists
        self._complex_item_prefab = nanome.ui.LayoutNode()
        self._complex_item_prefab.layout_orientation = nanome.ui.LayoutNode.LayoutTypes.horizontal
        child = self._complex_item_prefab.create_child_node()
        child.name = "complex_button"
        prefabButton = child.add_new_button()
        prefabButton.text.active = True

        # import the json file of the new UI
        menu = nanome.ui.Menu.io.from_json(os.path.join(os.path.dirname(__file__), 'rmsd_pluginator.json'))
        self._plugin.menu = menu

        # create the Run button
        self._run_button = menu.root.find_node("Run", True).get_content()
        self._run_button.register_pressed_callback(run_button_pressed_callback)

        # create the Refresh button
        refresh_button = menu.root.find_node("Refresh", True).get_content()
        refresh_button.register_pressed_callback(refresh_button_pressed_callback)

        # create the List 
        self._show_list = menu.root.find_node("List", True).get_content()
        self._mobile_list = []
        self._target_list = []

        # create the Receptor tab
        receptor_tab = menu.root.find_node("Receptor_tab",True).get_content()
        receptor_tab.register_pressed_callback(receptor_tab_pressed_callback)

        # create the Target tab
        target_tab = menu.root.find_node("Target_tab",True).get_content()
        target_tab.register_pressed_callback(target_tab_pressed_callback)

        # create the no hydrogen button
        no_hydrogen_button = menu.root.find_node("No Hydrogen btn",True).get_content()
        no_hydrogen_button.register_pressed_callback(no_hydrogen_button_pressed_callback)

        # create the use reflection button
        # use_reflections_button = menu.root.find_node("Use Reflection btn",True).get_content()
        # use_reflections_button.register_pressed_callback(use_reflections_button_pressed_callback)

        # create the backbone only button
        backbone_only_button = menu.root.find_node("Backbone only btn",True).get_content()
        backbone_only_button.register_pressed_callback(backbone_only_button_pressed_callback)

        # create the selected only button
        selected_only_button =  menu.root.find_node("Selected Only btn",True).get_content()
        selected_only_button.register_pressed_callback(selected_only_button_pressed_callback)

        # create the reorder button
        reorder_button = menu.root.find_node("Reorder menu",True).get_content()
        reorder_button.register_pressed_callback(reorder_button_pressed_callback)

        # create the roation "drop down"
        rotation_button = menu.root.find_node("Rotation menu",True).get_content()
        rotation_button.register_pressed_callback(rotation_button_pressed_callback)

        # create the rmsd score
        self.rmsd_score_label = menu.root.find_node("RMSD number",True).get_content()

        # create the receptor text
        self.receptor_text = menu.root.find_node("Receptor").get_content()
        
        # create the target text
        self.target_text = menu.root.find_node("Target").get_content()
        
        self._menu = menu
        

        # self._request_refresh()