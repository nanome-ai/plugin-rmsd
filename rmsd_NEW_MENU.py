import nanome

class RMSDMenu():
    def __init__(self, rmsd_plugin):
        self._menu = rmsd_plugin.menu
        self._plugin = rmsd_plugin
        self._selected_mobile = None # button
        self._selected_target = None # button
        self._complex_list = []
        self._run_button = None
        self._current_tab = "receptor" #receptor = 0, target = 1
        self._drop_down_dict={"rotation":["None", "Kabsch","Quaternion"],"reorder_method":["None","Hungarian","Brute", "Distance"]}
        self._current_reorder = "None"
        self._current_rotation = "None"

    def _request_refresh(self):
        self._plugin.request_refresh()

    def _run_rmsd(self):
        if self._selected_mobile != None or self._selected_target != None:
            self._plugin.run_rmsd(self._selected_mobile.complex, self._selected_target.complex)
        else:
            # maybe an error message?
            pass

    # change the args 
    def update_args(self,arg,option):
        self._plugin.update_args(arg,option)

    def make_plugin_usable(self, state = True):
        self._run_button.unusable = not state
        self._plugin.update_button(self._run_button)

    def change_complex_list(self, complex_list):
        def mobile_pressed(button):
            if self._selected_mobile != None:
                self._selected_mobile.selected = False
            button.selected = True
            self._selected_mobile = button
            self._plugin.update_menu(self._menu)

        def target_pressed(button):
            if self._selected_target != None:
                self._selected_target.selected = False
            button.selected = True
            self._selected_target = button
            self._plugin.update_menu(self._menu)

        self._complex_list = complex_list
        self._selected_mobile = None
        self._selected_target = None
        self._mobile_list.items = []
        self._target_list.items = []

        for complex in complex_list:
            clone = self._complex_item_prefab.clone()
            ln_btn = clone.get_children()[0]
            btn = ln_btn.get_content()
            btn.set_all_text(complex.name)
            btn.register_pressed_callback(mobile_pressed)
            btn.complex = complex
            self._mobile_list.items.append(clone)
            
            clone1 = clone.clone()
            ln_btn = clone1.get_children()[0]
            btn = ln_btn.get_content()
            btn.register_pressed_callback(target_pressed)
            btn.complex = complex
            self._target_list.items.append(clone1)
        self._plugin.update_menu(self._menu)

    def build_menu(self):
        def refresh_button_pressed_callback(button):
            self._request_refresh()

        def run_button_pressed_callback(button):
            self.make_plugin_usable(False)
            self._run_rmsd()

        def receptor_tab_pressed_callback(button):
            self._current_tab="receptor"
            receptor_tab.selected = True
            target_tab.selected = False
            self._plugin.update_menu(self._menu)

        def target_tab_pressed_callback(button):
            self._current_tab="target"
            target_tab.selected = True
            receptor_tab.selected = False
            self._plugin.update_menu(self._menu)
            
        
        # no hydrogen = ! no hydrogen
        def no_hydrogen_button_pressed_callback(button):
            self.update_args("no_hydrogen",0)
            no_hydrogen_button.selected = not no_hydrogen_button.selected
            self._plugin.update_menu(self._menu)

        # use reflections = ! use reflections
        def use_reflections_button_pressed_callback(button):
            self.update_args("use_reflections",0)
            use_reflections_button.selected = not use_reflections_button.selected
            self._plugin.update_menu(self._menu)
        
        # selected only = ! selected only
        def selected_only_button_pressed_callback(button):
            self.update_args("selected_only",0)
            selected_only_button.selected = not selected_only_button.selected            
            self._plugin.update_menu(self._menu)

        # change Reorder to the next option
        def reorder_button_pressed_callback(button):
            temp_length=len(self._drop_down_dict["reorder_method"])
            temp_before = self._current_reorder
            if self._drop_down_dict["reorder_method"].index(self._current_reorder) == temp_length-1:
                self._current_reorder = self._drop_down_dict["reorder_method"][0]
            else:
                self._current_reorder = self._drop_down_dict["reorder_method"][(self._drop_down_dict["reorder_method"].index(self._current_reorder)+1)]
            # if from None to not None, selected = ! selected
            if temp_before == "None":
                reorder_button.selected = not reorder_button.selected
            # if from last not None to None, selected = ! selected,
            # selected text = first not None
            elif self._drop_down_dict["reorder_method"].index(temp_before) == temp_length-1:
                reorder_button.selected = not reorder_button.selected
                reorder_button.text.value_selected = self._drop_down_dict["reorder_method"][1]
                reorder_button.text.value_selected_highlighted = self._drop_down_dict["reorder_method"][1]

            # if from not None to not None, selected text = self._current_reorder
            else:
                reorder_button.text.value_selected = self._current_reorder
                reorder_button.text.value_selected_highlighted = self._current_reorder
                
            self.update_args("reorder_method",self._current_reorder)
            self._plugin.update_menu(self._menu)

        # change Rotation to the next option
        def rotation_button_pressed_callback(button):
            temp_length=len(self._drop_down_dict["rotation"])
            temp_before = self._current_rotation
            if self._drop_down_dict["rotation"].index(self._current_rotation) == temp_length-1:
                self._current_rotation = self._drop_down_dict["rotation"][0]
            else:
                self._current_rotation = self._drop_down_dict["rotation"][(self._drop_down_dict["rotation"].index(self._current_rotation)+1)]
            # if from None to not None, selected = ! selected
            if temp_before == "None":
                rotation_button.selected = not rotation_button.selected
            # if from last not None to None, selected = ! selected,
            # selected text = first not None
            elif self._drop_down_dict["rotation"].index(temp_before) == temp_length-1:
                rotation_button.selected = not rotation_button.selected
                rotation_button.text.value_selected = self._drop_down_dict["rotation"][1]
                rotation_button.text.value_selected_highlighted = self._drop_down_dict["rotation"][1]

            # if from not None to not None, selected text = self._current_rotation
            else:
                rotation_button.text.value_selected = self._current_rotation
                rotation_button.text.value_selected_highlighted = self._current_rotation
            # "text_value_idle":"None",
            # "text_value_selected":"Kabsch","Quaternion",
            # "text_value_highlighted":"None",
            # "text_value_selected_highlighted":"Kabsch","Quaternion",
            self.update_args("rotation",self._current_rotation)
            self._plugin.update_menu(self._menu)

        
        # Create a prefab that will be used to populate the lists
        self._complex_item_prefab = nanome.ui.LayoutNode()
        self._complex_item_prefab.layout_orientation = nanome.ui.LayoutNode.LayoutTypes.horizontal
        child = self._complex_item_prefab.create_child_node()
        child.name = "complex_button"
        prefabButton = child.add_new_button()
        prefabButton.text.active = True

        menu = nanome.ui.Menu.io.from_json("rmsd_pluginator.json")
        self._plugin.menu = menu

        self._run_button = menu.root.find_node("Run", True).get_content()
        self._run_button.register_pressed_callback(run_button_pressed_callback)

        refresh_button = menu.root.find_node("Refresh", True).get_content()
        refresh_button.register_pressed_callback(refresh_button_pressed_callback)
        self._mobile_list = menu.root.find_node("List", True).get_content()
        self._target_list = menu.root.find_node("List", True).get_content()

        receptor_tab = menu.root.find_node("Receptor_tab",True).get_content()
        receptor_tab.register_pressed_callback(receptor_tab_pressed_callback)
        target_tab = menu.root.find_node("Target_tab",True).get_content()
        target_tab.register_pressed_callback(target_tab_pressed_callback)
        no_hydrogen_button = menu.root.find_node("No Hydrogen btn",True).get_content()
        no_hydrogen_button.register_pressed_callback(no_hydrogen_button_pressed_callback)
        use_reflections_button = menu.root.find_node("Use Reflection btn",True).get_content()
        use_reflections_button.register_pressed_callback(use_reflections_button_pressed_callback)
        selected_only_button =  menu.root.find_node("Selected Only btn",True).get_content()
        selected_only_button.register_pressed_callback(selected_only_button_pressed_callback)
        reorder_button = menu.root.find_node("Reorder menu",True).get_content()
        reorder_button.register_pressed_callback(reorder_button_pressed_callback)
        rotation_button = menu.root.find_node("Rotation menu",True).get_content()
        rotation_button.register_pressed_callback(rotation_button_pressed_callback)
        rmsd_score_label = menu.root.find_node("RMSD number",True)
        
        self._menu = menu
        

        # self._request_refresh()