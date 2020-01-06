import nanome
from nanome.util import Logs
from nanome.util import Color

import os

SELECTED_COLOR = Color.from_int(0x00ecc4ff)
DESELECTED_COLOR = Color.from_int(0xffffffff)
CHECKICON =  "GreenCheck.png"
LOCKICON = "Lock.png"
UNLOCKICON = "Unlock.png"
REFRESHICON = "Refresh.png"
QUESTIONMARKICON = "QuestionMark.png"

class RMSDMenu():
    def __init__(self, rmsd_plugin):
        self._menu = rmsd_plugin.menu
        self._plugin = rmsd_plugin
        self._selected_mobile = [] # button
        self._selected_target = None # button
        self._run_button = None
        self._current_tab = "receptor" #receptor = 0, target = 1
        self._drop_down_dict={"rotation":["None", "Kabsch","Quaternion"],"reorder_method":["None","Hungarian","Brute", "Distance"],\
        "select":["None","Global"]} # select["Local"] in the future
        self._current_reorder = "None"
        self._current_rotation = "None"
        self._current_select = "None"

    def _request_refresh(self):
        self._plugin.request_refresh()

    # run the rmsd algorithm
    def _run_rmsd(self):
        if self.check_resolve_error():
            self._plugin.update_structures_deep(self._selected_mobile + [self._selected_target])
            self._plugin.run_rmsd([a.complex for a in self._selected_mobile], self._selected_target.complex)
        else:
            self.hide_loading_bar()

    # check the "unselect" and "select_same" error and call change_error
    def check_resolve_error(self,clear_only=False):
        if(not clear_only):
            if self._selected_mobile == None or self._selected_target == None:
                self.change_error("unselected")
                return False
            #elif self._selected_mobile.complex.index == self._selected_target.complex.index:
            elif self._selected_target.complex.index in list(map(lambda a:a.complex.index,self._selected_mobile)):
                self.change_error("select_same")
                return False
            else:
                self.change_error("clear")
                return True
        else:
            self.change_error("clear")
            return True

    def change_loading_percentage(self,percentage):
        self.loadingBar.percentage = percentage
        self._plugin.update_content(self.loadingBar)

    # show the error message texts, fromRun means if the it is called after Run is pressed
    def change_error(self,error_type):
        if(error_type == "unselected"):
            self.error_message.text_auto_size=False
            self.error_message.text_size = 0.198
            self.error_message.text_value = "Select both target and receptor"
            self.update_score(None)
            self._plugin.update_content(self.error_message)
            return True
        if(error_type == "select_same"):
            self.error_message.text_auto_size=False
            self.error_message.text_size = 0.22
            self.error_message.text_value = "Select different complexes"
            self.update_score()
            self._plugin.update_content(self.error_message)
            return True
        if(error_type == "different_size"):
            self.error_message.text_auto_size=False
            self.error_message.text_size = 0.159
            self.error_message.text_value = "Receptor and target have different sizes"
            self.update_score()
            self._plugin.update_content(self.error_message)
            return True
        if(error_type == "different_order"):
            self.error_message.text_auto_size=False
            self.error_message.text_size = 0.159
            self.error_message.text_value = "Receptor and target have different order"
            self.update_score()
            self._plugin.update_content(self.error_message)
            return True
        if(error_type == "zero_size"):
            self.error_message.text_auto_size=False
            self.error_message.text_size = 0.15
            self.error_message.text_value = "At least one complex has no atom selected"
            self.update_score()
            self._plugin.update_content(self.error_message)
            return True
        if(error_type == "selected_changed"):
            self.error_message.text_auto_size=False
            self.error_message.text_size = 0.176
            self.error_message.text_value = "Selected complexes have changed"
            self.update_score()
            self._plugin.update_content(self.error_message)
            return True
        if(error_type == "loading"):
            self.error_message.text_auto_size=False
            self.error_message.text_size = 0.2
            self.error_message.text_value = "Loading..."
            self.update_score()
            self._plugin.update_content(self.error_message)
            return True
        if(error_type == "clear"):
            self.error_message.text_value = ""
            self.error_message.text_auto_size = True
            self._plugin.update_content(self.error_message)


    # change the args in the plugin
    def update_args(self,arg,option):
        self._plugin.update_args(arg,option)

    def update_score(self,value=None):
        Logs.debug("update score called: ",value)
        if value == None:
            self.rmsd_score_label.text_value = "--"
        else:
            self.rmsd_score_label.text_value = str("%.3g"%value)
        self._plugin.update_content(self.rmsd_score_label)

    def make_plugin_usable(self, state = True):
        self._run_button.unusable = not state
        self._plugin.update_button(self._run_button)

    def show_loading_bar(self):
        self.ln_select_run.enabled = False
        self.ln_loading_bar.enabled = True
        self._plugin.update_menu(self._menu)

    def hide_loading_bar(self):
        self.ln_select_run.enabled = True
        self.ln_loading_bar.enabled = False
        self._plugin.update_menu(self._menu)

    # change the complex list
    def change_complex_list(self, complex_list):
        
        # a button in the receptor list is pressed
        def mobile_pressed(button):

            # selecting button
            if button.complex.index not in [ a.complex.index for a in self._selected_mobile]:
                button.selected = True
                self._selected_mobile.append(button)
                self._plugin.update_mobile([x.complex for x in self._selected_mobile])
                if len(self._selected_mobile) == 1:
                    self.receptor_text.text_value = "Receptor: "+button.complex.name
                else:
                    self.receptor_text.text_value = "Receptor: multiple receptors"
                self.receptor_check.add_new_image(file_path = os.path.join(os.path.dirname(__file__), CHECKICON))
               

            # deselecting button
            else:
                button.selected = False
                # self._selected_mobile = [i for i in self._selected_mobile if i.copmplex.index != button.complex.index]
                for x in self._selected_mobile:
                    if x.complex.index == button.complex.index:
                        self._selected_mobile.remove(x)
                if len(self._selected_mobile) == 1:
                    self.receptor_text.text_value = "Receptor: "+self._selected_mobile[0].complex.name
                elif len(self._selected_mobile) == 0:
                    self.receptor_text.text_value = "Receptor: Unselected"
                    self.receptor_check.add_new_image(file_path = os.path.join(os.path.dirname(__file__),QUESTIONMARKICON))

                else:
                    self.receptor_text.text_value = "Receptor: multiple receptors"
                    self.receptor_check.add_new_image(file_path = os.path.join(os.path.dirname(__file__), CHECKICON))

            # tell the plugin and update the menu
            self._current_select = "None"
            # self.update_args("select", "None")
            self._plugin.update_mobile([x.complex for x in self._selected_mobile])
            self._plugin.update_menu(self._menu)

        # a button in the target list is pressed
        def target_pressed(button):
            if self._selected_target != None:
                self._selected_target.selected = False 
                if self._selected_target.complex != button.complex: 
                    button.selected = True
                    self._selected_target = button
                    self._plugin.update_target(self._selected_target.complex)
                    self.target_text.text_value ="Target: "+ button.complex.name
                    self.target_check.add_new_image(file_path = os.path.join(os.path.dirname(__file__), CHECKICON))
                  

                else: 
                    self._selected_target = None
                    self.target_text.text_value = "Target: Unselected"
                    self.target_check.add_new_image(file_path = os.path.join(os.path.dirname(__file__), QUESTIONMARKICON))

            else: 
                button.selected = True
                self._selected_target = button
                self._plugin.update_target(self._selected_target.complex)
             
                self.target_text.text_value ="Target: "+ button.complex.name
                # still setting the image just in case theres a bug
                self.target_check.add_new_image(file_path = os.path.join(os.path.dirname(__file__), CHECKICON))

            self.check_resolve_error(clear_only=True)
            self._current_select = "None"
            # self.update_args("select", "None")
            if self._selected_target == None:
                self._plugin.update_target(None)
            else:
                self._plugin.update_target(self._selected_target.complex)
           
            self._plugin.update_menu(self._menu)

        self._mobile_list = []
        self._target_list = []

        if self._selected_mobile == None:
            self._selected_mobile = []

        if len(self._selected_mobile) != 0:
            for x in self._selected_mobile:
                if x.complex.index not in [i.index for i in complex_list]:
                    self._selected_mobile.remove()
            self._selected_mobile =None
        if self._selected_target != None and \
           self._selected_target.complex.index not in [i.index for i in complex_list]:
            self._selected_target = None
        
        for complex in complex_list:
            clone = self._complex_item_prefab.clone()
            ln_btn = clone.get_children()[0]
            btn = ln_btn.get_content()
            btn.text.value.set_all(complex.name)
            btn.complex = complex
            btn.register_pressed_callback(mobile_pressed)
            self._mobile_list.append(clone)
            if self._selected_mobile != None and btn.complex.index in [a.complex.index for a in self._selected_mobile]:
                btn.selected = True

            #clone1 = clone.clone()
            clone1 = self._complex_item_prefab.clone()
            ln_btn = clone1.get_children()[0]
            btn = ln_btn.get_content()
            btn.text.value.set_all(complex.name)
            btn.complex = complex
            btn.register_pressed_callback(target_pressed)
            self._target_list.append(clone1)
            if self._selected_target != None and btn.complex.index == self._selected_target.complex.index:
                btn.selected = True

        if len(self._selected_mobile) == 0:
            self.receptor_text.text_value ="Receptor: Unselected"
        if self._selected_target == None:
            self.target_text.text_value ="Target: Unselected "
        if self._current_tab == "receptor":
            self._show_list.items=self._mobile_list
        else:
            self._show_list.items=self._target_list
            
        self._plugin.update_menu(self._menu)

    # change the lock image to Lock
    def lock_image(self):
        self.lock_img.add_new_image(file_path = os.path.join(os.path.dirname(__file__), LOCKICON))
        self._plugin.update_menu(self._menu)

    # build the menu
    def build_menu(self):
        # refresh the lists
        def refresh_button_pressed_callback(button):
            self._request_refresh()
            
        # press the run button and run the algorithm
        def run_button_pressed_callback(button):
            # self.show_loading_bar()
            if self._selected_mobile != None and self._selected_target != None:
                self.show_loading_bar()
                self._plugin.select([x.complex for x in self._selected_mobile],self._selected_target.complex)
               
            else:
                self.check_resolve_error()
        
        # press the lock button and lock/unlock the complexes
        def lock_button_pressed_callback(button):
            def toggle_lock(complex_list):
                new_locked = not all(elem.locked for elem in complex_list)

                for x in complex_list:
                    x.locked = new_locked
                    # x.boxed = new_locked
                
                self.lock_img.add_new_image(os.path.join(os.path.dirname(__file__), (LOCKICON if new_locked else UNLOCKICON)))
                self._plugin.update_menu(self._menu)
                self._plugin.update_structures_shallow(complex_list)
                
            if self._selected_target != None and len(self._selected_mobile) != 0:
                complex_list = self._plugin._mobile + [self._plugin._target]
                complex_indexes = [complex.index for complex in complex_list]
                self._plugin.request_complexes(complex_indexes, toggle_lock)

            else:
                self.change_error("unselected")
            
        # show the target list when the receptor tab is pressed
        def receptor_tab_pressed_callback(button):
            self._current_tab="receptor"
            receptor_tab.selected = True
            target_tab.selected = False
            self._show_list.items = self._mobile_list
            self._plugin.update_content(receptor_tab)
            self._plugin.update_content(target_tab)
            self._plugin.update_content(self._show_list)

        # show the target list when the target tab is pressed
        def target_tab_pressed_callback(button):
            self._current_tab="target"
            target_tab.selected = True
            receptor_tab.selected = False
            self._show_list.items = self._target_list
            self._plugin.update_content(receptor_tab)
            self._plugin.update_content(target_tab)
            self._plugin.update_content(self._show_list)

        # align sequence = ! align sequence
        def align_sequence_button_pressed_callback(button):
            align_sequence_button.selected = not align_sequence_button.selected
            if align_sequence_button.selected:
                align_sequence_text.text_color = SELECTED_COLOR
            else:
                align_sequence_text.text_color = DESELECTED_COLOR
            self.update_args("align_sequence", align_sequence_button.selected)
            self._current_select = "global" if align_sequence_button.selected else "none"
            self._plugin.update_content(align_sequence_button)
            self._plugin.update_content(align_sequence_text)

        # no hydrogen = ! no hydrogen
        def no_hydrogen_button_pressed_callback(button):
            button.selected = not button.selected
            if button.selected:
                no_hydrogen_text.text_color = SELECTED_COLOR
            else:
                no_hydrogen_text.text_color = DESELECTED_COLOR     
            self.update_args("no_hydrogen", button.selected)
            self._plugin.update_content(button)
            self._plugin.update_content(no_hydrogen_text)

        # backbone only = ! backbone only
        def backbone_only_button_pressed_callback(button):
            button.selected = not button.selected
            if button.selected:
                backbone_only_text.text_color = SELECTED_COLOR
            else:
                backbone_only_text.text_color = DESELECTED_COLOR     
            self.update_args("backbone_only", button.selected)
            self._plugin.update_content(button)
            self._plugin.update_content(backbone_only_text)
        
        # selected only = ! selected only
        def selected_only_button_pressed_callback(button):
            button.selected = not button.selected 
            if button.selected:
                selected_only_text.text_color = SELECTED_COLOR
            else:
                selected_only_text.text_color = DESELECTED_COLOR     
            self.update_args("selected_only", button.selected)
            self._plugin.update_content(button)
            self._plugin.update_content(selected_only_text)

        # align box = ! align box
        def align_box_button_pressed_callback(button):
            button.selected = not button.selected
            if button.selected:
                align_box_text.text_color = SELECTED_COLOR
            else:
                align_box_text.text_color = DESELECTED_COLOR     
            self.update_args("align_box",button.selected)
            self._plugin.update_content(button)
            self._plugin.update_content(align_box_text)

        # global <=> local
        def global_local_button_pressed_callback(button):
            if self._plugin.args.select == "global":
                self.update_args("select","local")
                button.text.value.set_all("Local")
            else:
                self.update_args("select","global")
                button.text.value.set_all("Global")
            self._plugin.update_content(button)

        # change Reorder to the next option
        def reorder_button_pressed_callback(button):
            drop_down  = self._drop_down_dict["reorder_method"]
            temp_length=len(drop_down)
            
            pre_index = drop_down.index(self._current_reorder)
            post_index = (pre_index + 1) % temp_length

            post_option = drop_down[post_index]

            button.selected = post_option != "None"
            button.text.value.set_all(post_option)
            
            if post_option == "None":
                reorder_text.text_color = DESELECTED_COLOR     
            else:
                reorder_text.text_color = SELECTED_COLOR     

            # tell the plugin and update the menu
            self._current_reorder = post_option
            self.update_args("reorder_method", post_option)
            self.update_args("reorder", post_option != "None")
            self._plugin.update_content(button)
            self._plugin.update_content(reorder_text)

        # change Rotation to the next option
        def rotation_button_pressed_callback(button):
            drop_down  = self._drop_down_dict["rotation"]
            temp_length=len(drop_down)
            
            pre_index = drop_down.index(self._current_rotation)
            post_index = (pre_index + 1) % temp_length

            post_option = drop_down[post_index]

            button.selected = post_option != "None"
            button.text.value.set_all(post_option)
            
            if post_option == "None":
                rotation_text.text_color = DESELECTED_COLOR     
            else:
                rotation_text.text_color = SELECTED_COLOR 

            # tell the plugin and update the menu
            self._current_rotation = post_option
            self.update_args("rotation_method", post_option)
            self._plugin.update_content(button)
            self._plugin.update_content(button)
            self._plugin.update_content(rotation_text)


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

        # add the refresh icon
        refresh_img = menu.root.find_node("Refresh Image", True)
        refresh_img.add_new_image(file_path = os.path.join(os.path.dirname(__file__), REFRESHICON))

        # add the receptor check icon
        self.receptor_check = menu.root.find_node("Receptor Check", True)
        self.receptor_check.add_new_image(file_path = os.path.join(os.path.dirname(__file__), QUESTIONMARKICON))

        # add the target check icon
        self.target_check = menu.root.find_node("Target Check", True)
        self.target_check.add_new_image(file_path = os.path.join(os.path.dirname(__file__), QUESTIONMARKICON))

        # create the layout node that contains select and run and refresh
        self.ln_select_run = menu.root.find_node("Refresh Run",True)

        # create the Run button
        self._run_button = menu.root.find_node("Run", True).get_content()
        self._run_button.register_pressed_callback(run_button_pressed_callback)

        # create the Refresh button
        refresh_button = menu.root.find_node("Refresh Button", True).get_content()
        refresh_button.register_pressed_callback(refresh_button_pressed_callback)

        # create the lock button
        lock_button = menu.root.find_node("Lock Button",True).get_content()
        lock_button.register_pressed_callback(lock_button_pressed_callback)

        # add the lock icon,
        self.lock_img = menu.root.find_node("Lock Image",True)
        mobile_locked = False
        for x in self._selected_mobile:
            if x.locked:
                mobile_locked = True

        if len(self._selected_mobile) != 0 and self._selected_target != None and self._selected_target.locked and mobile_locked:
            self.lock_img.add_new_image(file_path = os.path.join(os.path.dirname(__file__), LOCKICON))
        else:
            self.lock_img.add_new_image(file_path = os.path.join(os.path.dirname(__file__), UNLOCKICON))


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

        # create the align seqeuence button
        align_sequence_button = menu._root.find_node("Align sequence btn",True).get_content()
        align_sequence_button.register_pressed_callback(align_sequence_button_pressed_callback)
        align_sequence_text = menu.root.find_node("Align sequence txt",True).get_content()

        # create the no hydrogen button
        no_hydrogen_button = menu.root.find_node("No Hydrogen btn",True).get_content()
        no_hydrogen_button.register_pressed_callback(no_hydrogen_button_pressed_callback)
        no_hydrogen_text = menu.root.find_node("No Hydrogen txt",True).get_content()

        # create the use reflection button
        # use_reflections_button = menu.root.find_node("Use Reflection btn",True).get_content()
        # use_reflections_button.register_pressed_callback(use_reflections_button_pressed_callback)

        # create the backbone only button
        backbone_only_button = menu.root.find_node("Backbone only btn",True).get_content()
        backbone_only_button.register_pressed_callback(backbone_only_button_pressed_callback)
        backbone_only_text = menu.root.find_node("Backbone only txt",True).get_content()


        # create the selected only button
        selected_only_button =  menu.root.find_node("Selected Only btn",True).get_content()
        selected_only_button.register_pressed_callback(selected_only_button_pressed_callback)
        selected_only_text = menu.root.find_node("Selected Only txt",True).get_content()

        # create the global/local button
        global_local_button =  menu.root.find_node("Global Local menu",True).get_content()
        global_local_button.register_pressed_callback(global_local_button_pressed_callback)
        global_local_text = menu.root.find_node("Global Local txt",True).get_content()


        # create the align box button
        align_box_button =  menu.root.find_node("Box btn",True).get_content()
        align_box_button.register_pressed_callback(align_box_button_pressed_callback)
        align_box_text = menu.root.find_node("Box txt",True).get_content()

        # create the reorder button
        reorder_button = menu.root.find_node("Reorder menu",True).get_content()
        reorder_button.register_pressed_callback(reorder_button_pressed_callback)
        reorder_text = menu.root.find_node("Reorder txt",True).get_content()

        # create the roation "drop down"
        rotation_button = menu.root.find_node("Rotation menu",True).get_content()
        rotation_button.register_pressed_callback(rotation_button_pressed_callback)
        rotation_text = menu.root.find_node("Rotation txt",True).get_content()


        self.ln_loading_bar = menu.root.find_node("Loading Bar",True)
        self.ln_loading_bar.forward_dist = .003
        self.loadingBar = self.ln_loading_bar.add_new_loading_bar()
        self.loadingBar.description = "      Loading...          "

        # create the rmsd score
        self.rmsd_score_label = menu.root.find_node("RMSD number",True).get_content()

        # create the receptor text
        self.receptor_text = menu.root.find_node("Receptor").get_content()
        
        # create the target text
        self.target_text = menu.root.find_node("Target").get_content()
        
        # create the error message text
        error_node = menu.root.find_node("Error Message")
        self.error_message = error_node.get_content()

        self._menu = menu
        

        # self._request_refresh()