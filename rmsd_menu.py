import nanome

class RMSDMenu():
    def __init__(self, rmsd_plugin):
        self._menu = rmsd_plugin.menu
        self._plugin = rmsd_plugin
        self._selected_mobile = None # Tuple (button/complex)
        self._selected_target = None # Tuple (button/complex)
        self._complex_list = []
        self._run_button = None

    def _request_refresh(self):
        self._plugin.request_refresh()

    def _run_rmsd(self):
        self._plugin.run_rmsd(self._selected_mobile[1], self._selected_target[1])

    def _find_complex(self, name):
        for complex in self._complex_list:
            if complex.name == name:
                return complex
        return None

    def make_plugin_usable(self, state = True):
        self._run_button.unusable = not state
        self._plugin.update_button(self._run_button)

    def change_complex_list(self, complex_list):
        def mobile_pressed(button):
            complex = self._find_complex(button.text.value_idle)
            if complex == None:
                nanome.util.Logs.error("Couldn't retrieve a complex from its button text")
                return

            if self._selected_mobile != None:
                self._selected_mobile[0].selected = False
            button.selected = True
            self._selected_mobile = (button, complex)
            self._plugin.update_menu(self._menu)

        def target_pressed(button):
            complex = self._find_complex(button.text.value_idle)
            if complex == None:
                nanome.util.Logs.error("Couldn't retrieve a complex from its button text")
                return

            if self._selected_target != None:
                self._selected_target[0].selected = False
            button.selected = True
            self._selected_target = (button, complex)
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
            self._mobile_list.items.append(clone)
            
            clone1 = clone.clone()
            ln_btn = clone1.get_children()[0]
            btn = ln_btn.get_content()
            btn.register_pressed_callback(target_pressed)
            self._target_list.items.append(clone1)
        self._plugin.update_menu(self._menu)

    def build_menu(self):
        def refresh_button_pressed_callback(button):
            self._request_refresh()

        def run_button_pressed_callback(button):
            self.make_plugin_usable(False)
            self._run_rmsd()

        # Create a prefab that will be used to populate the lists
        self._complex_item_prefab = nanome.ui.LayoutNode()
        self._complex_item_prefab.layout_orientation = nanome.ui.LayoutNode.LayoutTypes.horizontal
        child = self._complex_item_prefab.create_child_node()
        child.name = "complex_button"
        prefabButton = child.add_new_button()
        prefabButton.text.active = True

        menu = nanome.ui.Menu.io.from_json("rmsd_menu.json")
        self._plugin.menu = menu

        self._run_button = menu.root.find_node("Run", True).get_content()
        self._run_button.register_pressed_callback(run_button_pressed_callback)

        refresh_button = menu.root.find_node("Refresh", True)
        refresh_button.get_content().register_pressed_callback(refresh_button_pressed_callback)

        self._mobile_list = menu.root.find_node("Mobile_List", True).get_content()
        self._target_list = menu.root.find_node("Target_List", True).get_content()

        self._menu = menu
        # self._request_refresh()