# -----------------------------------------------------------------------------
# Class to register the ADP EM plugin in Chimera.
# The plugin will be located in EM Fitting/Volume Data menu
#
import chimera.extension

# -----------------------------------------------------------------------------
#
class ADPEM_EMO(chimera.extension.EMO):

    def name(self):
        return 'ADP EM'
    def description(self):
        return 'Fitting density map in an exhaustive way through ADP_EM Algorithm'
    def categories(self):
        return ['EM Fitting']
    def icon(self):
        return None
    def activate(self):
        self.module('adpgui').show_adp_em_dialog()
        return None

chimera.extension.manager.registerExtension(ADPEM_EMO(__file__))