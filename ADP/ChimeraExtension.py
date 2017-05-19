# -----------------------------------------------------------------------------
#
from chimera.extension import EMO, manager

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


# -----------------------------------------------------------------------------
#
# def fit_map(cmdname, args):
#     from FitMap import fitcmd
#     fitcmd.fitmap_command(cmdname, args)
# def unfit_map(cmdname, args):
#     from FitMap.move import position_history
#     position_history.undo()
#
# from Midas import midas_text as mt
# mt.addCommand('fitmap', fit_map, unfit_map, help = True)

# -----------------------------------------------------------------------------
# #
# def fit_map_cb():
#     from FitMap import fitmap as F
#     F.move_selected_atoms_to_maximum()
# def fit_map_rotation_only_cb():
#     from FitMap import fitmap as F
#     F.move_selected_atoms_to_maximum(optimize_translation = False)
# def fit_map_shift_only_cb():
#     from FitMap import fitmap as F
#     F.move_selected_atoms_to_maximum(optimize_rotation = False)
# def move_atoms_to_maxima():
#     from FitMap import fitmap as F
#     F.move_atoms_to_maxima()
#
# from Accelerators import add_accelerator
# add_accelerator('ft', 'Move model to maximize density at selected atoms',
#                 fit_map_cb)
# add_accelerator('fr', 'Rotate model to maximize density at selected atoms',
#                 fit_map_rotation_only_cb)
# add_accelerator('fs', 'Shift model to maximize density at selected atoms',
#                 fit_map_shift_only_cb)
# add_accelerator('mX', 'Move selected atoms to local maxima',
#                 move_atoms_to_maxima)
