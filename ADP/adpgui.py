# ---------------------------------------------------------------------------------
# Dialog to perform rigid fitting through ADP EM algorithm.
#

import chimera
from chimera.baseDialog import ModelessDialog
import ADP

# ---------------------------------------------------------------------------------
# ADP EM Plugin Dialog
#
class ADP_EM_Dialog(ModelessDialog):

  # Title of APD EM plugin
  title = 'ADP EM Exhaustive Fitting'
  # Name of ADP EM plugin
  name = 'ADP EM'
  # Buttons of ADP EM GUI
  buttons = ('Fit', 'Options', 'Results', 'Close')
  # Path of help guide of ADP EM plugin
  help = ('adp_em.html', ADP)
  # Name of the folder where ADP EM plugin is located
  plugin_folder = 'ADP/'
  # Path of the process of ADP EM
  adp_em = plugin_folder + 'adp_em'
  # Steps for the process progress bar
  steps = ["Step 1", "Step 2", "Step 3"]
  # Steps marks for the process progress bar
  marks = ["ADP_EM", "Trans limit", "Interpolation", "Total Time"]
  # Variable to represent the dropdown menu that will show all the ADP EM Solutions
  mb = None
  # Array to store the solutions generated by ADP EM
  solutions_chimera = []
  # Variables to handle a possible movement of the molecule by the user
  # They will store an Xform object (gives the rotation and traslation for a model)
  # relative to the pdb, the map, its inverse or the Xform of the last solution chosen
  # It works as a movements traceback
  bos = None
  xf = None
  xfC = None
  xform_last_solution = None

  #-----------------------
  # ADP EM Chimera Commands

  # ADP in Chimera indicator
  adp_em_chimera_opt = "--chimera"

  # No pdb solutions files generated
  adp_em_chimera_no_save = "--no_save"

  # Fitting criterion
  adp_em_chimera_laplacian = "-f"
  adp_em_chimera_laplacian_val = "1"

  # Saved solutions
  adp_em_chimera_saved_solutions = "-n"
  adp_em_chimera_saved_solutions_val = "50"

  # Translational sampling
  adp_em_chimera_sampling = "-t"
  adp_em_chimera_sampling_val = "2"

  # Translational scan
  adp_em_chimera_scan = "-s"
  adp_em_chimera_scan_val = "2"

  # Peaks explored per docking
  adp_em_chimera_peaks_explored = "--ne"
  adp_em_chimera_peaks_explored_val = "30"

  # Peaks stored per iteration
  adp_em_chimera_peaks_iteration = "--np"
  adp_em_chimera_peaks_iteration_val = "20"

  # Peaks stored per search
  adp_em_chimera_peaks_search = "--nr"
  adp_em_chimera_peaks_search_val = "100"

  # Peaks stored per multidocking search
  adp_em_chimera_peaks_msearch = "--nrm"
  adp_em_chimera_peaks_msearch_val = "500"

  # Translational threshold
  adp_em_chimera_translational_threshold = "--rt"
  adp_em_chimera_translational_threshold_val = "2.0"

  # Rotational threshold
  adp_em_chimera_rotational_threshold = "--rc"
  adp_em_chimera_rotational_threshold_val = None

  # Width between spherical layers
  adp_em_chimera_width_layers = "--lw"
  adp_em_chimera_width_layers_val = "1.0"

  # Density cutoff of simulated map
  adp_em_chimera_cutoff_simulated = "--cutoff2"
  adp_em_chimera_cutoff_simulated_val = "0.0"


  #-------------------------------------------------
  # Master function for dialog contents
  #
  def fillInUI(self, parent):

    t = parent.winfo_toplevel()
    self.toplevel_widget = t
    t.withdraw()

    parent.columnconfigure(0, weight = 1)
    row = 0

    import Tkinter
    from CGLtk import Hybrid
    from VolumeViewer import Volume_Menu

    ff = Tkinter.Frame(parent)
    ff.grid(row = row, column = 0, sticky = 'w')
    row = row + 1

    # Files selection (only Molecules)
    from chimera import Molecule
    mlist = [m for m in fit_object_models() if isinstance(m, Molecule)]
    fstart = mlist[0] if mlist else None
    from chimera.widgets import ModelOptionMenu
    om = ModelOptionMenu(ff, labelpos = 'w', label_text = 'Fit ',
                         initialitem = fstart,
                         listFunc = fit_object_models,
                         sortFunc = compare_fit_objects,
                         command = self.object_chosen_cb)
    om.grid(row = 0, column = 0, sticky = 'w')
    self.object_menu = om

    # Maps selection (only Volumes)
    fm = Volume_Menu(ff, ' in map ')
    fm.frame.grid(row = 0, column = 1, sticky = 'w')
    self.map_menu = fm

    hf = Tkinter.Frame(parent)
    hf.grid(row=row, column=0, sticky='w')
    row += 1

    # Bandwidth
    from CGLtk import Hybrid
    bandwidths = ["16", "24", "32", "48", "64"]
    self.bandwidth = apply(Hybrid.Option_Menu, (hf, 'Bandwidth ')  + tuple(bandwidths))
    self.bandwidth.variable.set("16")
    self.bandwidth.frame.grid(row=0, column=0, sticky='w')

    # Resolution
    rs = Hybrid.Entry(hf, 'Resolution ', 5)
    rs.frame.grid(row=0, column=1, sticky='w')
    self.resolution = rs.variable

    # Cut-off
    co = Hybrid.Entry(hf, 'Cut-off level ', 5)
    co.frame.grid(row=1, column=0, sticky='w')
    self.cutoff = co.variable

    self.save_button = Tkinter.Button(hf, text="Get from Volume Viewer", command=self.get_cutoff)
    self.save_button.grid(row=1, column=1, sticky='w')

    # Space
    msg = Tkinter.Label(hf, anchor='w', justify='left')
    msg.grid(row=2, column=0, sticky='ew')
    row = row + 1
    self.message_label = msg

    # Advanced options
    ostext = 'Advanced options:'
    osh = Tkinter.Label(hf, font="Verdana 10 bold italic", fg="navy", text=ostext)
    osh.grid(row=3, column=0, sticky='w')

    # Fitting criterion
    opt = Hybrid.Radiobutton_Row(hf, 'Fitting criterion ', ('Laplacian', 'Linear'))
    opt.frame.grid(columnspan=2, sticky='w')
    self.fitting_criterion = opt.variable
    self.opt_widget = opt

    # Saved solutions
    sa = Hybrid.Entry(hf, 'Saved solutions ', 5)
    sa.frame.grid(row=5, column=0, sticky='w')
    self.saved_solutions = sa.variable

    # Translational sampling
    ts = Hybrid.Entry(hf, 'Translational sampling ', 5)
    ts.frame.grid(row=5, column=1, sticky='w')
    self.translational_sampling = ts.variable

    # Translational scan strategy
    tss = Hybrid.Radiobutton_Row(hf, 'Translational scan ', ('Masking', 'Full', 'Limited'))
    tss.frame.grid(columnspan=2, sticky='w')
    self.translational_scan = tss.variable
    self.tss_widget = tss

    # Options panel
    rowPanel = row
    op = Hybrid.Popup_Panel(parent)
    opf = op.frame
    opf.grid(row = rowPanel, column = 0, sticky = 'news')
    opf.grid_remove()
    opf.columnconfigure(0, weight=1)
    self.options_panel = op.panel_shown_variable
    row += 1
    orow = 0

    hf = Tkinter.Frame(opf)
    hf.grid(row=0, column=0, sticky='ew')
    hf.columnconfigure(1, weight=1)

    msg = Tkinter.Label(hf, width = 40, anchor = 'w', justify = 'left')
    msg.grid(row=0, column = 0, sticky = 'ew')
    row = row + 1
    self.message_label = msg

    # Options panel close button
    cb = op.make_close_button(hf)
    cb.grid(row=1, column=1, sticky='e')

    osf = Tkinter.Frame(opf)
    osf.grid(row=2, column=0, sticky='w', columnspan=2)

    # Expert options
    advtext = 'Other options for expert users:'
    advh = Tkinter.Label(osf, font="Verdana 10 bold italic", fg="navy", text=advtext)
    advh.grid(row=5, column=0, sticky='w')

    # Peaks explored per docking
    ne = Hybrid.Entry(osf, 'Number peaks explored per docking ', 5)
    ne.frame.grid(row=6, column=0, sticky='w')
    self.peaks_docking = ne.variable

    # Peaks stored per iteration
    np = Hybrid.Entry(osf, 'Number peaks stored per iteration ', 5)
    np.frame.grid(row=7, column=0, sticky='w')
    self.peaks_iteration = np.variable

    # Peaks stored in the search
    nr = Hybrid.Entry(osf, 'NNumber peaks stored in the search ', 5)
    nr.frame.grid(row=8, column=0, sticky='w')
    self.peaks_search = nr.variable

    # Peaks stored in the multi docking search
    nrm = Hybrid.Entry(osf, 'Number peaks stored in the multi docking search ', 5)
    nrm.frame.grid(row=9, column=0, sticky='w')
    self.peaks_multid_search = nrm.variable

    # Translational threshold in grid units
    rt = Hybrid.Entry(osf, 'Translational threshold in grid units ', 5)
    rt.frame.grid(row=10, column=0, sticky='w')
    self.translational_threshold = rt.variable

    # Rotational threshold in degrees
    rc = Hybrid.Entry(osf, 'Rotational threshold in degrees ', 5)
    rc.frame.grid(row=11, column=0, sticky='w')
    self.rotational_threshold = rc.variable

    # Width between spherical layers
    lw = Hybrid.Entry(osf, 'Width between spherical layers ', 5)
    lw.frame.grid(row=12, column=0, sticky='w')
    self.width_spherical = lw.variable

    # Density cut-off of simulated map
    co2 = Hybrid.Entry(osf, 'Density cut-off of simulated map ', 5)
    co2.frame.grid(row=13, column=0, sticky='w')
    self.cutoff_2 = co2.variable

    # Specify a label width so dialog is not resized for long messages.
    msg = Tkinter.Label(parent, width = 40, anchor = 'w', justify = 'left')
    msg.grid(row = row, column = 0, sticky = 'ew')
    row = row + 1
    self.message_label = msg

    # Results panel
    row = rowPanel
    re = Hybrid.Popup_Panel(parent)
    ref = re.frame
    ref.grid(row=rowPanel, column=0, sticky='news')
    ref.grid_remove()
    ref.columnconfigure(0, weight=1)
    self.results_panel = re.panel_shown_variable
    row += 1
    orow = 0

    resf = Tkinter.Frame(ref)
    resf.grid(row=0, column=0, sticky='ew')
    resf.columnconfigure(1, weight=1)

    msgR = Tkinter.Label(resf, width=40, anchor='w', justify='left')
    msgR.grid(row=0, column=0, sticky='ew')
    row = row + 1
    self.messageR_label = msgR

    # Results label
    ostextR = 'Choose a solution to visualize ADP EM results:'
    oshR = Tkinter.Label(resf, font="Verdana 10 bold italic", fg="navy", text=ostextR)
    oshR.grid(row=1, column=0, sticky='w')

    # Results panel close button
    cbR = re.make_close_button(resf)
    cbR.grid(row=1, column=1, sticky='e')

    self.mmf = Tkinter.Frame(ref)
    self.mmf.grid(row=rowPanel, column=0, sticky='w')
    row = row + 1

    # Disable Results panel at first
    self.results_button = self.buttonWidgets['Results']
    self.results_button['state'] = 'disabled'

    # Update internal state of components
    self.update_gray_out()

  # ---------------------------------------------------------------------------
  # Shows a message in ADP EM Plugin
  #
  def message(self, text):

    self.message_label['text'] = text
    self.message_label.update_idletasks()

  # ---------------------------------------------------------------------------
  # Map chosen to fit into base map.
  #
  def fit_map(self):

    m = self.object_menu.getvalue()
    from VolumeViewer import Volume
    if isinstance(m, Volume):
      return m
    
    return None

  # ---------------------------------------------------------------------------
  # Atoms chosen in dialog for fitting.
  #
  def fit_atoms(self):

    m = self.object_menu.getvalue()
    if m == 'selected atoms':
      from chimera import selection
      atoms = selection.currentAtoms()
      return atoms

    from VolumeViewer import Volume
    from chimera import Molecule
    if isinstance(m, Molecule):
      return m.atoms
    
    return []

  # ---------------------------------------------------------------------------
  # Updates the interla state of the parameter when a molecule is chosen
  #
  def object_chosen_cb(self, obj_name):

    self.update_gray_out()

  # ---------------------------------------------------------------------------
  # Updates the internal state of the parameters
  #
  def update_gray_out(self):

    state ='normal'
    for c in self.opt_widget.frame.winfo_children():
      c['state'] = state

  # ---------------------------------------------------------------------------
  # Gets the parameter values introduced by the user
  #
  def get_options_chimera (self):

    # Fitting criterion
    if 'Laplacian' in self.fitting_criterion.get():
      self.adp_em_chimera_laplacian_val = "1"
    else:
      self.adp_em_chimera_laplacian_val = "0"

    # Saved solutions
    if len(self.saved_solutions.get()) > 0:
      self.adp_em_chimera_saved_solutions_val = self.saved_solutions.get()

    # Translational sampling
    if len(self.translational_sampling.get()) > 0:
      self.adp_em_chimera_sampling_val = self.translational_sampling.get()

    # Translational scan
    if 'Masking' in self.translational_scan.get():
      self.adp_em_chimera_scan_val = "2"
    elif 'Limited' in self.translational_scan.get():
      self.adp_em_chimera_scan_val = "1"
    else:
      self.adp_em_chimera_scan_val = "0"

    # Peaks explored per docking
    if len(self.peaks_docking.get()) > 0:
      self.adp_em_chimera_peaks_explored_val = self.peaks_docking.get()

    # Peaks stored per iteration
    if len(self.peaks_iteration.get()) > 0:
      self.adp_em_chimera_peaks_iteration_val = self.peaks_iteration.get()

    # Peaks stored in the search
    if len(self.peaks_search.get()) > 0:
      self.adp_em_chimera_peaks_search_val = self.peaks_search.get()

    # Peaks stored in the multi docking search
    if len(self.peaks_multid_search.get()) > 0:
      self.adp_em_chimera_peaks_msearch_val = self.peaks_multid_search.get()

    # Translational threshold
    if len(self.translational_threshold.get()) > 0:
      self.adp_em_chimera_translational_threshold_val = self.translational_threshold.get()

    # Rotational threshold
    if len(self.rotational_threshold.get()) > 0:
      self.adp_em_chimera_rotational_threshold_val = self.rotational_threshold.get()
    else:
      def_val = 360/int(self.bandwidth.variable.get())
      self.adp_em_chimera_rotational_threshold_val = str(def_val)

    # Width between spherical layers
    if len(self.width_spherical.get()) > 0:
      self.adp_em_chimera_width_layers = self.width_spherical.get()

    # Density cutoff of simulated map
    if len(self.cutoff_2.get()) > 0:
      self.adp_em_chimera_cutoff_simulated_val = self.cutoff_2.get()


  # ---------------------------------------------------------------------------
  # Performs the ADP EM process
  #
  def Fit(self):

    # If a fitting is performed when Results panel is active, close it
    for widget in self.mmf.winfo_children():
      widget.destroy()
      # and clean the array which store the solutions (previous fitting)
      self.solutions_chimera = []
    self.results_panel.set(False)

    # Validation of the parameters introduced by the user
    if self.check_models() is False:
      return

    # Disable Fit, Options and Close buttons when ADP EM process is performed
    self.disable_process_buttons()

    # Retrieve plugin path
    self.plugin_path = __file__[:__file__.index(self.plugin_folder)]

    #-----------------------
    # Calling ADPEM process
    #-----------------------
    from subprocess import STDOUT, PIPE, Popen
    import os, sys

    # Get the full path of ADP EM process
    command = self.plugin_path + self.adp_em
    # Set the workspace
    cwd = self.plugin_path + self.plugin_folder
    # Set the file name that will be generated by ADP EM with all the solutions parameters:
    # center of mass, Euler Angles and traslations of the PDB solutions
    self.filename_chimera = cwd + "chimera.bin"

    # PDB selected in the menu
    pdbSelected = self.object_menu.getvalue()
    # Map selected in the menu
    mapSelected = self.map_menu.volume()
    # Temporal pdb that will be wrote when performing the process.
    # It is necessary because if the user moves the molecule, the camera changes, and the
    # coordinates and the internal states of the pdb and the map remain inconsistent
    # So saving the pdb ensures that the ADP EM process is going to be done with consistent
    #  values relative to the map
    pdb_path = cwd + "temporal_" + pdbSelected.name

    # Save the pdb xform just before performing the fitting
    self.saved_pdb = pdbSelected.openState.xform
    # Save the pdb xform premultiplied by the map inverse xform
    # This is necessary to be able to move the molecule to the origin regardless the ADP EM process
    self.xf = pdbSelected.openState.xform
    self.xf.premultiply(mapSelected.openState.xform.inverse())

    # Back to origin (Ensure that pdb and map have the same internal state when performing ADP EM process
    pdbSelected.openState.xform = mapSelected.openState.xform

    # Save pdb relative to the map
    from Midas import write
    write(pdbSelected, mapSelected, pdb_path)

    # Update user view
    pdbSelected.openState.xform = self.saved_pdb

    # Record position state
    self.record_position_state()

    # Get options values
    self.get_options_chimera()

    # Variable to move to the center of mass when moving the pdb to the first ADP EM solution
    self.fitting_center_mass = True

    # Retrieve the full command to perform the fitting: ap_em + arguments
    cmd = [command, mapSelected.openedAs[0], pdb_path,
           self.bandwidth.variable.get(), self.cutoff.get(), self.resolution.get(),
           self.adp_em_chimera_laplacian, self.adp_em_chimera_laplacian_val,
           self.adp_em_chimera_saved_solutions, self.adp_em_chimera_saved_solutions_val,
           self.adp_em_chimera_sampling, self.adp_em_chimera_sampling_val,
           self.adp_em_chimera_scan, self.adp_em_chimera_scan_val,
           self.adp_em_chimera_peaks_explored, self.adp_em_chimera_peaks_explored_val,
           self.adp_em_chimera_peaks_iteration, self.adp_em_chimera_peaks_iteration_val,
           self.adp_em_chimera_peaks_search, self.adp_em_chimera_peaks_search_val,
           self.adp_em_chimera_peaks_msearch, self.adp_em_chimera_peaks_msearch_val,
           self.adp_em_chimera_translational_threshold, self.adp_em_chimera_translational_threshold_val,
           self.adp_em_chimera_rotational_threshold, self.adp_em_chimera_rotational_threshold_val,
           self.adp_em_chimera_width_layers, self.adp_em_chimera_width_layers_val,
           self.adp_em_chimera_cutoff_simulated, self.adp_em_chimera_cutoff_simulated_val,
           self.adp_em_chimera_no_save, self.adp_em_chimera_opt]

    # Execute the command with the respective arguments creating pipes between the process and Chimera
    # Pipes will be associated to the standard output and standard error required to show the process log
    # in the window
    adp_em_process = Popen(cmd, stdout=PIPE, stderr=PIPE, cwd=cwd, universal_newlines=True)

    # Text widget for process log that will showthe standard output of the process
    from Tkinter import *
    root = Tk()
    root.wm_title("ADP EM Process Log")
    import ttk
    self.var_det = IntVar(root)
    self.pbar_det = ttk.Progressbar(root, orient="horizontal", length=400, mode="determinate", variable=self.var_det, maximum=100)
    self.pbar_det.pack(side=TOP, fill=X)
    S = Scrollbar(root)
    T = Text(root, height=30, width=85)
    S.pack(side=RIGHT, fill=Y)
    T.pack(side=LEFT, fill=Y)
    S.config(command=T.yview)
    T.config(yscrollcommand=S.set)

    # Read first line
    line = adp_em_process.stdout.readline()
    # Variable to check the process status and show its output in a friendly format to the user
    process_progress = False

    # Continue reading the standard output until ADP EM is finished
    # If the current line is an iteration for a model, replace in the widget the last showed
    # If it is a new model or is part of the ADP EM process, inserts the line at the end of the widget
    while line:
        if process_progress is False:
          index_before_last_print = T.index(END)
          T.insert(END, line)
        else:
          T.delete(index_before_last_print + "-1c linestart", index_before_last_print)
          T.insert(END, line)
        T.update()
        T.yview(END)
        self.check_steps_marks(line)
        line = adp_em_process.stdout.readline()
        if '[' in line and '%' in line:
          process_progress = True
        elif '100% Finish!':
          process_progress = False
        sys.stdout.flush()
    T.insert(END, "\n --> ADP EM Process has finished. Check 'Results' button to visualize solutions. <--\n")
    T.update()
    T.yview(END)

    # When ADP EM process is finished, the results are set into the Results panel...
    self.fill_results()
    # and the plugin buttons are enabled again
    self.enable_process_buttons()

    # Remove the temporal pdb
    os.remove(pdb_path)

  # ---------------------------------------------------------------------------
  # Fills the Results panel with the corresponding components
  #
  def fill_results(self):

    # Creates a list with all the solutions generated by ADP EM
    a = range(1, self.adp_number_solutions+1)
    SOLUTIONS = ['ADP EM Solution [%s]' % s for s in a]
    SOLUTIONS.insert(0, '---')

    # Creates a drop-down menu with the previous list and a callback to show the chosen one in Chimera
    from CGLtk import Hybrid
    self.mb = apply(Hybrid.Option_Menu, (self.mmf, 'Choose solution ')  + tuple(SOLUTIONS))
    self.mb.variable.set(SOLUTIONS[0])
    self.mb.frame.grid(row=0, column=0, sticky='w')
    self.mb.add_callback(self.show_adp_solution)

    # Button to copy to the Model Panel the fitted molecule
    import Tkinter
    self.save_button =  Tkinter.Button (self.mmf, text="Save solution", command=self.save_solution_adp)
    self.save_button.grid(row=0, column=1, sticky='w')

    self.init_button = Tkinter.Button(self.mmf, text="Back to initial", command=self.back_initial_adp)
    self.init_button.grid(row=1, column=0, sticky='w')

    # Reading solutions for the first time and allocate them in memory
    if len(self.solutions_chimera) == 0:
      self.read_solutions_chimera()

  # ---------------------------------------------------------------------------
  # Records position state between pdb and map
  #
  def record_position_state(self):

    self.bos = self.map_menu.volume().openState
    bxfinv = self.bos.xform.inverse()
    self.xf = self.object_menu.getvalue().openState.xform
    self.xf.premultiply(bxfinv)

  # ---------------------------------------------------------------------------
  # Checks the process log to update the ADP EM process progress bar
  #
  def check_steps_marks(self, line):

    if any(substring in line for substring in self.steps):
      self.var_det.set(self.var_det.get() + 20)
    elif any(substring in line for substring in self.marks):
      self.var_det.set(self.var_det.get() + 10)
    if 'Saved in adpEM[' in line:
      # Get the number of solutions generated by ADP EM from the process log
      self.get_number_solutions(line)

  # ---------------------------------------------------------------------------
  # Gets the number of solutions generated by ADP EM
  #
  def get_number_solutions(self, line):

    index1 = line.index('-')
    index2 = line.index(']')
    self.adp_solutions = line[index1+1:index2]
    self.adp_number_solutions = int(self.adp_solutions)

  # ------------------------------------------------------------------------------
  # Shows a solution generated by ADP EM from the drop-down in the Results Panel
  #
  def show_adp_solution (self):

    # If an ADP EM solution is chosen...
    if not '---' in self.mb.variable.get():

      # Import adphandle to handle pdb transformations
      import adphandle as adph

      # Get the molecule
      m = self.object_menu.getvalue()

      # Roll back previous state. Set the molecule to the previous state
      self.back_current_state_adp()

      # Go back to origin...
      if self.xform_last_solution is None:
        # if a solution is chosen fot the first time, it is necessary to move the molecule to the center of mass
        self.move_center_mass()
      else:
        # if not, move the molecule to the previous position to the one that was in the last solution chosen
        adph.transform_atom_coordinates_adp(m.atoms, self.xform_last_solution.inverse())

      # -----------------------
      # Get selected solution
      # -----------------------
      # Solution_X  = [Euler Angles(Phi, Theta, Psi), Traslation(X, Y, Z)]
      self.adp_chosed_solution = self.get_adp_chosed_solution(self.mb.variable.get())
      so = self.solutions_chimera[self.adp_chosed_solution - 1]

      # Get the Euler Angles associated with the solution
      ea = map(float, so[len(so) / 2:])

      # Get the traslation associated with the solution
      t = map(float, so[:len(so) / 2])

      # -----------------------
      # Xform of the solution
      # -----------------------
      # First apply rotation, then traslation and finally get the related Xform
      xform_solution = adph.euler_xform_adp(ea, t)

      # Updates the coordinates of the molecule with the xform of the chosen solution
      # This moves it into the position of the chosen ADP EM solution
      adph.transform_atom_coordinates_adp(m.atoms, xform_solution)

      # Save current state because user may move the solution
      self.save_current_state_adp()
      self.xform_last_solution = xform_solution

  # ---------------------------------------------------------------------------
  # Moves the molecule to the center of mas
  #
  def move_center_mass(self):

    # Get the molecule
    m = self.object_menu.getvalue()
    # Import adphandle to handle pdb transformations
    import adphandle as adph

    # To move the center of mass, Euler Angles should be (0, 0, 0)
    ea = map(float, "0 0 0".split())
    # The position (traslation) to the pdb center of mass is given by ADP EM
    t = map(float, self.com)

    # First apply rotation, then traslation and finally get the related Xform
    self.xf_center_mass = adph.euler_xform_adp(ea, t)
    # Updates the coordinates of the molecule with the xform of the center of mass
    adph.transform_atom_coordinates_adp(m.atoms, self.xf_center_mass)

  # ---------------------------------------------------------------------------
  # Saves the current state of the molecule relative to the map
  # It is necessary because user may move the solution
  #
  def save_current_state_adp(self, event = None):

    m = self.object_menu.getvalue()
    v = self.map_menu.volume()
    if m is None or v is None:
      return

    # Applies the inverse of the map xform to the molecule xform to save the current
    # state and movement made
    bxfinvC = v.openState.xform.inverse()
    self.xfC = m.openState.xform
    self.xfC.premultiply(bxfinvC)

  # ---------------------------------------------------------------------------
  # Roll back previous state. Set the molecule to the previous state
  #
  def back_current_state_adp(self, event = None):

    m = self.object_menu.getvalue()
    v = self.map_menu.volume()
    if m is None or v is None:
      return

    # This only applies to the first movement after process is finished
    if self.xfC is None:
      # Back to origin
      m.openState.xform = v.openState.xform
      return

    # Update the molecule xform with xform that stores the last state and movement made
    oxfC = v.openState.xform
    oxfC.multiply(self.xfC)
    m.openState.xform = oxfC

  # ---------------------------------------------------------------------------
  # Gets the size of the file that will be generated by ADP EM
  # with all the solutions parameters (chimera.bin)
  #
  def get_size_chimera_bin(self, fileobject):

    # Move the cursor to the end of the file
    fileobject.seek(0, 2)
    size = fileobject.tell()
    # Move the cursor back to the begin of the file
    fileobject.seek(0, 0)
    return size

  # ---------------------------------------------------------------------------
  # Reads all the solutions parameters from the file that will be generated
  # by ADP EM, chimera.bin (center of mass, Euler Angles and traslations)
  #
  def read_solutions_chimera(self):

    import struct
    import os
    # Get the file size
    f = open(self.filename_chimera, "rb")
    size = self.get_size_chimera_bin(f)
    n = (size / 4) - 3;

    # PDB Center of mass
    self.com = struct.unpack('f' * 3, f.read(4 * 3))
    # Data: Traslation(xyz), 3 Ruler Angles (ZXZ convention)
    num = struct.unpack('f' * n, f.read(4 * n))

    solution = []
    for x in range(0, n):
      if x % 6 == 0 and x > 0 and len(solution) > 0:
        self.solutions_chimera.append(solution)
        solution = []
      data = "{:10.6f}".format(num[x])
      solution.append(data.strip())

    # Store all the solutions data intro the main array
    self.solutions_chimera.append(solution)

    # Remove the temporal chimera.bin file generated
    os.remove(self.filename_chimera)

  # ---------------------------------------------------------------------------
  # Gives the index of the selected solution
  #
  def get_adp_chosed_solution(self, s):

    index1 = s.index('[')
    index2 = s.index(']')
    self.chosen_solution = s[index1 + 1:index2]
    return int(self.chosen_solution)


  # ---------------------------------------------------------------------------
  # Makes a copy to the Model Panel of the fitted molecule
  #
  def save_solution_adp(self):

    # Get opened molecule (the one selected in the menu)
    m = self.object_menu.getvalue()

    # Make copy using the copy_molecule native functionality from Chimera
    from Molecule import copy_molecule
    mc = copy_molecule(m)

    # Set copy name
    mc.name = m.name.split('.')[0] + '_' + self.chosen_solution + '.pdb'

    # Add copy to list of open models
    chimera.openModels.add([mc])

  # ---------------------------------------------------------------------------
  # Moves the molecule to the initial position.
  # Just before the fitting has been made. The idea is to undo all
  # the movements, xforms, and updates made between the original molecule,
  # the solutions and the map
  #
  def back_initial_adp(self):

    # Get map and pdb
    m = self.object_menu.getvalue()
    v = self.map_menu.volume()

    # Roll back previous state. Set the molecule to the previous state
    self.back_current_state_adp()

    # Restore last solution xform and center of mass
    import adphandle as adph
    adph.transform_atom_coordinates_adp(m.atoms, self.xform_last_solution.inverse())
    adph.transform_atom_coordinates_adp(m.atoms, self.xf_center_mass.inverse())

    # Restore all stored xforms and internal states
    self.xfC = None
    self.xform_last_solution = None
    oxf = v.openState.xform
    oxf.multiply(self.xf)
    m.openState.xform = oxf

  # ---------------------------------------------------------------------------
  # Disables the the ADP EM GUI Fit, Options, Close and Results buttons
  #
  def disable_process_buttons(self):
    self.fit_button = self.buttonWidgets['Fit']
    self.fit_button['state'] = 'disabled'
    self.options_button = self.buttonWidgets['Options']
    self.options_button['state'] = 'disabled'
    self.close_ch_button = self.buttonWidgets['Close']
    self.close_ch_button['state'] = 'disabled'
    self.results_button['state'] = 'disabled'

  # ---------------------------------------------------------------------------
  # Enables the the ADP EM GUI Fit, Options, Close and Results buttons
  #
  def enable_process_buttons(self):
    self.fit_button = self.buttonWidgets['Fit']
    self.fit_button['state'] = 'normal'
    self.options_button = self.buttonWidgets['Options']
    self.options_button['state'] = 'normal'
    self.close_ch_button = self.buttonWidgets['Close']
    self.close_ch_button['state'] = 'normal'
    self.results_button['state'] = 'normal'

  # ---------------------------------------------------------------------------
  # When Options button is pressed, the Result panel is hidden and Options
  # panel is shown
  #
  def Options(self):
    self.results_panel.set(False)
    self.options_panel.set(not self.options_panel.get())

  # ---------------------------------------------------------------------------
  # When Results button is pressed, the Options panel is hidden and Results
  # panel is shown
  #
  def Results(self):
    self.options_panel.set(False)
    self.results_panel.set(not self.results_panel.get())

  # -----------------------------------------------------------------------------
  # Gets the cut-off level value from the Volume Viewer dialog
  # Useful when the user does not know an appropiate resolution and plays
  # with the map density in this dialog.
  #
  def get_cutoff (self):
    bmap = self.map_menu.data_region()
    if bmap is None:
      self.message('Choose map.')
      return
    from chimera import dialogs
    vdlg = dialogs.find("volume viewer")
    cutoff_panel = vdlg.thresholds_panel.threshold.get()
    self.cutoff.set(cutoff_panel)
    self.message("")

  # -----------------------------------------------------------------------------
  # Validates the values of the parameteres introduced by the users.
  # Moreover, check if molecule and maps are loaded
  #
  def check_models (self):
    fatoms = self.fit_atoms()
    fmap = self.fit_map()
    bmap = self.map_menu.data_region()
    if (len(fatoms) == 0 and fmap is None) or bmap is None:
      self.message('Choose model and map.')
      return False
    if fmap == bmap:
      self.message('Chosen maps are the same.')
      return False
    if len(self.cutoff.get()) == 0:
      self.message('Cutoff must be defined.')
      return False
    if len(self.resolution.get()) == 0:
      self.message('Resolution must be defined.')
      return False
    if float(self.resolution.get()) < 0 or float(self.resolution.get()) >= 60:
      self.message('Resolution must be between 0 and 59.')
      return False
    self.message("")
    return True


# -----------------------------------------------------------------------------
# Returns a list of molecules from the models opened in Chimera
# to be selectables for the fitting
#
def fit_object_models():

  from VolumeViewer import Volume
  from chimera import openModels as om, Molecule
  mlist = om.list(modelTypes = [Molecule])
  folist = ['selected atoms'] + mlist
  return folist

# -----------------------------------------------------------------------------
# Put 'selected atoms' first, then all molecules, then all volumes.
#
def compare_fit_objects(a, b):

  if a == 'selected atoms':
    return -1
  if b == 'selected atoms':
    return 1
  from VolumeViewer import Volume
  from chimera import Molecule
  if isinstance(a,Molecule) and isinstance(b,Volume):
    return -1
  if isinstance(a,Volume) and isinstance(b,Molecule):
    return 1
  return cmp((a.id, a.subid), (b.id, b.subid))

# -----------------------------------------------------------------------------
# Shows the ADP EM Dialog in Chimera when it is registered
#
def show_adp_em_dialog():

  from chimera import dialogs
  return dialogs.display(ADP_EM_Dialog.name)

# -----------------------------------------------------------------------------
# Registers the ADP EM Dialog in Chimera
#
from chimera import dialogs
dialogs.register(ADP_EM_Dialog.name, ADP_EM_Dialog, replace = True)
