# ---------------------------------------------------------------------------------
# Dialog for fitting molecules or maps in maps through ADP_EM exhaustive algorithm.
#

import chimera
from chimera.baseDialog import ModelessDialog
import ADP

# ---------------------------------------------------------------------------------
#
class ADP_EM_Dialog(ModelessDialog):

  title = 'ADP EM Exhaustive Fitting'
  name = 'ADP EM'
  buttons = ('Fit', 'Options', 'Results', 'Close')
  help = ('adp_em.html', ADP)
  plugin_folder = 'ADP/'
  adp_em = plugin_folder + 'adp_em'
  steps = ["Step 1", "Step 2", "Step 3"]
  marks = ["ADP_EM", "Trans limit", "Interpolation", "Total Time"]
  mb = None
  solutions_chimera = []
  bos = None
  xf = None
  xfC = None
  xform_last_solution = None

  #-----------------------
  # ADP EM Chimera Commands

  # ADP in Chimera indicator
  adp_em_chimera_opt = "--chimera"

  # No solutions
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
  def fillInUI(self, parent):

    #self.requested_halt = False
    self.xform_handler = None
    self.last_relative_xform = None

    self.max_steps = 2000
    self.ijk_step_size_min = 0.01
    self.ijk_step_size_max = 0.5
    self.last_status_time = 0
    self.status_interval = 0.5    # seconds

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

    #
    # Specify a label width so dialog is not resized for long messages.
    #
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

    ostextR = 'Choose a solution to visualize ADP EM results:'
    oshR = Tkinter.Label(resf, font="Verdana 10 bold italic", fg="navy", text=ostextR)
    oshR.grid(row=1, column=0, sticky='w')

    cbR = re.make_close_button(resf)
    cbR.grid(row=1, column=1, sticky='e')

    self.mmf = Tkinter.Frame(ref)
    self.mmf.grid(row=rowPanel, column=0, sticky='w')
    row = row + 1

    self.results_button = self.buttonWidgets['Results']
    self.results_button['state'] = 'disabled'

    self.update_gray_out()

  # ---------------------------------------------------------------------------
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
  #
  def object_chosen_cb(self, obj_name):

    self.update_gray_out()

  # ---------------------------------------------------------------------------
  #
  def update_gray_out(self):

    state ='normal'
    for c in self.opt_widget.frame.winfo_children():
      c['state'] = state

  # ---------------------------------------------------------------------------
  #
  def fitting_atoms(self):

    return not (self.fit_map())
      

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
  #
  def Fit(self):

    for widget in self.mmf.winfo_children():
      widget.destroy()
      self.solutions_chimera = []
    self.results_panel.set(False)
    self.results_button['state'] = 'disabled'
    if self.check_models() is False:
      return

    self.disable_process_buttons()

    # Retrieve plugin path
    self.plugin_path = __file__[:__file__.index(self.plugin_folder)]

    #-----------------------
    # Calling ADP_EM process
    #-----------------------
    from subprocess import STDOUT, PIPE, Popen
    import os, sys
    command = self.plugin_path + self.adp_em
    cwd = self.plugin_path + self.plugin_folder
    self.filename_chimera = cwd + "chimera.bin"

    # PDB selected in the menu
    pdbSelected = self.object_menu.getvalue()
    # Map selected in the menu
    mapSelected = self.map_menu.volume()
    # Temporal pdb
    pdb_path = cwd + "temporal_" + pdbSelected.name

    self.saved_pdb = pdbSelected.openState.xform
    #self.saved_map = mapSelected.openState.xform

    self.xf = pdbSelected.openState.xform
    self.xf.premultiply(mapSelected.openState.xform.inverse())

    # Back to origin
    pdbSelected.openState.xform = mapSelected.openState.xform

    # Save pdb relative to the map
    from Midas import write
    write(pdbSelected, mapSelected, pdb_path)

    # Update user view
    pdbSelected.openState.xform = self.saved_pdb

    # Record position State
    self.record_position_state()

    self.get_options_chimera()
    self.fitting_center_mass = True

    cmd = [command, mapSelected.openedAs[0], pdb_path, self.bandwidth.variable.get(), self.cutoff.get(), self.resolution.get(),
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

    adp_em_process = Popen(cmd, stdout=PIPE, stderr=PIPE, cwd=cwd, universal_newlines=True)

    # Text widget for process log
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

    line = adp_em_process.stdout.readline()
    process_progress = False
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
    self.enable_process_buttons()
    self.fill_results()
    os.remove(pdb_path)

  # ---------------------------------------------------------------------------
  #
  def fill_results(self):

    a = range(1, self.adp_number_solutions+1)
    SOLUTIONS = ['ADP EM Solution [%s]' % s for s in a]
    SOLUTIONS.insert(0, '---')

    from CGLtk import Hybrid
    self.mb = apply(Hybrid.Option_Menu, (self.mmf, 'Choose solution ')  + tuple(SOLUTIONS))
    self.mb.variable.set(SOLUTIONS[0])
    self.mb.frame.grid(row=0, column=0, sticky='w')
    self.mb.add_callback(self.show_adp_solution)

    import Tkinter
    self.save_button =  Tkinter.Button (self.mmf, text="Save solution", command=self.save_solution_adp)
    self.save_button.grid(row=0, column=1, sticky='w')

   # self.init_button = Tkinter.Button(self.mmf, text="Back to initial", command=self.back_initial_adp)
   # self.init_button.grid(row=1, column=0, sticky='w')

    self.results_button['state'] = 'normal'

    # Reading solutions for the first time
    if len(self.solutions_chimera) == 0:
      self.read_solutions_chimera()

  # ---------------------------------------------------------------------------
  # Record position state
  #
  def record_position_state(self):
    self.bos = self.map_menu.volume().openState
    bxfinv = self.bos.xform.inverse()
    self.xf = self.object_menu.getvalue().openState.xform
    self.xf.premultiply(bxfinv)

  # ---------------------------------------------------------------------------
  #
  def check_steps_marks(self, line):
    if any(substring in line for substring in self.steps):
      self.var_det.set(self.var_det.get() + 20)
    elif any(substring in line for substring in self.marks):
      self.var_det.set(self.var_det.get() + 10)
    if 'Saved in adpEM[' in line:
      self.get_number_solutions(line)

  # ---------------------------------------------------------------------------
  #
  def get_number_solutions(self, line):
    index1 = line.index('-')
    index2 = line.index(']')
    self.adp_solutions = line[index1+1:index2]
    self.adp_number_solutions = int(self.adp_solutions)

  # ---------------------------------------------------------------------------
  #
  def show_adp_solution (self):
    if not '---' in self.mb.variable.get():

      import adphandle as adph
      m = self.object_menu.getvalue()

      # Roll back previous state
      self.back_current_state_adp()

      # Go back to origin
      if self.xform_last_solution is None:
        self.move_center_mass()
      else:
        adph.transform_atom_coordinates_adp(m.atoms, self.xform_last_solution.inverse())

      # Get selected solution
      self.adp_chosed_solution = self.get_adp_chosed_solution(self.mb.variable.get())
      so = self.solutions_chimera[self.adp_chosed_solution - 1]
      ea = map(float, so[len(so) / 2:])
      t = map(float, so[:len(so) / 2])

      # Xform solution
      xform_solution = adph.euler_xform_adp(ea, t)
      adph.transform_atom_coordinates_adp(m.atoms, xform_solution)

      # Save current state because user may move the solution and xform solution
      self.save_current_state_adp()
      self.xform_last_solution = xform_solution

  # ---------------------------------------------------------------------------
  #
  def move_center_mass(self):
    m = self.object_menu.getvalue()
    import adphandle as adph
    ea = map(float, "0 0 0".split())
    t = map(float, self.com)
    self.xf_center_mass = adph.euler_xform_adp(ea, t)
    adph.transform_atom_coordinates_adp(m.atoms, self.xf_center_mass)
    #m.openState.xform = xf
    #self.record_xform_adp(m, xf)

  # ---------------------------------------------------------------------------
  #
  def restore_center_mass(self):
    m = self.object_menu.getvalue()
    import adphandle as adph
    ea = map(float, "0 0 0".split())
    t = map(float, self.com)
    t = [x * -1 for x in t]
    xf = adph.euler_xform_adp(ea, t)
    adph.transform_atom_coordinates_adp(m.atoms, xf)
    # self.record_xform_adp(m, xf)

  # ---------------------------------------------------------------------------
  #
  def save_current_state_adp(self, event = None):

    # m = self.object_menu.getvalue()
    # if m is None or not hasattr(m, 'applied_xform'):
    #   return
    #
    # xf = m.applied_xform.inverse()
    # import adphandle as adph
    # adph.transform_atom_coordinates_adp(m.atoms, xf)
    #
    # self.record_xform_adp(m, None)

    # m = self.object_menu.getvalue()
    # if m is None:
    #   return
    #
    # import adphandle as adph
    # adph.transform_atom_coordinates_adp(m.atoms, self.reset_xform.inverse())
    #
    # if self.fitting_center_mass is False:
    #   self.restore_center_mass()
    # else:
    #   self.fitting_center_mass = False

    m = self.object_menu.getvalue()
    v = self.map_menu.volume()
    if m is None or v is None:
      return

    bxfinvC = v.openState.xform.inverse()
    self.xfC = m.openState.xform
    self.xfC.premultiply(bxfinvC)

  # ---------------------------------------------------------------------------
  #
  def back_current_state_adp(self, event = None):

    # m = self.object_menu.getvalue()
    # if m is None or not hasattr(m, 'applied_xform'):
    #   return
    #
    # xf = m.applied_xform.inverse()
    # import adphandle as adph
    # adph.transform_atom_coordinates_adp(m.atoms, xf)
    #
    # self.record_xform_adp(m, None)

    # m = self.object_menu.getvalue()
    # if m is None:
    #   return
    #
    # import adphandle as adph
    # adph.transform_atom_coordinates_adp(m.atoms, self.reset_xform.inverse())
    #
    # if self.fitting_center_mass is False:
    #   self.restore_center_mass()
    # else:
    #   self.fitting_center_mass = False

    m = self.object_menu.getvalue()
    v = self.map_menu.volume()
    if m is None or v is None:
      return

    if self.xfC is None:
      # Back to origin
      m.openState.xform = v.openState.xform
      return

    oxfC = v.openState.xform
    oxfC.multiply(self.xfC)
    m.openState.xform = oxfC

  # ---------------------------------------------------------------------------
  #
  def record_xform_adp(self, m, xf):

    mxf = getattr(m, 'applied_xform', None)
    if xf is None:
      if mxf:
        delattr(m, 'applied_xform')
    elif mxf:
      mxf.premultiply(xf)
      m.applied_xform = mxf
    else:
      m.applied_xform = xf

  # ---------------------------------------------------------------------------
  #
  def get_size_chimera_bin(self, fileobject):
    fileobject.seek(0, 2)  # move the cursor to the end of the file
    size = fileobject.tell()
    fileobject.seek(0, 0)  # move the cursor back to the begin of the file
    return size

  # ---------------------------------------------------------------------------
  #
  def read_solutions_chimera(self):
    import struct
    import os
    f = open(self.filename_chimera, "rb")
    size = self.get_size_chimera_bin(f)
    n = (size / 4) - 3;

    # PDB CoM
    self.com = struct.unpack('f' * 3, f.read(4 * 3))
    # Data: 3 CoM (xyz), 3 euler angles (ZXZ convention)
    num = struct.unpack('f' * n, f.read(4 * n))

    solution = []
    for x in range(0, n):
      if x % 6 == 0 and x > 0 and len(solution) > 0:
        self.solutions_chimera.append(solution)
        solution = []
      data = "{:10.6f}".format(num[x])
      solution.append(data.strip())
    self.solutions_chimera.append(solution)

    os.remove(self.filename_chimera)

  # ---------------------------------------------------------------------------
  #
  def get_adp_chosed_solution(self, s):
    index1 = s.index('[')
    index2 = s.index(']')
    self.chosed_solution = s[index1 + 1:index2]
    return int(self.chosed_solution)


  # ---------------------------------------------------------------------------
  #
  def save_solution_adp(self):

    import chimera
    m = self.object_menu.getvalue()

    # Make copy
    from Molecule import copy_molecule
    mc = copy_molecule(m)

    # Set copy name
    mc.name = m.name.split('.')[0] + '_' + self.chosed_solution + '.pdb'

    # # Change atom coordinates
    # xf = chimera.Xform.translation(chimera.Vector(20, 0, 0))
    # from MoleculeTransform import transform_atom_coordinates
    # transform_atom_coordinates(mc.atoms, xf)

    # Add copy to list of open models
    chimera.openModels.add([mc])

    # Make sure molecule copy gets same position as molecule.
    # By default it gets same position as model with lowest id number.
    #mc.openState.xform = m.openState.xform

  # ---------------------------------------------------------------------------
  #
  def back_initial_adp(self):

    m = self.object_menu.getvalue()
    v = self.map_menu.volume()

    # Roll back previous state
    self.back_current_state_adp()

    # Restore center of mass
    import adphandle as adph
    #adph.transform_atom_coordinates_adp(m.atoms, self.xform_last_solution.inverse())
    adph.transform_atom_coordinates_adp(m.atoms, self.xf_center_mass.inverse())
    #
    # self.xfC = None
    # self.xform_last_solution = chimera.Xform.xform(1.0, 0, 0, 0, 0, 1.0, 0, 0, 0, 0, 1.0, 0, orthogonalize = 1)
    # m.openState.xform = self.saved_pdb
    #
    # xfv = self.map_menu.volume().openState.xform
    # xfv.multiply()
    # mypdb.openState.xform = oxf
    #
    # bxfinv = mymap.openState.xform.inverse()
    # xf = mypdb.openState.xform
    # xf.premultiply(bxfinv)
    self.xfC = None
    self.xform_last_solution = None
    oxf = v.openState.xform
    oxf.multiply(self.xf)
    m.openState.xform = oxf

  # ---------------------------------------------------------------------------
  #
  def disable_process_buttons(self):
    self.fit_button = self.buttonWidgets['Fit']
    self.fit_button['state'] = 'disabled'
    self.options_button = self.buttonWidgets['Options']
    self.options_button['state'] = 'disabled'
    self.close_ch_button = self.buttonWidgets['Close']
    self.close_ch_button['state'] = 'disabled'

  # ---------------------------------------------------------------------------
  #
  def enable_process_buttons(self):
    self.fit_button = self.buttonWidgets['Fit']
    self.fit_button['state'] = 'normal'
    self.options_button = self.buttonWidgets['Options']
    self.options_button['state'] = 'normal'
    self.close_ch_button = self.buttonWidgets['Close']
    self.close_ch_button['state'] = 'normal'

  # ---------------------------------------------------------------------------
  #
  def Options(self):
    self.results_panel.set(False)
    self.options_panel.set(not self.options_panel.get())


  # ---------------------------------------------------------------------------
  #
  def Results(self):
    self.options_panel.set(False)
    self.results_panel.set(not self.results_panel.get())

  # -----------------------------------------------------------------------------
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
#
def fit_adp_em_dialog(create = False):

  from chimera import dialogs
  return dialogs.find(ADP_EM_Dialog.name, create=create)

# -----------------------------------------------------------------------------
#
def show_adp_em_dialog():

  from chimera import dialogs
  return dialogs.display(ADP_EM_Dialog.name)

# -----------------------------------------------------------------------------
#
def region_name_cb(self, event):
    name = self.region_name.get()
    self.show_named_region(name)
# -----------------------------------------------------------------------------
#
from chimera import dialogs
dialogs.register(ADP_EM_Dialog.name, ADP_EM_Dialog, replace = True)
