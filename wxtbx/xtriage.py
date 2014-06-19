
"""
Classes for displaying Xtriage results using wxPython (and matplotlib).
"""

from __future__ import division
import wxtbx.windows
import wxtbx.tables
import wxtbx.plots
import mmtbx.scaling
import wx.lib.scrolledpanel
import wx

TEXT_WIDTH = 800

class wx_panel (wx.lib.scrolledpanel.ScrolledPanel) :
  def OnChildFocus (self, evt) :
    pass

class wx_output (wxtbx.windows.ChoiceBook,
                 mmtbx.scaling.xtriage_output) :
  """
  Xtriage output handler implemented as a wxPython notebook.
  """
  gui_output = True
  def __init__ (self, *args, **kwds) :
    super(wx_output, self).__init__(*args, **kwds)
    self._current_panel = None
    self._current_sizer = None
    self._graphs = []
    self._in_box = False

  def SetupScrolling (self) :
    """
    Initialize scrollbars on child panels.
    """
    for i_page in range(self.GetPageCount()) :
      panel = self.GetPage(i_page)
      panel.Layout()
      panel.SetupScrolling(scrollToTop=True)

  def show_header (self, title) :
    """
    Creates a new notebook page with the specified title.  This will be the
    output panel for all subsequent method calls until show_header() is called
    again.
    """
    panel = wx_panel(parent=self, style=wx.SUNKEN_BORDER)
    szr = wx.BoxSizer(wx.VERTICAL)
    panel.SetSizer(szr)
    self._current_panel = panel
    self._current_sizer = szr
    self.AddPage(panel, title)

  def show_sub_header (self, title) :
    """
    Create a wx.StaticBox
    """
    assert (self._current_panel is not None)
    if (self._in_box) :
      self._current_sizer = self._current_panel.GetSizer()
    box = wx.StaticBox(parent=self._current_panel,
      label=title,
      style=wx.NO_BORDER)
    box_szr = wx.StaticBoxSizer(box, wx.VERTICAL)
    self._current_sizer.Add(box_szr, 0, wx.ALL|wx.EXPAND, 5)
    self._current_sizer = box_szr
    self._in_box = True

  def show_text (self, text) :
    """
    Create wx.StaticText object(s) with automatic wrapping.  Double-newlines
    will be treated as paragraph breaks, otherwise newlines are replaced by
    spaces.
    """
    assert (self._current_panel is not None)
    blocks = text.split("\n\n")
    for block in blocks :
      block = " ".join([ l.strip() for l in block.splitlines() ]).strip()
      wxtxt = wx.StaticText(parent=self._current_panel,
        label=block)
      wxtxt.Wrap(TEXT_WIDTH)
      self._current_sizer.Add(wxtxt, 0, wx.ALL, 5)

  def show_preformatted_text (self, text) :
    """
    Draw monospaced text, preserving the original formatting.
    """
    assert (self._current_panel is not None)
    if text.startswith("\n") :
      text = text[1:]
    wx_txt = wx.StaticText(parent=self._current_panel,
      label=text)
    font = wx_txt.GetFont()
    # FIXME this seems not to work on wxPython 3/Mac OS 10.9
    font.SetFamily(wx.FONTFAMILY_MODERN)
    wx_txt.SetFont(font)
    self._current_sizer.Add(wx_txt, 0, wx.ALL, 5)

  def warn (self, text) :
    """
    Create wx.StaticText object(s) with automatic wrapping.  Font will be
    boldface with red foreground.
    """
    assert (self._current_panel is not None)
    blocks = text.split("\n\n")
    for block in blocks :
      block = " ".join([ l.strip() for l in block.splitlines() ]).strip()
      wxtxt = wx.StaticText(parent=self._current_panel,
        label=block)
      wxtxt.Wrap(TEXT_WIDTH)
      font = wxtxt.GetFont()
      font.SetWeight(wx.FONTWEIGHT_BOLD)
      wxtxt.SetFont(font)
      wxtxt.SetForegroundColour((240,0,0))
      self._current_sizer.Add(wxtxt, 0, wx.ALL, 5)

  def show_lines (self, text) :
    assert (self._current_panel is not None)
    wxtxt = wx.StaticText(parent=self._current_panel,
      label=text)
    wxtxt.Wrap(TEXT_WIDTH)
    self._current_sizer.Add(wxtxt, 0, wx.ALL, 5)

  def show_paragraph_header (self, text) :
    """
    Draws left-aligned boldface text in a slightly larger font.
    """
    assert (self._current_panel is not None)
    wx_text = wx.StaticText(
      parent=self._current_panel,
      label=text)
    font = wx_text.GetFont()
    font.SetWeight(wx.FONTWEIGHT_BOLD)
    font.SetPointSize(14)
    wx_text.SetFont(font)
    self._current_sizer.Add(wx_text, 0, wx.ALL, 5)

  def show_table (self, table, indent=0, plot_button=False) :
    """
    Draw a wx.ListCtrl and populate from the table.  Can optionally include
    a button to launch a graph viewer window; this is used when the table
    contains more than one graph.
    """
    assert (self._current_panel is not None)
    height = min(400, max(80, 20 + (20*table.n_rows)))
    width = TEXT_WIDTH - 40
    wxtable = wxtbx.tables.TableView(
      parent=self._current_panel,
      style=wx.LC_REPORT|wx.SIMPLE_BORDER,
      size=(width, height))
    wxtable.SetTable(table)
    self._current_sizer.Add(wxtable, 0, wx.ALL, 5)
    if (plot_button) :
      assert hasattr(table, "get_graph")
      btn = wx.Button(parent=self._current_panel, label="Show graphs")
      self._current_sizer.Add(btn, 0, wx.ALL, 5)
      self.Bind(wx.EVT_BUTTON, wxtable.OnViewGraphs, btn)

  def show_plot (self, table) :
    """
    Create an inline matplotlib plot, using the wxtbx.plots.loggraph wrapper.
    """
    assert (self._current_panel is not None)
    graph = wxtbx.plots.iotbx_data_plot_base(
      parent=self._current_panel,
      tables=[table],
      size=(TEXT_WIDTH - 40, min(500, TEXT_WIDTH * 0.75)))
    graph.set_plot(table.only_plot())
    self._graphs.append((graph, table.title))
    self._current_sizer.Add(graph, 0, wx.ALL|wx.EXPAND, 5)

  def show_plots_row (self, tables) :
    pass

  def show_text_columns (self, rows, indent=0) :
    prefix = " "*indent
    n_cols = len(rows[0])
    sizer = wx.FlexGridSizer(rows=len(rows), cols=n_cols)
    self._current_sizer.Add(sizer)
    for row in rows :
      assert (len(row) == n_cols)
      txt1 = wx.StaticText(parent=self._current_panel, label=prefix+row[0])
      font = txt1.GetFont()
      font.SetWeight(wx.FONTWEIGHT_BOLD)
      txt1.SetFont(font)
      sizer.Add(txt1, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
      for cell in row[1:] :
        txt_cell = wx.StaticText(parent=self._current_panel, label=cell)
        sizer.Add(txt_cell, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)

  def newline (self) :
    pass

  def OnSaveImage (self, event) :
    graph_names = [ title for (graph, title) in self._graphs ]
    selected = wxtbx.misc_dialogs.ChoiceSelector(
      title="Save plot image",
      message="Please choose a plot to save:",
      choices=graph_names,
      defaultChoice=0)
    graph_index = graph_names.index(selected)
    (graph, title) = self._graphs[graph_index]
    graph.save_image()

class XtriageFrame (wx.Frame) :
  """
  Frame for displaying a wx_output notebook independently (used in AutoSol
  and similar apps).
  """
  def __init__ (self, *args, **kwds) :
    super(XtriageFrame, self).__init__(*args, **kwds)
    self.output_panel = wx_output(parent=self)
    self._result = None

  def SetResult (self, result) :
    self._result = result
    result.show(out=self.output_panel)
    self.output_panel.SetupScrolling()
    self.output_panel.SetPage(0)