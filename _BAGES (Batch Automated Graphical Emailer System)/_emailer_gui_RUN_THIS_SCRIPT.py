"""
1a) MAKE SURE YOU ARE IN THE CAGE ENVIRONMENT conda activate cage
1b) run import packages.py first
2) conda activate base
3) conda activate cage
-should reset the environment?  Might have to restart VS code
"""

import ttkbootstrap as tbs
from ttkbootstrap.constants import *
from ttkbootstrap.scrolled import ScrolledFrame
import sys
import os
from initial_emailer import Initial_Tab
""" from cell_drop_off_emailer import DropOff_Tab
from status_emailer_srm import Status_Tab_srm
from design_emailer import Design_Tab
from ngs_emailer import NGS_Tab
from tails_emailer import Tails_Tab
from ngs_emailer_sample_num import NGS_Tab_Sample_Num """
import emailer_functions


app = tbs.Window(
            title="Batch Emailer",
            themename = "vapor",
            size=(1600,800),
            resizable=(True,True),
        )
app.iconbitmap('misc/cage_icon.ico')
colors = app.style.colors
note_tab = tbs.Notebook(master=app,
                        height=1100)

main_frame = tbs.Frame(note_tab)

init_tab = tbs.Frame(note_tab)
""" drop_tab = tbs.Frame(note_tab)
stats_tab_srm = tbs.Frame(note_tab)
design_tab = tbs.Frame(note_tab)
ngs_tab = tbs.Frame(note_tab)
tails_tab = ScrolledFrame(note_tab)
ngs_sample_tab = tbs.Frame(note_tab) """

#call individual emailer scripts to fill in tabs
Initial_Tab(init_tab)
""" DropOff_Tab(drop_tab)
Status_Tab_srm(stats_tab_srm)
Design_Tab(design_tab)
NGS_Tab(ngs_tab)
Tails_Tab(tails_tab)
NGS_Tab_Sample_Num(ngs_sample_tab) """

note_tab.add(init_tab, text="Initial")
#note_tab.add(drop_tab, text="gRNA Validation")
#note_tab.add(stats_tab_srm, text="Donor Validation")
#ote_tab.add(design_tab, text="Mouse gRNA/donor")

note_tab.pack(pady=20)


app.mainloop()

