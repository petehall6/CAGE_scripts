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
from billing_emailer import Billing_Tab
from cell_drop_off_emailer import DropOff_Tab
from status_emailer import Status_Tab
from design_emailer import Design_Tab
from ngs_emailer import NGS_Tab
from tails_emailer import Tails_Tab


app = tbs.Window(
    title="CAGE Emailer",
    themename = "superhero",
    size=(1600,1200),
    resizable=(True,True),
    
)

colors = app.style.colors
note_tab = tbs.Notebook(master=app,
                        height=1100)

main_frame = tbs.Frame(note_tab)


bill_tab = tbs.Frame(note_tab)
drop_tab = tbs.Frame(note_tab)
stats_tab = tbs.Frame(note_tab)
design_tab = tbs.Frame(note_tab)
ngs_tab = tbs.Frame(note_tab)
tails_tab = ScrolledFrame(note_tab)

#call individual emailer scripts to fill in tabs
Billing_Tab(bill_tab)
DropOff_Tab(drop_tab)
Status_Tab(stats_tab)
Design_Tab(design_tab)
NGS_Tab(ngs_tab)
Tails_Tab(tails_tab)

note_tab.add(bill_tab, text="Billing")
note_tab.add(drop_tab, text="Cell Drop Off")
note_tab.add(stats_tab, text="Status")
note_tab.add(design_tab, text="gRNA Designs")
note_tab.add(ngs_tab, text="NGS Analsysis")
note_tab.add(tails_tab.container, text="Tails")
note_tab.pack(pady=20)


note_tab.select(tails_tab.container)
app.mainloop()

