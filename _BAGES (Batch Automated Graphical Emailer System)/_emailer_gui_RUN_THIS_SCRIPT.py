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
from batch_emailer import Batch_Tab

app = tbs.Window(
            title="Batch Emailer",
            themename = "vapor",
            size=(1600,800),
            resizable=(True,True),
        )
app.iconbitmap('misc/cage_icon.ico')
note_tab = tbs.Notebook(master=app,
                        height=1100)

main_frame = tbs.Frame(note_tab)

batch_tab = tbs.Frame(note_tab)

#call individual emailer scripts to fill in tabs
Batch_Tab(batch_tab)


note_tab.add(batch_tab, text="Batch")


note_tab.pack(pady=20)


app.mainloop()

