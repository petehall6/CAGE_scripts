"""
1a) MAKE SURE YOU ARE IN THE CAGE ENVIRONMENT conda activate cage
1b) run import packages.py first
2) conda activate base
3) conda activate cage
-should reset the environment?  Might have to restart VS code
"""

import ttkbootstrap as tbs
from ttkbootstrap.constants import *
from billing_emailer import Billing_Tab
from ngs_emailer import NGS_Tab
from cell_drop_off_emailer import DropOff_Tab


app = tbs.Window(
    title="CAGE Emailer",
    themename = "superhero",
    size=(1600,800),
    resizable=(True,True),
)

colors = app.style.colors
note_tab = tbs.Notebook(app)

main_frame = tbs.Frame(note_tab)
button1 = tbs.Button(main_frame, text='test')
button1.pack()

tab1 = tbs.Frame(note_tab)
tab2 = tbs.Frame(note_tab)
tab3 = tbs.Frame(note_tab)

#call individual emailer scripts to fill in tabs
Billing_Tab(tab1)
DropOff_Tab(tab2)
NGS_Tab(tab3)

note_tab.add(tab1, text="Billing")
note_tab.add(tab2, text="Cell Drop Off")
note_tab.add(tab3, text="NGS")
note_tab.pack(pady=20)



app.mainloop()