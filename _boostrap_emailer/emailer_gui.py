import ttkbootstrap as tbs
from ttkbootstrap.constants import *
from billing_emailer import Billing_Tab
from ngs_emailer import NGS_Tab
from other_emailer import Other_Tab


app = tbs.Window(
    title="CAGE Emailer",
    themename = "superhero",
    size=(1400,800),
    resizable=(False,False),
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
NGS_Tab(tab2)
Other_Tab(tab3)


note_tab.add(tab1, text="Billing")
note_tab.add(tab2, text="NGS")
note_tab.add(tab3, text="Other")
note_tab.pack(pady=20)



app.mainloop()