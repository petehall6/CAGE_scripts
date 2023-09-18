import os
import tkinter as tk
from tkinter import ttk
from tkinter import filedialog as file_dialog
from tkinter.messagebox import showinfo



#individual scripts for each emailer with names that reflect what they do


#TODO
#possible packer scirpt for easier updating

#TODO
#*** SM Emailer
#*** Mouse Tail Emailer
#*** Design Emailer
#*** Shaina Emailer
#*** PCR 1 Emailer



INITIALS = {
    "PMH" : "Hall",
    "JK" : "Klein",
    "SM": "Miller"
}




def pick_file():
    filetypes = (('all','*.*'),('excel', '*.xls'))
    filename = file_dialog.askopenfilename(title='Open a file', initialdir='Z:\ResearchHome\Groups\millergrp\home\common\Python\_pete\emailer_gui', filetypes=filetypes)
    
    excel_title.insert(tk.END,filename)

def radio_select():
    print(f"Initals selected: {str(var.get())} ")


def parse_excel():
    return

def attach_email():
    return

#create root window
root = tk.Tk()
root.title('Emailer GUI')
root.resizable(True,True)
root.geometry('800x600')
tabControl = ttk.Notebook(root)

tab1 = ttk.Frame(tabControl)
tab2 = ttk.Frame(tabControl)


tabControl.add(tab1, text='Emailer 1')
tabControl.add(tab2, text='Emailer 2')


#radio button for initals

radio_strings = []

for initial in INITIALS:
    radio_strings.append(initial)
    

var = tk.StringVar()
var.set(radio_strings[0])

for i in radio_strings:
    radio_btn1 = tk.Radiobutton(tab1, text = i, variable=var, value=i, command=radio_select)
    radio_btn2 = tk.Radiobutton(tab2, text = i, variable=var, value=i, command=radio_select)
    radio_btn1.pack(anchor='w')
    radio_btn2.pack(anchor='w')



file_pick_btn = ttk.Button(tab1,text='Choose Excel Template', command = pick_file)
other_file_pick_btn = ttk.Button(tab2,text='Choose different Excel Template', command = pick_file)


excel_title = tk.Text(tab1, height=2, width = 25)
excel_other_title = tk.Text(tab2, height=2, width = '10')

#possible to pack with list.
#separate controls based on side, alignment and pack via different list
controls = [file_pick_btn, other_file_pick_btn, excel_title, excel_other_title]


tabControl.pack(expand = 1, fill = 'both')
for con in controls:
    con.pack(side="left")




#run loop
root.mainloop()