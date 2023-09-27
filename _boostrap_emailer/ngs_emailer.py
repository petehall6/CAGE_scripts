import ttkbootstrap as tbs
from ttkbootstrap.constants import *
from emailer_functions import *


class NGS_Tab(tbs.Frame):
    
    def __init__(self, master_window):
        super().__init__(master_window, padding=(20,10))
        
        self.frame = tbs.Frame(master_window, width=200, height=200, style=PRIMARY)
        
        self.excel_button = tbs.Button(self.frame, text='NGS Emailer').pack()

        self.excel_name = tbs.Label(self.frame, text="Excel Name", font=(18)).pack()
        self.pack(fill=BOTH, expand=YES)
        self.frame.pack()
        
if __name__ == '__main__':
    root = tk.Tk() # the app window
    main_frame = tk.Frame(root, height=200, width=200, bg='blue', bd=2) # main frame
    Tab(main_frame) # instatiate Tab(), sending main_frame as the parent_widget
    tk.Button(main_frame, text='only if class', command=root.destroy).pack()
    main_frame.pack() # display main frame on window
    tk.mainloop()      
        


