import tkinter as tk
from tkinter import ttk

def on_select(event):
    item_id = tree.selection()[0]
    column = tree.identify_column(event.x)
    
    print(item_id)
    print(column)
    
    values = ["Option 1", "Option 2", "Option 3"]
    if column == "#1":  # Check if the click is in the first column
        x, y, width, height = tree.bbox(item_id, column)
        combo = ttk.Combobox(root, values=values)
        combo.place(x=x, y=y, width=width, height=height)
        combo.current(0)  # Set the initial value
        combo.bind("<<ComboboxSelected>>", lambda e: update_cell(item_id, column, combo.get()))
        combo.focus_set()

def update_cell(item_id, column, value):
    tree.item(item_id, values=(value, tree.item(item_id, "values")[1]))

root = tk.Tk()
tree = ttk.Treeview(root, columns=("Column 1", "Column 2"))
tree.heading("Column 1", text="Column 1")
tree.heading("Column 2", text="Column 2")
tree.pack()

tree.insert("", "end", values=("Item 1", "Value 1"))
tree.insert("", "end", values=("Item 2", "Value 2"))

tree.bind("<ButtonRelease-1>", on_select)  # Bind the click event

root.mainloop()