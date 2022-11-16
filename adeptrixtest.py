from tkinter import *
from tkinterdnd2 import *
import threading

def drop(event):
    var.set(event.data)

ws = TkinterDnD.Tk()
ws.title('PythonGuides')
ws.geometry('300x200')
ws.config(bg='#fcba03')

var = StringVar()
e_box = Entry(ws, textvar=var, width=80)
e_box.pack(fill=X, padx=10)
e_box.drop_target_register(DND_FILES)
e_box.dnd_bind('<<Drop>>', drop)
i = 0
def test():
    while i == 0:
        if var.get() != '':
            Adeptrix.negcontrolfile = filename
            break
threading.Thread(target=test).start()
ws.mainloop()