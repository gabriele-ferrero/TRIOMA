import tkinter as tk
from tkinter import simpledialog


class Component:
    def __init__(self, canvas, x, y):
        self.canvas = canvas
        self.x, self.y = x, y
        self.id = canvas.create_rectangle(
            x - 25, y - 25, x + 25, y + 25, fill="blue", tags="draggable"
        )
        self.canvas.tag_bind(self.id, "<Button-1>", self.on_click)
        self.canvas.tag_bind(self.id, "<B1-Motion>", self.on_drag)

    def on_click(self, event):
        new_color = simpledialog.askstring(
            "Edit Component", "Enter new color (e.g., red, blue):"
        )
        if new_color:
            self.canvas.itemconfig(self.id, fill=new_color)

    def on_drag(self, event):
        dx, dy = event.x - self.x, event.y - self.y
        self.canvas.move(self.id, dx, dy)
        self.x, self.y = event.x, event.y


class Port:
    def __init__(self, canvas, x, y):
        self.canvas = canvas
        self.id = canvas.create_oval(x - 5, y - 5, x + 5, y + 5, fill="red")
        self.canvas.tag_bind(self.id, "<Button-1>", self.on_click)

    def on_click(self, event):
        # This method can be expanded to handle port-specific interactions
        print("Port clicked at", event.x, event.y)


class Application(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("Component Connection UI")
        self.canvas = tk.Canvas(self, width=600, height=400, bg="white")
        self.canvas.pack(padx=10, pady=10)

        # Bindings for creating ports
        self.canvas.bind("<Double-Button-1>", self.create_port)

        # Temporary storage for line start point
        self.temp_line_start = None

    def create_port(self, event):
        port = Port(self.canvas, event.x, event.y)
        if not self.temp_line_start:
            self.temp_line_start = (event.x, event.y)
        else:
            x0, y0 = self.temp_line_start
            x1, y1 = event.x, event.y
            self.canvas.create_line(x0, y0, x1, y1, arrow=tk.LAST, fill="green")
            self.temp_line_start = None

    def run(self):
        self.mainloop()


if __name__ == "__main__":
    app = Application()
    app.run()
