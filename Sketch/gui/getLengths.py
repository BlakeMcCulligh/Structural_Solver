import tkinter as tk

def get_text_content():
    """Retrieves the text from the Text widget and prints it."""
    content = text_widget.get("1.0", "end-1c") # "1.0" for first char, "end-1c" to exclude trailing newline
    print("Text content:", content)

# Create the main window
root = tk.Tk()
root.title("Text Input Example")

# Create a Text widget
text_widget = tk.Text(root, height=10, width=40, wrap="word") # wrap="word" for word wrapping
text_widget.pack(pady=10, padx=10)

# Optional: Add some initial text
text_widget.insert(tk.END, "Length: ")

# Optional: Add a button to retrieve the text
get_text_button = tk.Button(root, text="Get Text", command=get_text_content)
get_text_button.pack(pady=5)

# Run the Tkinter event loop
root.mainloop()