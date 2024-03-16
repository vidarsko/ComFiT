import os
import nbformat
from nbformat.v4 import new_notebook

def clear_outputs_in_notebook(notebook_path):
    with open(notebook_path, 'r', encoding='utf-8') as f:
        nb = nbformat.read(f, as_version=4)
    for cell in nb.cells:
        if 'outputs' in cell:
            cell['outputs'] = []
        if 'execution_count' in cell:
            cell['execution_count'] = None
    with open(notebook_path, 'w', encoding='utf-8') as f:
        nbformat.write(nb, f)

def clear_outputs_in_directory(directory_path):
    for file in os.listdir(directory_path):
        if file.endswith(".ipynb"):
            notebook_path = os.path.join(directory_path, file)
            clear_outputs_in_notebook(notebook_path)
            print(f"Cleared outputs in {file}")

# Replace '/path/to/notebooks' with the path to your folder containing .ipynb files
directory_path = 'tutorial/'
clear_outputs_in_directory(directory_path)