import os
import nbformat
from nbconvert.preprocessors import ExecutePreprocessor
import tempfile
from pathlib import Path
import time

def open_notebook(notebook_path):
    with open(notebook_path, 'r') as f:
        nb = nbformat.read(f, as_version=4)
    return nb

def run_notebook(notebook_path):
    nb = open_notebook(notebook_path)
    ep = ExecutePreprocessor(timeout=600, kernel_name='python3')

    # Create a temporary directory to save the modified notebook
    with tempfile.TemporaryDirectory() as tmpdirname:
        temp_notebook_path = os.path.join(tmpdirname, 'temp_notebook.ipynb')
        with open(temp_notebook_path, 'w') as f:
            nbformat.write(nb, f)
        
        try:
            ep.preprocess(nb, {'metadata': {'path': os.path.dirname(notebook_path)}})
            print(f"Successfully ran {notebook_path}")
        except Exception as e:
            print(f"Error running {notebook_path}: {e}")

def run_all_notebooks(folder_path):
    total_start_time = time.time()
    for filename in os.listdir(folder_path):
        if filename.endswith(".ipynb"):
            notebook_start_time = time.time()
            run_notebook(os.path.join(folder_path, filename))
            print(f"Notebook {filename} took {time.time() - notebook_start_time} seconds to test")
    print(f"Total time taken: {time.time() - total_start_time} seconds")

# Specify the folder containing the notebooks
folder_path = 'tutorial/'

run_all_notebooks(folder_path)
