import glob
import re
import random
from pathlib import Path

# Print current path
print(Path().absolute())

# Path to your tests folder
tests_folder_path = 'C:/Users/vidar/Desktop/ComFiT/tests/'

# Generate a new random seed
def new_seed():
    return str(random.randint(10000000, 99999999))

# Regular expression pattern to find np.random.seed() lines
pattern = re.compile(r'np\.random\.seed\(\d{8}\)')

# List all Python files in the tests folder
test_files = glob.glob(f'{tests_folder_path}/*.py')

print(test_files)

for file_path in test_files:
    with open(file_path, 'r') as file:
        content = file.read()
    
    # Replace the seed number
    new_content = pattern.sub(lambda match: f'np.random.seed({new_seed()})', content)
    
    # Write the changes back to the file
    with open(file_path, 'w') as file:
        file.write(new_content)
