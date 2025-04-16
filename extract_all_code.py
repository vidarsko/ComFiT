import os

# This script aggregates all Python code, documentation, and test files from the 'comfit' directory and its subdirectories.
# It creates three separate output files: one for code, one for documentation, and one for tests.

# --- Code Aggregation ---

# Define the base directory where 'comfit' is located.
# Assumes the script is run from the directory containing 'comfit'.
# If not, provide the full path to the 'comfit' directory.
base_dir = 'comfit'

# List of subdirectories within 'comfit' to process for code
code_subdirs = [
    '', # To include files directly within 'comfit'
    'core',
    'plot',
    'quantum_mechanics',
    'bose_einstein_condensate',
    'nematic_liquid_crystal',
    'phase_field_crystal'
]

# Output file name for code
code_output_filename = 'all_code.txt'

# Check if the base directory exists
if not os.path.isdir(base_dir):
    print(f"Error: Base directory '{base_dir}' not found.")
    print("Please ensure the script is run from the correct location or update the 'base_dir' variable.")
    exit()

try:
    # Open the code output file in write mode with UTF-8 encoding
    with open(code_output_filename, 'w', encoding='utf-8') as outfile:
        print(f"Starting code extraction into {code_output_filename}...")

        # Iterate through the specified subdirectories for code
        for subdir in code_subdirs:
            current_dir_path = os.path.join(base_dir, subdir)

            # Check if the subdirectory exists
            if not os.path.isdir(current_dir_path):
                print(f"Warning: Code directory not found, skipping: {current_dir_path}")
                continue

            print(f"Processing code directory: {current_dir_path}")

            # List all items (files and directories) in the current directory
            try:
                for item_name in os.listdir(current_dir_path):
                    item_path = os.path.join(current_dir_path, item_name)

                    # Check if it's a Python file
                    if os.path.isfile(item_path) and item_name.endswith('.py'):
                        print(f"  Adding code file: {item_path}")
                        try:
                            # Open the source file in read mode
                            with open(item_path, 'r', encoding='utf-8', errors='ignore') as infile:
                                content = infile.read()

                            # Write a header and the content to the output file
                            outfile.write(f"--- Start of file: {item_path} ---\n")
                            outfile.write(content)
                            outfile.write(f"\n--- End of file: {item_path} ---\n\n")
                        except Exception as e:
                            print(f"    Error reading file {item_path}: {e}")
            except Exception as e:
                 print(f"Error listing directory {current_dir_path}: {e}")

    print(f"Successfully aggregated code into {code_output_filename}")

except Exception as e:
    print(f"An error occurred during the code aggregation process: {e}")


# --- Documentation Aggregation ---

# Define the documentation directory (relative to the script's location or base_dir)
# Assuming 'docs' is at the same level as 'comfit' or provide full path
docs_dir = 'docs' # Adjust if 'docs' is inside 'comfit' or elsewhere
doc_output_filename = 'all_documentation.txt'

# Check if the documentation directory exists
if not os.path.isdir(docs_dir):
    print(f"\nWarning: Documentation directory '{docs_dir}' not found. Skipping documentation aggregation.")
else:
    try:
        # Open the documentation output file in write mode with UTF-8 encoding
        with open(doc_output_filename, 'w', encoding='utf-8') as outfile:
            print(f"\nStarting documentation extraction into {doc_output_filename}...")
            print(f"Processing documentation directory: {docs_dir}")

            # Walk through the documentation directory
            for root, dirs, files in os.walk(docs_dir):
                # Exclude specific directories if needed, e.g., build artifacts
                dirs[:] = [d for d in dirs if d not in ['_build', '.ipynb_checkpoints']] # Example exclusion

                for filename in files:
                    # Check if the file is a Markdown file
                    if filename.endswith('.md'):
                        file_path = os.path.join(root, filename)
                        print(f"  Adding documentation file: {file_path}")
                        try:
                            # Open the source file in read mode
                            with open(file_path, 'r', encoding='utf-8', errors='ignore') as infile:
                                content = infile.read()

                            # Write a header and the content to the output file
                            outfile.write(f"--- Start of file: {file_path} ---\n")
                            outfile.write(content)
                            outfile.write(f"\n--- End of file: {file_path} ---\n\n")
                        except Exception as e:
                            print(f"    Error reading file {file_path}: {e}")

        print(f"Successfully aggregated documentation into {doc_output_filename}")

    except Exception as e:
        print(f"An error occurred during the documentation aggregation process: {e}")


# --- Test Aggregation ---

# Define the tests directory (relative to the script's location or base_dir)
tests_dir = 'tests' # Adjust if 'tests' is located elsewhere
test_output_filename = 'all_tests.txt'

# Check if the tests directory exists
if not os.path.isdir(tests_dir):
    print(f"\nWarning: Tests directory '{tests_dir}' not found. Skipping test aggregation.")
else:
    try:
        # Open the test output file in write mode with UTF-8 encoding
        with open(test_output_filename, 'w', encoding='utf-8') as outfile:
            print(f"\nStarting test extraction into {test_output_filename}...")
            print(f"Processing tests directory: {tests_dir}")

            # Walk through the tests directory
            for root, dirs, files in os.walk(tests_dir):
                 # Exclude specific directories if needed, e.g., cache
                dirs[:] = [d for d in dirs if d not in ['__pycache__', '.pytest_cache']] # Example exclusion

                for filename in files:
                    # Check if the file is a Python file
                    if filename.endswith('.py'):
                        file_path = os.path.join(root, filename)
                        print(f"  Adding test file: {file_path}")
                        try:
                            # Open the source file in read mode
                            with open(file_path, 'r', encoding='utf-8', errors='ignore') as infile:
                                content = infile.read()

                            # Write a header and the content to the output file
                            outfile.write(f"--- Start of file: {file_path} ---\n")
                            outfile.write(content)
                            outfile.write(f"\n--- End of file: {file_path} ---\n\n")
                        except Exception as e:
                            print(f"    Error reading file {file_path}: {e}")

        print(f"Successfully aggregated tests into {test_output_filename}")

    except Exception as e:
        print(f"An error occurred during the test aggregation process: {e}")


print("\nScript finished.")