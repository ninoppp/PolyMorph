import os

def generate_expected_files():
    ns = [2, 4, 6, 8, 10, 13, 17, 20, 25, 30, 35, 40]
    ms = range(100)  # This generates numbers from 0 to 99
    expected_files = [f"{n}_{m}.off" for n in ns for m in ms]
    return expected_files

def find_missing_files(directory):
    expected_files = set(generate_expected_files())
    actual_files = set(os.listdir(directory))
    
    # print number of expected and actual files
    print(f"Expected files: {len(expected_files)}")
    print(f"Actual files: {len(actual_files)}")

    missing_files = expected_files - actual_files
    if missing_files:
        print("Missing files:")
        for file in sorted(missing_files):
            print(file)
    else:
        print("No files are missing.")

# Change 'your_directory_path' to the path of the directory containing the files.
directory_path = 'tissues_varwidth'
find_missing_files(directory_path)
