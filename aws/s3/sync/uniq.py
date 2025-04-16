# Function to extract the file names without the folder path
def extract_filename(file_path):
    return file_path.split('/')[-1]

# Read the contents of both files
with open('files.txt', 'r') as f:
    files_list_1 = [line.strip() for line in f.readlines()]

with open('morefiles.txt', 'r') as f:
    files_list_2 = [line.strip() for line in f.readlines()]

# Extract only the file names (excluding the folder path)
files_set_1 = set(extract_filename(file) for file in files_list_1)
files_set_2 = set(extract_filename(file) for file in files_list_2)

# Find unique files in both sets
unique_files_1 = files_set_1 - files_set_2  # Files in files.txt but not in morefiles.txt
unique_files_2 = files_set_2 - files_set_1  # Files in morefiles.txt but not in files.txt

# Output the unique files
print("Unique files in files.txt:")
for file in unique_files_1:
    print(file)

print("\nUnique files in morefiles.txt:")
for file in unique_files_2:
    print(file)
