import os

# This is to get the directory that the program
# is currently running in.
path = input("Enter the path:\n")
dir_path = os.path.dirname(os.path.realpath(path))
count = 0
for root, dirs, files in os.walk(dir_path):
    for file in files:
        if file.endswith('.fasta'):
            count += 1
print(count)
