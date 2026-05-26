import os
from pathlib import Path
from datetime import datetime, timedelta

def manage_ligand_files_in_all_folders(base_dir):
    # Iterate over all directories that start with "Human papillomavirus"
    for virus_dir in Path(base_dir).glob("Human papillomavirus*"):
        if not virus_dir.is_dir():
            continue
        
        # Check for with_ligands.txt directly within the virus directory
        with_ligands_path = virus_dir / 'with_ligands.txt'
        
        # Check if the with_ligands.txt file exists and is older than 10 minutes
        if with_ligands_path.exists():
            file_age = datetime.now() - datetime.fromtimestamp(with_ligands_path.stat().st_mtime)
            if file_age > timedelta(minutes=1):
                print(f"Removing old {with_ligands_path}")
                with_ligands_path.unlink()
            else:
                print(f"{with_ligands_path} is less than 10 minutes old, not removing.")
                continue
        
        # Find the most recent text file in the directory (excluding the one just deleted)
        text_files = [f for f in virus_dir.glob('with_ligands_*.txt') if f != with_ligands_path]
        if not text_files:
            print(f"No text files found in {virus_dir} after removing the old file.")
            continue
        
        most_recent_file = max(text_files, key=os.path.getmtime)
        
        # Rename the most recent text file to with_ligands.txt
        new_with_ligands_path = virus_dir / 'with_ligands.txt'
        print(f"Renaming {most_recent_file} to {new_with_ligands_path}")
        most_recent_file.rename(new_with_ligands_path)

# Example usage
base_directory = './Database_DATA/Sorted_PDB_Files'
manage_ligand_files_in_all_folders(base_directory)
