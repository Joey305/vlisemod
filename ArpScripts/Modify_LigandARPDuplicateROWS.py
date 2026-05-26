import sqlite3

# Connect to the database
conn = sqlite3.connect('viral_data.db')
cursor = conn.cursor()

# Query to count duplicate rows
count_query = '''
SELECT COUNT(*)
FROM Ligand_Arp_Diagram
WHERE rowid NOT IN (
    SELECT MIN(rowid)
    FROM Ligand_Arp_Diagram
    GROUP BY virus_name, pdb_id, ligand, chain, residue_id
);
'''

# Execute the count query and fetch the result
cursor.execute(count_query)
duplicate_count = cursor.fetchone()[0]

# Print the number of duplicates found
print(f"Number of duplicate rows found: {duplicate_count}")

# If there are duplicates, proceed to delete them
if duplicate_count > 0:
    # Query to remove duplicate rows
    delete_query = '''
    DELETE FROM Ligand_Arp_Diagram
    WHERE rowid NOT IN (
        SELECT MIN(rowid)
        FROM Ligand_Arp_Diagram
        GROUP BY virus_name, pdb_id, ligand, chain, residue_id
    );
    '''
    
    # Execute the delete query
    cursor.execute(delete_query)
    
    # Commit the changes
    conn.commit()
    print("Duplicate rows removed successfully.")

# Close the connection
conn.close()
