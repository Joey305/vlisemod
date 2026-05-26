"""
This Python script connects to the 'viral_data.db' SQLite database and removes duplicate rows from all tables.
It performs the following steps:
1. Connects to the specified SQLite database.
2. Retrieves the names of all tables in the database.
3. Iterates through each table, removing any duplicate rows based on all column values.
   - Duplicate rows are identified by comparing all columns and keeping the row with the lowest 'rowid'.
4. Commits the changes and closes the database connection.

Note:
- This script assumes that all tables in the database contain a 'rowid' column (which SQLite creates by default unless explicitly disabled).
- It is recommended to back up the database before running this script in case you need to restore the original data.
"""

import sqlite3



def remove_duplicates_from_table(cursor, table_name):
    # Fetch all columns from the table
    cursor.execute(f"PRAGMA table_info({table_name})")
    columns = cursor.fetchall()
    
    # Get the column names
    column_names = [col[1] for col in columns]
    
    # Formulate a query to delete duplicate rows based on all columns
    # We use rowid to preserve one unique row and delete duplicates
    column_str = ', '.join(column_names)
    delete_query = f"""
    DELETE FROM {table_name}
    WHERE rowid NOT IN (
        SELECT MIN(rowid)
        FROM {table_name}
        GROUP BY {column_str}
    )
    """
    
    # Execute the delete query
    cursor.execute(delete_query)

def clean_database(db_path):
    try:
        # Connect to the database
        conn = sqlite3.connect(db_path)
        cursor = conn.cursor()

        # Get all table names
        cursor.execute("SELECT name FROM sqlite_master WHERE type='table'")
        tables = cursor.fetchall()

        for table in tables:
            table_name = table[0]
            print(f"Cleaning table: {table_name}")
            remove_duplicates_from_table(cursor, table_name)

        # Commit the changes
        conn.commit()

        print("Duplicate rows removed from all tables.")

    except sqlite3.Error as e:
        print(f"An error occurred: {e}")

    finally:
        if conn:
            conn.close()

# Path to your viral_data.db
db_path = 'viral_data.db'

# Call the clean_database function
clean_database(db_path)
