import sqlite3

# Path to your database
db_path = "viral_data.db"

# Connect to the database
conn = sqlite3.connect(db_path)
cursor = conn.cursor()

# SQL query to drop the table
drop_table_query = "DROP TABLE IF EXISTS Arpeggio_Contacts_Data"

try:
    # Execute the query
    cursor.execute(drop_table_query)
    conn.commit()
    print("Table 'Arpeggio_Contacts_Data' has been removed successfully.")
except sqlite3.Error as e:
    print(f"Error occurred: {e}")
finally:
    # Close the database connection
    conn.close()
