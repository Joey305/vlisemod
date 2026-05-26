import sqlite3

# Connect to the database
db_name = 'viral_data.db'
conn = sqlite3.connect(db_name)
cursor = conn.cursor()

# Fetch all table names
cursor.execute("SELECT name FROM sqlite_master WHERE type='table';")
tables = cursor.fetchall()

# Function to fetch and display columns for each table
def display_table_info():
    for table in tables:
        table_name = table[0]
        print(f"Table: {table_name}")
        
        # Fetch column names for the current table
        cursor.execute(f"PRAGMA table_info({table_name});")
        columns = cursor.fetchall()
        
        # Display column headers
        column_names = [column[1] for column in columns]
        print("Headers:", column_names)
        print("-" * 40)

# Display table names and headers
display_table_info()

# Close the connection
conn.close()
