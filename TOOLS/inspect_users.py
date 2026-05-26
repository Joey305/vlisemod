#!/usr/bin/env python
import sqlite3
from pathlib import Path

db_path = Path(__file__).resolve().parent.parent / "users.db"

if not db_path.exists():
    print(f"✗ DB not found: {db_path}")
    raise SystemExit(1)

print(f"DB: {db_path}")

con = sqlite3.connect(db_path)
cur = con.cursor()

print("\n---- TABLES ----")
for (name,) in cur.execute("SELECT name FROM sqlite_master WHERE type='table' ORDER BY name;"):
    print(name)

print("\n---- users SCHEMA ----")
for (schema,) in cur.execute("SELECT sql FROM sqlite_master WHERE type='table' AND name='users';"):
    print(schema)

print("\n---- users COLUMNS ----")
for col in cur.execute("PRAGMA table_info(users);"):
    print(col)

print("\n---- FIRST 50 ROWS ----")
for row in cur.execute("SELECT * FROM users LIMIT 50;"):
    print(row)

print("\n---- DIAGNOSTIC: id/email/hash_len ----")
for row in cur.execute("""
    SELECT id, typeof(id) AS id_type, email, length(password_hash) AS hash_len
    FROM users
    ORDER BY id
    LIMIT 50;
"""):
    print(row)

con.close()
