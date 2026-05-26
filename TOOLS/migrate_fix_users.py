#!/usr/bin/env python3
import sqlite3
from pathlib import Path
import sys

DB = Path(__file__).resolve().parent.parent / "users.db"

DDL_NEW = """
CREATE TABLE users_new (
  id            INTEGER PRIMARY KEY AUTOINCREMENT,
  email         TEXT NOT NULL UNIQUE,
  password_hash TEXT NOT NULL,
  role          TEXT NOT NULL DEFAULT 'user',
  is_active     INTEGER NOT NULL DEFAULT 1
);
"""

COPY_SQL = """
INSERT INTO users_new (id, email, password_hash, role, is_active)
SELECT
  CASE WHEN typeof(id)='integer' THEN id ELSE NULL END,   -- preserve good ids, autogen for NULL
  lower(email),
  password_hash,
  COALESCE(role, 'user'),
  COALESCE(is_active, 1)
FROM users;
"""

CLEAN_SESSIONS = """
DELETE FROM sessions
WHERE user_id IS NULL OR user_id NOT IN (SELECT id FROM users);
"""

def main():
    if not DB.exists():
        print(f"✗ DB not found: {DB}")
        sys.exit(1)

    con = sqlite3.connect(DB)
    con.isolation_level = None  # autocommit off for BEGIN/COMMIT

    cur = con.cursor()

    def q1(sql, params=()):
        cur.execute(sql, params)
        return cur.fetchone()

    def qall(sql, params=()):
        cur.execute(sql, params)
        return cur.fetchall()

    print(f"DB: {DB}")
    print("Current users schema:")
    for (schema,) in qall("SELECT sql FROM sqlite_master WHERE type='table' AND name='users';"):
        print("  ", schema)

    try:
        cur.execute("BEGIN;")

        # 1) Create new table with correct schema
        cur.execute("DROP TABLE IF EXISTS users_new;")
        cur.execute(DDL_NEW)

        # 2) Copy data (preserve existing integer ids; fix NULL ids)
        cur.execute(COPY_SQL)

        # 3) Replace old table
        cur.execute("DROP TABLE users;")
        cur.execute("ALTER TABLE users_new RENAME TO users;")

        # 4) (Optional) Clean sessions that point to non-existent users
        try:
            cur.execute(CLEAN_SESSIONS)
        except sqlite3.OperationalError:
            # sessions table might not exist; ignore
            pass

        cur.execute("COMMIT;")
    except Exception as e:
        cur.execute("ROLLBACK;")
        raise
    finally:
        # VACUUM to rebuild file & sqlite_sequence
        try:
            cur.execute("VACUUM;")
        except Exception:
            pass

    # Show new schema and a quick diagnostic
    print("\nNew users schema:")
    for (schema,) in qall("SELECT sql FROM sqlite_master WHERE type='table' AND name='users';"):
        print("  ", schema)

    print("\nSample rows (id/email):")
    for row in qall("SELECT id, typeof(id) AS id_type, email FROM users ORDER BY id LIMIT 20;"):
        print("  ", row)

    con.close()
    print("\n✅ Migration complete.")

if __name__ == "__main__":
    main()
