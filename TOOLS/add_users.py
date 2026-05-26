#!/usr/bin/env python3
"""
Add or update users in users.db (SQLite) with hashed passwords.

DEFAULT BEHAVIOR:
  If you run with NO --email and NO --file, it will enter INTERACTIVE PROMPT MODE automatically.

It will also:
  - Create the 'users' table if it doesn't exist.
  - ADD missing columns 'role', 'is_active', 'created_at' if they are not present
    (so you won't hit "no column named role" again).

Supported modes (run order: single -> file -> interactive):
  - Single:      --email alice@example.com --password 'Secret123!'
  - Bulk file:   --file users.csv          (columns: email,password)
                 --file users.txt          ("email,password" OR just "email")
                 --delimiter ','           (default ',')
                 --same-password           (prompt once if password missing)
  - Interactive: --prompt   (or just run with no --email/--file)

Options:
  --db PATH                 Override users.db path (auto-detected by default)
  --role ROLE               Default 'user'
  --inactive                Create/update users as inactive (is_active=0)
  --update-existing         If email exists, update password/role/is_active
  --dry-run                 Print actions only (no writes)
"""

import argparse
import csv
import getpass
import os
import re
import sqlite3
from pathlib import Path
from typing import Optional, Tuple

# Use Werkzeug hasher (same family your Flask app uses)
try:
    from werkzeug.security import generate_password_hash
except Exception as e:
    raise SystemExit("Werkzeug is required (pip install flask or werkzeug).") from e

# ---------- locate repo root robustly ----------
def find_repo_root(start: Path) -> Path:
    for p in [start, *start.parents]:
        if (p / "app.py").exists() or (p / "templates").is_dir() or (p / "requirements.txt").exists():
            return p
    return start  # fallback

HERE = Path(__file__).resolve().parent
ROOT = find_repo_root(HERE)
DEFAULT_DB = ROOT / "users.db"

EMAIL_RE = re.compile(r"^[A-Za-z0-9_.+\-]+@[A-Za-z0-9\-.]+\.[A-Za-z]{2,}$")

def log(msg: str):  # simple logger
    print(msg, flush=True)

def connect(db_path: Path) -> sqlite3.Connection:
    db_path = Path(db_path)
    db_path.parent.mkdir(parents=True, exist_ok=True)
    con = sqlite3.connect(db_path)
    con.execute("PRAGMA journal_mode=WAL;")
    con.execute("PRAGMA synchronous=NORMAL;")
    return con

def table_columns(con: sqlite3.Connection, table: str) -> list[str]:
    rows = con.execute(f'PRAGMA table_info("{table}")').fetchall()
    return [r[1] for r in rows]

def ensure_users_table_and_columns(con: sqlite3.Connection):
    # Create correct schema if missing
    con.execute(
        """
        CREATE TABLE IF NOT EXISTS users (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            email TEXT UNIQUE NOT NULL,
            password_hash TEXT NOT NULL
        );
        """
    )
    con.commit()

    cols = {c for c in table_columns(con, "users")}
    if "role" not in cols:
        con.execute("ALTER TABLE users ADD COLUMN role TEXT NOT NULL DEFAULT 'user';")
    if "is_active" not in cols:
        con.execute("ALTER TABLE users ADD COLUMN is_active INTEGER NOT NULL DEFAULT 1;")
    con.commit()


def valid_email(email: str) -> bool:
    return bool(EMAIL_RE.match(email or ""))

def hash_password(password: str) -> str:
    return generate_password_hash(password)

def upsert_user(
    con: sqlite3.Connection,
    email: str,
    password: Optional[str],
    role: str,
    is_active: int,
    update_existing: bool,
    dry_run: bool,
) -> Tuple[bool, str]:
    """
    Returns (changed, message).
    Will only reference columns that exist (post-migration).
    """
    email = email.strip().lower()
    if not valid_email(email):
        return False, f"✗ Invalid email: {email}"

    cols = set(table_columns(con, "users"))
    has_role = "role" in cols
    has_active = "is_active" in cols

    cur = con.execute("SELECT id FROM users WHERE lower(email)=?", (email,))
    row = cur.fetchone()

    if row:
        if not update_existing:
            return False, f"↷ Exists, skipping: {email} (use --update-existing to modify)"
        sets = []
        params = []
        if password is not None:
            sets.append("password_hash=?")
            params.append(hash_password(password))
        if has_role:
            sets.append("role=?")
            params.append(role)
        if has_active:
            sets.append("is_active=?")
            params.append(int(is_active))
        if not sets:
            return False, f"↷ Nothing to update for {email}"
        params.append(email)

        sql = f"UPDATE users SET {', '.join(sets)} WHERE lower(email)=?"
        if dry_run:
            return False, f"[dry-run] UPDATE {email} ({', '.join(sets)})"
        con.execute(sql, params)
        con.commit()
        return True, f"✓ Updated: {email}"
    else:
        if password is None:
            return False, f"✗ New user {email} needs a password"
        pw_hash = hash_password(password)

        # Build INSERT with only existing columns
        cols_list = ["email", "password_hash"]
        params = [email, pw_hash]
        if has_role:
            cols_list.append("role")
            params.append(role)
        if has_active:
            cols_list.append("is_active")
            params.append(int(is_active))

        placeholders = ", ".join(["?"] * len(cols_list))
        sql = f'INSERT INTO users ({", ".join(cols_list)}) VALUES ({placeholders})'
        if dry_run:
            return False, f"[dry-run] INSERT {email} (cols: {', '.join(cols_list)})"
        con.execute(sql, params)
        con.commit()
        return True, f"✓ Inserted: {email}"

def prompt_one_user(default_role: str, default_active: bool) -> Optional[Tuple[str, str, str, int]]:
    while True:
        email = input("Email: ").strip()
        if not email:
            print("Empty email. Aborting this entry.")
            return None
        if not valid_email(email):
            print("Invalid email format. Try again.")
            continue
        print(f"You entered: {email}")
        if input("Is that correct? [y/N]: ").strip().lower() != "y":
            continue
        pw1 = getpass.getpass("Password: ")
        pw2 = getpass.getpass("Confirm password: ")
        if pw1 != pw2:
            print("Passwords do not match. Try again.")
            continue
        role = input(f"Role [{default_role}]: ").strip() or default_role
        act_in = input(f"Active? (y/N) [{'y' if default_active else 'n'}]: ").strip().lower()
        is_active = default_active if act_in == "" else (act_in.startswith("y"))
        return (email, pw1, role, int(is_active))

def load_file_rows(path: Path, delimiter: str) -> list[tuple[str, Optional[str]]]:
    rows: list[tuple[str, Optional[str]]] = []
    with open(path, "r", newline="", encoding="utf-8") as f:
        sample = f.read(2048)
        f.seek(0)
        if delimiter in sample:
            reader = csv.DictReader(f, delimiter=delimiter)
            if reader.fieldnames and any(fn.lower() == "email" for fn in reader.fieldnames):
                for r in reader:
                    email = (r.get("email") or r.get("Email") or "").strip()
                    pwd = (r.get("password") or r.get("Password") or r.get("pass") or "").strip() or None
                    if email:
                        rows.append((email, pwd))
            else:
                f.seek(0)
                reader2 = csv.reader(f, delimiter=delimiter)
                for parts in reader2:
                    if not parts:
                        continue
                    if len(parts) == 1:
                        rows.append((parts[0].strip(), None))
                    else:
                        rows.append((parts[0].strip(), (parts[1] or "").strip() or None))
        else:
            f.seek(0)
            for line in f:
                line = line.strip()
                if not line:
                    continue
                if "," in line:
                    email, pwd = line.split(",", 1)
                    rows.append((email.strip(), pwd.strip() or None))
                else:
                    rows.append((line, None))
    return rows

def main():
    ap = argparse.ArgumentParser(description="Add/update users in users.db with hashed passwords.")
    ap.add_argument("--db", type=Path, default=DEFAULT_DB, help=f"Path to users.db (default: {DEFAULT_DB})")

    # modes
    ap.add_argument("--prompt", action="store_true", help="Interactive mode (also the default if no --email/--file)")
    ap.add_argument("--email", help="Single user email")
    ap.add_argument("--password", help="Single user password (not recommended—use --prompt)")
    ap.add_argument("--file", type=Path, help="Bulk file (csv/tsv/txt) with emails and optional passwords")
    ap.add_argument("--delimiter", default=",", help="Delimiter for --file (default ',')")
    ap.add_argument("--same-password", action="store_true", help="Prompt once for a password for all rows missing one")

    # behavior
    ap.add_argument("--role", default="user", help="Default role for new/updated users")
    ap.add_argument("--inactive", action="store_true", help="Create/update users as inactive")
    ap.add_argument("--update-existing", action="store_true", help="Update existing users instead of skipping")
    ap.add_argument("--dry-run", action="store_true", help="Print actions without writing")

    args = ap.parse_args()

    con = connect(args.db)
    # Ensure table exists and columns are present (this fixes your error)
    ensure_users_table_and_columns(con)

    changed_total = 0

    # Mode 1: Single user via flags
    if args.email and (args.password or args.update_existing):
        email = args.email.strip()
        password = args.password
        ch, msg = upsert_user(
            con,
            email=email,
            password=password,
            role=args.role,
            is_active=(0 if args.inactive else 1),
            update_existing=args.update_existing,
            dry_run=args.dry_run,
        )
        log(msg)
        changed_total += int(ch)

    # Mode 2: Bulk file
    if args.file:
        rows = load_file_rows(args.file, args.delimiter)
        if not rows:
            log(f"✗ No usable rows found in {args.file}")
        default_password = None
        if args.same_password:
            pw1 = getpass.getpass("Enter password to use for rows without a password: ")
            pw2 = getpass.getpass("Confirm: ")
            if pw1 != pw2:
                log("Passwords do not match; aborting.")
                return
            default_password = pw1

        for email, pwd in rows:
            use_pwd = pwd or default_password
            ch, msg = upsert_user(
                con,
                email=email,
                password=use_pwd,
                role=args.role,
                is_active=(0 if args.inactive else 1),
                update_existing=args.update_existing,
                dry_run=args.dry_run,
            )
            log(msg)
            changed_total += int(ch)

    # Mode 3: Interactive — DEFAULT if neither --email nor --file was specified
    interactive_mode = args.prompt or (not args.email and not args.file)
    if interactive_mode:
        default_role = args.role
        default_active = not args.inactive
        while True:
            packed = prompt_one_user(default_role, default_active)
            if not packed:
                pass
            else:
                email, pw, role, is_active = packed
                ch, msg = upsert_user(
                    con,
                    email=email,
                    password=pw,
                    role=role,
                    is_active=is_active,
                    update_existing=args.update_existing,
                    dry_run=args.dry_run,
                )
                log(msg)
                changed_total += int(ch)

            cont = input("Add another user? [y/N]: ").strip().lower()
            if cont != "y":
                break

    if changed_total and not args.dry_run:
        log(f"✅ Done. Changed rows: {changed_total}")
    elif not changed_total:
        log("↷ No changes (nothing to do or dry-run only).")

if __name__ == "__main__":
    main()
