#!/usr/bin/env python3
"""
Delete / deactivate users from users.db (SQLite) — with UNDO support

RUN MODES
---------
• No flags → Interactive Q&A:
  - Offer UNDO last change
  - Offer to list users
  - Select by INDEX (supports 2,4-6) or by EMAIL
  - Choose Delete / Deactivate / Reactivate
  - Optional Backup (default: Yes) ⇒ enables Undo
  - Double confirmation

• Flags (non-interactive):
  Selection:
    --email EMAIL            (repeatable)
    --by-id "1,2,3"
    --file users.csv|txt     (auto-detect 'email' column; CSV/TSV/lines)
    --delimiter ","          (default ',')
    --domain med.miami.edu   (repeatable)

  Action:
    --deactivate             (sets is_active=0)
    --reactivate             (sets is_active=1)
    (default action is DELETE)

  Behavior:
    --list                   list users and exit
    --dry-run                show what would change; no writes
    --yes                    skip prompts
    --backup                 force backup before changing (enables Undo)
    --db PATH                path to users.db (default: repo_root/users.db)

  UNDO / Restore:
    --undo-last              restore users.db from the most recent backup created by this tool
    --list-backups           list available backup files
    --restore PATH           restore from a specific backup path

BACKUPS
-------
Stored under: backups/sqlite/
  e.g. backups/sqlite/users.undo-20250903-140501.sqlite
We also write a marker file pointing to the latest backup:
  backups/sqlite/.last_delete_users_backup
"""

import argparse
import csv
import os
import re
import sqlite3
import sys
from datetime import datetime
from pathlib import Path
from typing import Iterable, List, Dict, Tuple, Optional

# ---------- locate repo root robustly ----------
def find_repo_root(start: Path) -> Path:
    for p in [start, *start.parents]:
        if (p / "app.py").exists() or (p / "templates").is_dir() or (p / "requirements.txt").exists():
            return p
    return start

HERE = Path(__file__).resolve().parent
ROOT = find_repo_root(HERE)
DEFAULT_DB = ROOT / "users.db"
BACKUP_DIR = ROOT / "backups" / "sqlite"
LAST_MARK = BACKUP_DIR / ".last_delete_users_backup"

EMAIL_RE = re.compile(r"^[A-Za-z0-9_.+\-]+@[A-Za-z0-9\-.]+\.[A-Za-z]{2,}$")

def log(msg: str) -> None:
    print(msg, flush=True)

def err(msg: str) -> None:
    print(msg, file=sys.stderr, flush=True)

def connect(db_path: Path) -> sqlite3.Connection:
    db_path = Path(db_path)
    if not db_path.exists():
        err(f"✗ DB not found: {db_path}")
        sys.exit(1)
    con = sqlite3.connect(db_path)
    con.execute("PRAGMA journal_mode=WAL;")
    con.execute("PRAGMA synchronous=NORMAL;")
    return con

def table_exists(con: sqlite3.Connection, name: str) -> bool:
    row = con.execute(
        "SELECT 1 FROM sqlite_master WHERE type='table' AND name=?",
        (name,)
    ).fetchone()
    return bool(row)

def table_columns(con: sqlite3.Connection, table: str) -> List[str]:
    rows = con.execute(f'PRAGMA table_info("{table}")').fetchall()
    return [r[1] for r in rows]

def fetch_users(con: sqlite3.Connection) -> List[dict]:
    cols = set(table_columns(con, "users"))
    has_role = "role" in cols
    has_active = "is_active" in cols
    has_created = "created_at" in cols

    select_cols = ["id", "email"]
    if has_role:   select_cols.append("role")
    if has_active: select_cols.append("is_active")
    if has_created: select_cols.append("created_at")

    sql = f"SELECT {', '.join(select_cols)} FROM users ORDER BY id"
    cur = con.execute(sql)
    rows = cur.fetchall()

    return [dict(zip(select_cols, r)) for r in rows]

def print_users(users: List[dict], show_limit: Optional[int] = None) -> None:
    if show_limit is not None:
        users = users[:show_limit]
    if not users:
        log("(no users)")
        return

    keys = ["id", "email"]
    if any("role" in u for u in users): keys.append("role")
    if any("is_active" in u for u in users): keys.append("is_active")
    if any("created_at" in u for u in users): keys.append("created_at")

    header = ["#"] + keys
    widths = [max(len(str(i+1)) for i in range(len(users)))]
    for k in keys:
        widths.append(max(len(k), max(len(str(u.get(k, ""))) for u in users)))
    fmt = "  " + "  ".join(f"{{:{w}}}" for w in widths)

    log(fmt.format(*header))
    for idx, u in enumerate(users, start=1):
        row = [idx] + [u.get(k, "") for k in keys]
        log(fmt.format(*[str(x) for x in row]))

def parse_id_csv(csv_like: str) -> List[int]:
    out: List[int] = []
    for part in csv_like.split(","):
        part = part.strip()
        if not part: continue
        try:
            out.append(int(part))
        except ValueError:
            err(f"✗ Bad id: {part!r} (skipping)")
    return out

def parse_index_input(spec: str, max_index: int) -> List[int]:
    result: List[int] = []
    for chunk in spec.split(","):
        chunk = chunk.strip()
        if not chunk:
            continue
        if "-" in chunk:
            a, b = chunk.split("-", 1)
            try:
                start = int(a.strip()); end = int(b.strip())
            except ValueError:
                continue
            if start > end: start, end = end, start
            for i in range(start, end+1):
                if 1 <= i <= max_index:
                    result.append(i)
        else:
            try:
                i = int(chunk)
                if 1 <= i <= max_index:
                    result.append(i)
            except ValueError:
                continue
    seen = set()
    out = []
    for i in result:
        if i not in seen:
            seen.add(i); out.append(i)
    return out

def load_file_emails(path: Path, delimiter: str) -> List[str]:
    emails: List[str] = []
    with open(path, "r", newline="", encoding="utf-8") as f:
        sample = f.read(2048)
        f.seek(0)
        if delimiter in sample:
            reader = csv.DictReader(f, delimiter=delimiter)
            if reader.fieldnames and any(fn.lower() == "email" for fn in reader.fieldnames):
                for r in reader:
                    e = (r.get("email") or r.get("Email") or "").strip().lower()
                    if e: emails.append(e)
            else:
                f.seek(0)
                reader2 = csv.reader(f, delimiter=delimiter)
                for parts in reader2:
                    if not parts: continue
                    e = (parts[0] or "").strip().lower()
                    if e: emails.append(e)
        else:
            for line in f:
                line = line.strip()
                if not line: continue
                if "," in line:
                    e = line.split(",", 1)[0].strip().lower()
                    emails.append(e)
                else:
                    emails.append(line.lower())
    return emails

def gather_targets(
    con: sqlite3.Connection,
    emails: Iterable[str],
    ids: Iterable[int],
    domains: Iterable[str]
) -> Dict[int, str]:
    targets: Dict[int, str] = {}

    email_list = [e.strip().lower() for e in emails if e and EMAIL_RE.match(e.strip().lower())]
    if email_list:
        q_marks = ",".join("?" for _ in email_list)
        sql = f"SELECT id, email FROM users WHERE lower(email) IN ({q_marks})"
        for uid, em in con.execute(sql, email_list):
            targets[int(uid)] = em

    id_list = [int(i) for i in ids if isinstance(i, int)]
    if id_list:
        q_marks = ",".join("?" for _ in id_list)
        sql = f"SELECT id, email FROM users WHERE id IN ({q_marks})"
        for uid, em in con.execute(sql, id_list):
            targets[int(uid)] = em

    doms = [d.strip().lower().lstrip("@") for d in domains if d and d.strip()]
    for d in doms:
        like = f"%@{d}"
        sql = "SELECT id, email FROM users WHERE lower(email) LIKE ?"
        for uid, em in con.execute(sql, (like,)):
            targets[int(uid)] = em

    return targets

# -------------------- Backups / Undo --------------------
def ensure_backup_dir() -> None:
    BACKUP_DIR.mkdir(parents=True, exist_ok=True)

def backup_db(db_path: Path, mark_for_undo: bool = False) -> Path:
    ensure_backup_dir()
    ts = datetime.now().strftime("%Y%m%d-%H%M%S")
    out = BACKUP_DIR / f"{db_path.stem}.undo-{ts}{db_path.suffix}"
    out.write_bytes(Path(db_path).read_bytes())
    log(f"🗃  Backup written: {out}")
    if mark_for_undo:
        LAST_MARK.write_text(str(out), encoding="utf-8")
    return out

def list_backups() -> List[Path]:
    if not BACKUP_DIR.exists():
        return []
    return sorted([p for p in BACKUP_DIR.glob("*.undo-*.sqlite") if p.is_file()], key=lambda p: p.stat().st_mtime, reverse=True)

def restore_backup(backup_path: Path, db_path: Path) -> None:
    backup_path = Path(backup_path)
    db_path = Path(db_path)
    if not backup_path.exists():
        err(f"✗ Backup not found: {backup_path}")
        sys.exit(1)
    # safety backup of current DB
    ensure_backup_dir()
    ts = datetime.now().strftime("%Y%m%d-%H%M%S")
    safety = BACKUP_DIR / f"{db_path.stem}.safety-{ts}{db_path.suffix}"
    if db_path.exists():
        safety.write_bytes(db_path.read_bytes())
        log(f"🛟 Current DB saved as: {safety}")
    db_path.write_bytes(backup_path.read_bytes())
    log(f"✅ Restored {db_path} from {backup_path}")

# -------------------- Core change --------------------
def delete_or_toggle_users(
    con: sqlite3.Connection,
    user_ids: List[int],
    do_delete: bool,
    deactivate: Optional[bool],
    dry_run: bool
) -> Tuple[int, int]:
    if not user_ids:
        return (0, 0)

    if dry_run:
        action = "DELETE" if do_delete else ("DEACTIVATE" if deactivate else "REACTIVATE")
        log(f"[dry-run] Would {action} users: {user_ids}")
        return (0, 0)

    cur = con.cursor()
    cur.execute("BEGIN;")
    affected_sessions = 0
    affected_users = 0
    try:
        if do_delete and table_exists(con, "sessions"):
            q_marks = ",".join("?" for _ in user_ids)
            cur.execute(f"DELETE FROM sessions WHERE user_id IN ({q_marks})", user_ids)
            affected_sessions = cur.rowcount or 0

        if do_delete:
            q_marks = ",".join("?" for _ in user_ids)
            cur.execute(f"DELETE FROM users WHERE id IN ({q_marks})", user_ids)
            affected_users = cur.rowcount or 0
        else:
            cols = set(table_columns(con, "users"))
            if "is_active" not in cols:
                raise sqlite3.OperationalError("users.is_active column not found; cannot toggle activation.")
            val = 0 if deactivate else 1
            q_marks = ",".join("?" for _ in user_ids)
            cur.execute(f"UPDATE users SET is_active=? WHERE id IN ({q_marks})", [val, *user_ids])
            affected_users = cur.rowcount or 0

        cur.execute("COMMIT;")
    except Exception:
        cur.execute("ROLLBACK;")
        raise

    return (affected_users, affected_sessions)

# -------------------- Interactive Q&A --------------------
def interactive_flow(db_path: Path) -> None:
    # Offer UNDO first
    undo = input("Undo last deletion/deactivation? [y/N]: ").strip().lower()
    if undo == "y":
        if LAST_MARK.exists():
            p = Path(LAST_MARK.read_text(encoding="utf-8").strip())
            if p.exists():
                yn = input(f"Restore from {p.name}? [y/N]: ").strip().lower()
                if yn == "y":
                    restore_backup(p, db_path)
                    return
            else:
                log("Previous undo marker points to a missing file.")
        backups = list_backups()
        if not backups:
            log("No backups available to restore.")
        else:
            log("\nAvailable backups (newest first):")
            for i, b in enumerate(backups, 1):
                log(f"  {i}. {b.name}")
            spec = input("Choose a number to restore (or Enter to cancel): ").strip()
            if spec:
                try:
                    idx = int(spec)
                    if 1 <= idx <= len(backups):
                        restore_backup(backups[idx-1], db_path)
                        return
                except ValueError:
                    pass
            log("Restore canceled.\n")

    con = connect(db_path)
    users = fetch_users(con)

    show = input("Would you like to see the list of users? [Y/n]: ").strip().lower()
    if show in ("", "y", "yes"):
        print_users(users)

    while True:
        pick = input(
            "\nSelect users by (1) index/range or (2) email?\n"
            "Enter 1 or 2 [1]: "
        ).strip()
        if pick == "" or pick == "1":
            if not users:
                log("No users to select.")
                return
            spec = input("Enter indexes (e.g., 2,4-6): ").strip()
            idxs = parse_index_input(spec, max_index=len(users))
            if not idxs:
                log("No valid indexes chosen.")
                continue
            targets = {users[i-1]["id"]: users[i-1]["email"] for i in idxs}
            break
        elif pick == "2":
            emails_raw = input("Enter email(s) separated by commas: ").strip()
            emails = [e.strip().lower() for e in emails_raw.split(",") if e.strip()]
            targets = gather_targets(con, emails=emails, ids=[], domains=[])
            if not targets:
                log("No matching emails found.")
                continue
            break
        else:
            log("Please enter 1 or 2.")

    while True:
        action = input(
            "\nAction? (d)elete, (x) deactivate, (r) reactivate [d]: "
        ).strip().lower()
        if action in ("", "d", "delete"):
            do_delete = True; deactivate = None; verb = "DELETE"
            break
        elif action in ("x", "deactivate"):
            do_delete = False; deactivate = True; verb = "DEACTIVATE"
            break
        elif action in ("r", "reactivate"):
            do_delete = False; deactivate = False; verb = "REACTIVATE"
            break
        else:
            log("Choose d / x / r.")

    # Backup default Yes (enables Undo)
    do_backup = input("Backup DB before changes? [Y/n]: ").strip().lower()
    if do_backup in ("", "y", "yes"):
        backup_db(db_path, mark_for_undo=True)

    log("\nTargets:")
    for uid, em in sorted(targets.items()):
        log(f"  id={uid}\t{em}")

    confirm = input(f"\nAre you sure you want to {verb} {len(targets)} user(s)? [y/N]: ").strip().lower()
    if confirm != "y":
        log("Aborted.")
        return

    second = input(f"Type '{verb}' to confirm: ").strip().upper()
    if second != verb:
        log("Double-confirmation failed. Aborted.")
        return

    affected_users, affected_sessions = delete_or_toggle_users(
        con,
        sorted(targets.keys()),
        do_delete=do_delete,
        deactivate=deactivate,
        dry_run=False
    )

    if do_delete:
        log(f"\n✅ Deleted users: {affected_users}")
        if table_exists(con, "sessions"):
            log(f"🧹 Removed sessions: {affected_sessions}")
    elif deactivate:
        log(f"\n✅ Deactivated users: {affected_users}")
    else:
        log(f"\n✅ Reactivated users: {affected_users}")

# -------------------- CLI --------------------
def main():
    ap = argparse.ArgumentParser(description="Delete / deactivate users from users.db (with undo)")
    ap.add_argument("--db", type=Path, default=DEFAULT_DB, help=f"Path to users.db (default: {DEFAULT_DB})")

    # selection
    ap.add_argument("--email", action="append", help="Email to target (repeatable)")
    ap.add_argument("--by-id", help="Comma-separated user IDs")
    ap.add_argument("--file", type=Path, help="File of emails (csv/tsv/txt)")
    ap.add_argument("--delimiter", default=",", help="Delimiter for --file (default ',')")
    ap.add_argument("--domain", action="append", help="Target emails ending in this domain (repeatable)")

    # actions
    mode = ap.add_mutually_exclusive_group()
    mode.add_argument("--deactivate", action="store_true", help="Set is_active=0 instead of deleting")
    mode.add_argument("--reactivate", action="store_true", help="Set is_active=1 instead of deleting")
    # default action is DELETE

    # behavior
    ap.add_argument("--list", action="store_true", help="List users and exit")
    ap.add_argument("--dry-run", action="store_true", help="Show what would change")
    ap.add_argument("--yes", action="store_true", help="No confirmation")
    ap.add_argument("--backup", action="store_true", help="Backup DB before changing (enables Undo)")
    ap.add_argument("--interactive", action="store_true", help="Force interactive Q&A")

    # undo / restore
    ap.add_argument("--undo-last", action="store_true", help="Restore DB from most recent undo backup created by this tool")
    ap.add_argument("--list-backups", action="store_true", help="List available backups")
    ap.add_argument("--restore", type=Path, help="Restore from a specific backup path")

    args = ap.parse_args()

    # Undo / restore operations take priority
    if args.list_backups:
        backups = list_backups()
        if not backups:
            log("No backups found.")
            return
        log("Available backups (newest first):")
        for b in backups:
            log(f"  {b}")
        return

    if args.undo_last:
        if LAST_MARK.exists():
            p = Path(LAST_MARK.read_text(encoding="utf-8").strip())
            if p.exists():
                yn = "y" if args.yes else input(f"Restore from {p}? [y/N]: ").strip().lower()
                if yn == "y":
                    restore_backup(p, args.db)
                else:
                    log("Aborted.")
                return
        # fallback to newest
        backups = list_backups()
        if not backups:
            err("✗ No undo backups available.")
            sys.exit(1)
        p = backups[0]
        yn = "y" if args.yes else input(f"Restore from newest {p}? [y/N]: ").strip().lower()
        if yn == "y":
            restore_backup(p, args.db)
        else:
            log("Aborted.")
        return

    if args.restore:
        p = Path(args.restore)
        yn = "y" if args.yes else input(f"Restore from {p}? [y/N]: ").strip().lower()
        if yn == "y":
            restore_backup(p, args.db)
        else:
            log("Aborted.")
        return

    # Decide mode: interactive by default when NO selection flags and no --list
    no_targets = not any([args.email, args.by_id, args.file, args.domain])
    if args.interactive or (no_targets and not args.list):
        interactive_flow(args.db)
        return

    # Non-interactive paths
    con = connect(args.db)

    if args.list:
        print_users(fetch_users(con))
        return

    emails: List[str] = []
    ids: List[int] = []
    domains: List[str] = []

    if args.email:
        emails.extend([e.strip().lower() for e in args.email if e.strip()])
    if args.by_id:
        ids.extend(parse_id_csv(args.by_id))
    if args.file:
        emails.extend(load_file_emails(args.file, args.delimiter))
    if args.domain:
        domains.extend(args.domain)

    targets = gather_targets(con, emails=emails, ids=ids, domains=domains)
    if not targets:
        log("↷ No matching users found.")
        return

    do_delete = not (args.deactivate or args.reactivate)
    deactivate = True if args.deactivate else (False if args.reactivate else None)

    log("\nTargets:")
    for uid, em in sorted(targets.items()):
        log(f"  id={uid}\t{em}")

    # Ensure we can UNDO if requested
    if not args.dry_run and (args.backup or do_delete or deactivate is not None):
        # If user asked --backup, or any mutating op, create a backup if --backup set
        if args.backup:
            backup_db(args.db, mark_for_undo=True)

    if not args.yes:
        action = "DELETE" if do_delete else ("DEACTIVATE" if deactivate else "REACTIVATE")
        yn = input(f"\nAbout to {action} {len(targets)} user(s). Proceed? [y/N]: ").strip().lower()
        if yn != "y":
            log("Aborted.")
            return

    affected_users, affected_sessions = delete_or_toggle_users(
        con, sorted(targets.keys()), do_delete, deactivate, args.dry_run
    )

    if args.dry_run:
        log("\n[dry-run] No changes made.")
    else:
        if do_delete:
            log(f"\n✅ Deleted users: {affected_users}")
            if table_exists(con, "sessions"):
                log(f"🧹 Removed sessions: {affected_sessions}")
        elif deactivate:
            log(f"\n✅ Deactivated users: {affected_users}")
        else:
            log(f"\n✅ Reactivated users: {affected_users}")

if __name__ == "__main__":
    main()
