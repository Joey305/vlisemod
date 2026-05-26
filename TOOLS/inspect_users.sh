#!/usr/bin/env bash
set -euo pipefail

# Usage:
#   TOOLS/inspect_users.sh                 # uses ./users.db
#   TOOLS/inspect_users.sh /path/to/db.db  # custom path

# Resolve repo root (git) or fall back to script directory
if git_root=$(git rev-parse --show-toplevel 2>/dev/null); then
  ROOT="$git_root"
else
  ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
fi

DB="${1:-$ROOT/users.db}"

if ! command -v sqlite3 >/dev/null 2>&1; then
  echo "✗ sqlite3 not found in PATH" >&2
  exit 1
fi

if [[ ! -f "$DB" ]]; then
  echo "✗ DB not found: $DB" >&2
  exit 1
fi

echo "DB: $DB"
echo "---- TABLES ----"
sqlite3 "$DB" "SELECT name FROM sqlite_master WHERE type='table' ORDER BY name;"

echo
echo "---- users SCHEMA (.schema users) ----"
sqlite3 "$DB" ".schema users"

echo
echo "---- users COLUMNS (PRAGMA table_info) ----"
sqlite3 "$DB" <<'SQL'
.headers on
.mode column
PRAGMA table_info(users);
SQL

echo
echo "---- FIRST 50 ROWS (SELECT *) ----"
sqlite3 "$DB" <<'SQL'
.headers on
.mode column
SELECT * FROM users LIMIT 50;
SQL

echo
