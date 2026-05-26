import sqlite3
from pathlib import Path
from werkzeug.security import generate_password_hash, check_password_hash

DB_PATH = Path(__file__).resolve().parent.parent / "users.db"

def init_db():
    DB_PATH.touch(exist_ok=True)
    with sqlite3.connect(DB_PATH) as con:
        con.execute("""
        CREATE TABLE IF NOT EXISTS users(
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            email TEXT UNIQUE NOT NULL,
            password_hash TEXT NOT NULL
        )""")
        con.commit()

class User:
    def __init__(self, id, email, password_hash):
        self.id = id
        self.email = email
        self.password_hash = password_hash

    # Flask-Login required props/methods:
    @property
    def is_authenticated(self): return True
    @property
    def is_active(self): return True
    @property
    def is_anonymous(self): return False
    def get_id(self): return str(self.id)

def get_user_by_email(email):
    with sqlite3.connect(DB_PATH) as con:
        row = con.execute("SELECT id,email,password_hash FROM users WHERE email=?", (email,)).fetchone()
        return User(*row) if row else None

def get_user_by_id(user_id):
    with sqlite3.connect(DB_PATH) as con:
        row = con.execute("SELECT id,email,password_hash FROM users WHERE id=?", (user_id,)).fetchone()
        return User(*row) if row else None

def create_user(email, password):
    hash_ = generate_password_hash(password)
    with sqlite3.connect(DB_PATH) as con:
        con.execute("INSERT INTO users(email,password_hash) VALUES(?,?)", (email, hash_))
        con.commit()

def verify_password(user: User, password: str) -> bool:
    return check_password_hash(user.password_hash, password)
