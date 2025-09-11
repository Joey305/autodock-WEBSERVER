#!/usr/bin/env python3
"""
Interactive user creator for the AutoDock Prep portal.

Usage:
  python scripts/create_user.py                 # create a new user
  python scripts/create_user.py --update        # update password for existing user
  python scripts/create_user.py --confirm=no    # mark as *not* confirmed
"""

import sys, re, getpass, argparse
from datetime import datetime

# Adjust these imports if your app structure differs
try:
    from app import create_app
    from models import db, User
except Exception as e:
    print("Could not import app/models. Make sure you're running from the project root.")
    raise

from werkzeug.security import generate_password_hash

EMAIL_RE = re.compile(r"^[^@\s]+@[^@\s]+\.[^@\s]+$")

def prompt_email() -> str:
    while True:
        email = input("Email: ").strip()
        if not EMAIL_RE.match(email):
            print("  -> Please enter a valid email address.")
            continue
        return email

def prompt_confirm_flag(default_yes: bool = True) -> bool:
    default_str = "Y/n" if default_yes else "y/N"
    while True:
        ans = input(f"Mark email as confirmed? [{default_str}]: ").strip().lower()
        if ans == "" and default_yes: return True
        if ans == "" and not default_yes: return False
        if ans in ("y", "yes"): return True
        if ans in ("n", "no"): return False
        print("  -> Please answer y or n.")

def prompt_password_twice() -> str:
    while True:
        p1 = getpass.getpass("Password: ")
        if len(p1) < 8:
            print("  -> Password must be at least 8 characters.")
            continue
        p2 = getpass.getpass("Repeat password: ")
        if p1 != p2:
            print("  -> Passwords do not match; try again.")
            continue
        return p1

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--update", action="store_true", help="Update password for an existing user")
    parser.add_argument("--confirm", choices=["yes","no"], default="yes",
                        help="Set initial email-confirmed state (default: yes)")
    args = parser.parse_args()

    app = create_app()
    with app.app_context():
        email = prompt_email()
        user = db.session.query(User).filter_by(email=email).first()

        if args.update:
            if not user:
                print("User not found; use without --update to create.")
                sys.exit(1)
            pwd = prompt_password_twice()
            user.password_hash = generate_password_hash(pwd)
            db.session.commit()
            print(f"Updated password for {email}")
            return

        if user:
            print("A user with this email already exists. Use --update to reset password.")
            sys.exit(1)

        # If your User model does not have `is_confirmed`, this will be ignored below.
        is_confirmed = (args.confirm == "yes") and prompt_confirm_flag(default_yes=True) \
                       or (args.confirm == "no" and prompt_confirm_flag(default_yes=False))

        pwd = prompt_password_twice()

        # Build user
        u = User(
            email=email,
            password_hash=generate_password_hash(pwd),
        )

        # Optional fields if present on your model
        for field, value in (
            ("is_confirmed", is_confirmed),
            ("created_at", datetime.utcnow()),
            ("updated_at", datetime.utcnow()),
        ):
            if hasattr(u, field):
                setattr(u, field, value)

        db.session.add(u)
        db.session.commit()
        print(f"Created user {email} (confirmed={getattr(u, 'is_confirmed', 'n/a')})")

if __name__ == "__main__":
    main()
