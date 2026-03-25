#!/usr/bin/env python3
import getpass
from app import create_app
from models import db, User

def prompt_email():
    while True:
        e1 = input("Enter email: ").strip().lower()
        e2 = input("Re-enter email: ").strip().lower()

        if not e1 or not e2:
            print("Email cannot be blank.\n")
            continue

        if e1 != e2:
            print("Emails do not match.\n")
            continue

        if "@" not in e1 or "." not in e1.split("@")[-1]:
            print("That does not look like a valid email.\n")
            continue

        return e1

def prompt_password():
    while True:
        p1 = getpass.getpass("Enter password: ")
        p2 = getpass.getpass("Re-enter password: ")

        if not p1 or not p2:
            print("Password cannot be blank.\n")
            continue

        if p1 != p2:
            print("Passwords do not match.\n")
            continue

        if len(p1) < 8:
            print("Password must be at least 8 characters.\n")
            continue

        return p1

def main():
    app = create_app()

    with app.app_context():
        email = prompt_email()
        password = prompt_password()

        user = User.query.filter_by(email=email).first()
        if user is None:
            user = User(email=email)
            db.session.add(user)
            action = "created"
        else:
            action = "updated"

        user.set_password(password)
        db.session.commit()

        # immediate verification
        fresh_user = User.query.filter_by(email=email).first()
        ok = fresh_user.check_password(password)

        print(f"\nUser {action} successfully: {email}")
        print(f"Password verification immediately after save: {ok}\n")

if __name__ == "__main__":
    main()