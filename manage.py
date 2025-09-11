# manage.py
import getpass
from app import create_app
from models import db, User
from werkzeug.security import generate_password_hash

app = create_app()

def _model_has_column(model, name: str) -> bool:
    try:
        # works for real SQLAlchemy mapped columns
        return name in model.__table__.c
    except Exception:
        # fall back (e.g., hybrid properties)
        return hasattr(model, name)

def create_user():
    with app.app_context():
        email = input("Email: ").strip().lower()
        if not email:
            print("Email required"); return
        if User.query.filter_by(email=email).first():
            print("User already exists"); return

        pw1 = getpass.getpass("Password: ")
        pw2 = getpass.getpass("Confirm password: ")
        if pw1 != pw2 or not pw1:
            print("Passwords do not match / empty"); return

        kwargs = {"email": email, "password_hash": generate_password_hash(pw1)}
        if _model_has_column(User, "is_confirmed"):
            confirm = input("Mark as confirmed? [y/N]: ").strip().lower() == "y"
            kwargs["is_confirmed"] = confirm

        u = User(**kwargs)
        db.session.add(u)
        db.session.commit()
        print(f"Created user {email}")

def list_users():
    with app.app_context():
        rows = User.query.all()
        if not rows:
            print("No users yet."); return
        has_confirmed = _model_has_column(User, "is_confirmed")
        for u in rows:
            status = getattr(u, "is_confirmed", None) if has_confirmed else "-"
            print(u.id, u.email, status)

def delete_user():
    with app.app_context():
        email = input("Email to delete: ").strip().lower()
        u = User.query.filter_by(email=email).first()
        if not u:
            print("No such user"); return
        db.session.delete(u)
        db.session.commit()
        print("Deleted", email)

def edit_user():
    with app.app_context():
        email = input("Email to edit: ").strip().lower()
        u = User.query.filter_by(email=email).first()
        if not u:
            print("No such user"); return

        print(f"Editing {u.email} (id={u.id})")
        print("1) Change password")
        print("2) Change email")
        has_confirmed = _model_has_column(User, "is_confirmed")
        if has_confirmed:
            print("3) Set confirmed/unconfirmed")
        choice = input("Select option [1-3]: ").strip()

        if choice == "1":
            pw1 = getpass.getpass("New password: ")
            pw2 = getpass.getpass("Confirm new password: ")
            if pw1 != pw2 or not pw1:
                print("Passwords do not match / empty"); return
            u.password_hash = generate_password_hash(pw1)

        elif choice == "2":
            new_email = input("New email: ").strip().lower()
            if not new_email:
                print("Email required"); return
            if User.query.filter_by(email=new_email).first():
                print("Another user already has that email"); return
            u.email = new_email

        elif choice == "3" and has_confirmed:
            val = input("Mark as confirmed? [y/N]: ").strip().lower() == "y"
            setattr(u, "is_confirmed", val)
        else:
            print("No changes made."); return

        db.session.commit()
        print("Saved.")

if __name__ == "__main__":
    print("Commands: create-user | list-users | delete-user | edit-user")
    cmd = input("> ").strip()
    if cmd == "create-user":
        create_user()
    elif cmd == "list-users":
        list_users()
    elif cmd == "delete-user":
        delete_user()
    elif cmd == "edit-user":
        edit_user()
    else:
        print("Unknown command")
