# ==============================
# auth.py
# ==============================
from __future__ import annotations
from flask import Blueprint, redirect, url_for
from flask_login import logout_user

auth_bp = Blueprint("auth", __name__, url_prefix="/auth")

@auth_bp.route("/login", methods=["GET","POST"])
def login():
    return redirect(url_for("home"))

@auth_bp.route("/logout")
def logout():
    logout_user()
    return redirect(url_for("home"))

