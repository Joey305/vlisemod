from flask import render_template, redirect, url_for, flash, request
from flask_login import login_user, logout_user, current_user
from . import auth_bp
from .models import get_user_by_email, verify_password
from .forms import LoginForm

@auth_bp.route("/login", methods=["GET", "POST"])
def login():
    if current_user.is_authenticated:
        return redirect(url_for("home"))
    form = LoginForm()
    if form.validate_on_submit():
        user = get_user_by_email(form.email.data.strip().lower())
        if user and verify_password(user, form.password.data):
            login_user(user)
            next_url = request.args.get("next") or url_for("home")
            return redirect(next_url)
        flash("Invalid credentials", "danger")
    return render_template("login.html", form=form)

@auth_bp.route("/logout")
def logout():
    if current_user.is_authenticated:
        logout_user()
    return redirect(url_for("auth.login"))
