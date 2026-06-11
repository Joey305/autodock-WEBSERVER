AutoDock-Vina PrepServer SEO content pack

Drop-in replacements:
- templates/about.html
- templates/documentation.html
- templates/contact.html
- templates/modules.html

New templates to add:
- templates/docking_problems.html
- templates/open_source.html
- templates/api.html
- templates/troubleshooting.html

Patch snippets:
- app_routes_patch.py
- nav_patch.html

Repository URL used:
https://github.com/Joey305/autodock-WEBSERVER

After copying files:
1. Patch app.py with REPOSITORY_URL, new routes, repository_url context variable, and expanded nav_links.
2. Patch templates/partials/nav.html using nav_patch.html if you want external GitHub in the nav.
3. Run:
   python -c "from app import create_app; app=create_app(); print(app.url_map)"
   gunicorn "app:create_app()" --log-level debug
4. Commit and deploy.
