# Add near the top of app.py with your other constants:
REPOSITORY_URL = "https://github.com/Joey305/autodock-WEBSERVER"
SCHURER_LAB_URL = "https://schurerlab.org"

# Inside inject_site_links(), add repository_url to the returned dictionary:
return {
    "contact_email": SITE_CONTACT_EMAIL,
    "contact_mailto": mailto,
    "repository_url": REPOSITORY_URL,
    "tool_links": TOOL_LINKS,
    "lab_link": LAB_LINK,
    "ecosystem_links": [*TOOL_LINKS, LAB_LINK],
    "nav_links": nav_links,
    "public_mode": current_app.config.get("PUBLIC_MODE", True),
}

# Replace nav_links with this expanded version:
nav_links = [
    {"name": "Home", "endpoint": "home"},
    {"name": "Build", "endpoint": "build"},
    {"name": "Workflow", "endpoint": "workflow"},
    {"name": "About", "endpoint": "about"},
    {"name": "Docs", "endpoint": "documentation"},
    {"name": "GitHub", "url": REPOSITORY_URL, "external": True},
    {"name": "Contact", "endpoint": "contact"},
]

# Add these routes under the existing page routes:
@app.get("/docking-preparation-problems")
def docking_problems():
    return render_template("docking_problems.html")

@app.get("/open-source")
def open_source():
    return render_template("open_source.html")

@app.get("/api")
def api_docs():
    return render_template("api.html")

@app.get("/troubleshooting")
def troubleshooting():
    return render_template("troubleshooting.html")

# Change contact route to pass repository_url:
@app.get("/contact")
def contact():
    return render_template("contact.html", repository_url=REPOSITORY_URL)
