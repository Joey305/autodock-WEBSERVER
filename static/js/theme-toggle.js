(() => {
  const STORAGE_KEY = "autodock-prepserver-theme";
  const THEMES = new Set(["light", "dark"]);

  function getInitialTheme() {
    try {
      const saved = localStorage.getItem(STORAGE_KEY);
      if (THEMES.has(saved)) return saved;
    } catch (error) {
      // Storage can fail in restricted browser contexts; fall back gracefully.
    }
    return window.matchMedia && window.matchMedia("(prefers-color-scheme: dark)").matches ? "dark" : "light";
  }

  function applyTheme(theme, persist = false) {
    const normalized = THEMES.has(theme) ? theme : "light";
    const nextTheme = normalized === "dark" ? "light" : "dark";
    document.documentElement.dataset.theme = normalized;

    if (persist) {
      try {
        localStorage.setItem(STORAGE_KEY, normalized);
      } catch (error) {
        // Keep the active in-memory theme even if persistence is unavailable.
      }
    }

    document.querySelectorAll("[data-theme-toggle]").forEach((button) => {
      const icon = button.querySelector("[data-theme-icon]");
      const label = button.querySelector("[data-theme-label]");
      button.setAttribute("aria-pressed", normalized === "dark" ? "true" : "false");
      button.setAttribute("aria-label", `Toggle theme. Switch to ${nextTheme} mode.`);
      button.dataset.currentTheme = normalized;
      if (icon) {
        icon.className = `bi ${normalized === "dark" ? "bi-sun" : "bi-moon-stars"}`;
      }
      if (label) {
        label.textContent = normalized === "dark" ? "Light" : "Dark";
      }
    });
  }

  function initThemeToggle() {
    applyTheme(document.documentElement.dataset.theme || getInitialTheme(), false);
    document.querySelectorAll("[data-theme-toggle]").forEach((button) => {
      button.addEventListener("click", () => {
        const current = THEMES.has(document.documentElement.dataset.theme)
          ? document.documentElement.dataset.theme
          : getInitialTheme();
        applyTheme(current === "dark" ? "light" : "dark", true);
      });
    });
  }

  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", initThemeToggle, { once: true });
  } else {
    initThemeToggle();
  }
})();
