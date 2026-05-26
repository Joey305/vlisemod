(function () {
    const STORAGE_KEY = "vlisemod_theme";

    function normalizeTheme(theme) {
        return theme === "dark" ? "dark" : "light";
    }

    function getSavedTheme() {
        try {
            return normalizeTheme(localStorage.getItem(STORAGE_KEY) || "light");
        } catch (e) {
            return "light";
        }
    }

    function saveTheme(theme) {
        try {
            localStorage.setItem(STORAGE_KEY, normalizeTheme(theme));
        } catch (e) {
            // Theme persistence is a convenience; the UI should still toggle if storage is unavailable.
        }
    }

    function updateToggle(theme) {
        const toggle = document.getElementById("theme-toggle");
        if (!toggle) return;

        const isDark = normalizeTheme(theme) === "dark";
        const icon = toggle.querySelector(".theme-toggle-icon");
        const label = toggle.querySelector(".theme-toggle-label");

        toggle.setAttribute("aria-pressed", isDark ? "true" : "false");
        toggle.setAttribute("aria-label", isDark ? "Switch to light mode" : "Switch to dark mode");

        if (icon) icon.textContent = isDark ? "☀️" : "🌙";
        if (label) label.textContent = isDark ? "Light" : "Dark";
    }

    function applyTheme(theme, persist) {
        const normalized = normalizeTheme(theme);
        document.documentElement.setAttribute("data-theme", normalized);
        updateToggle(normalized);
        if (persist) saveTheme(normalized);
    }

    function initThemeToggle() {
        applyTheme(getSavedTheme(), false);

        const toggle = document.getElementById("theme-toggle");
        if (!toggle) return;

        toggle.addEventListener("click", function () {
            const current = normalizeTheme(document.documentElement.getAttribute("data-theme"));
            applyTheme(current === "dark" ? "light" : "dark", true);
        });
    }

    if (document.readyState === "loading") {
        document.addEventListener("DOMContentLoaded", initThemeToggle);
    } else {
        initThemeToggle();
    }
})();
