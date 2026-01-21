// static/script.js

document.addEventListener("DOMContentLoaded", function() {
    const modal = document.getElementById("disclaimer-modal");
    const btn = document.getElementById("accept-btn");

    // Check Local Storage to see if they already accepted
    const hasAccepted = localStorage.getItem("prototype_accepted");

    // Note: We check if modal exists to prevent errors on pages without the modal
    if (modal && !hasAccepted) {
        modal.style.display = "flex";
    }

    if (btn) {
        btn.addEventListener("click", function() {
            modal.style.display = "none";
            // Save to Local Storage
            localStorage.setItem("prototype_accepted", "true");
        });
    }
});