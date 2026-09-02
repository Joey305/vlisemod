$(document).ready(function() {
    // Initialize Select2 for searchable dropdowns
    $('.searchable-select').select2();

    // Modal popup handling
    setupModal();

    // Fetch initial data
    fetchInitialData();

    // Set up event handlers for selections
    setupEventHandlers();

    // Animate PyMOL button if required
    checkForPyMOLAnimation();
});

function setupModal() {
    var modal = $('#betaModal');
    $('.close').on('click', function() {
        modal.hide();
    });

    $(window).on('click', function(event) {
        if ($(event.target).is(modal)) {
            modal.hide();
        }
    });
}

function fetchInitialData() {
    fetch('/get_viruses')
        .then(response => response.json())
        .then(data => {
            const virusSelect = $('#virus');
            data.viruses.forEach(virus => {
                virusSelect.append(new Option(virus, virus));
            });
        });
}

function setupEventHandlers() {
    $('#virus').on('change', handleVirusChange);
    $('#pdb_code').on('change', handlePDBCodeChange);
    $('#ligand').on('change', handleLigandChange);
}

function handleVirusChange() {
    const virusName = $(this).val();
    const pdbSelect = $('#pdb_code');
    pdbSelect.empty().append('<option value="">--Select PDB Code--</option>');
    if (virusName) {
        fetch(`/get_pdb_codes/${virusName}`)
            .then(response => response.json())
            .then(data => {
                data.pdb_codes.forEach(pdb => {
                    pdbSelect.append(new Option(pdb, pdb));
                });
                pdbSelect.prop('disabled', false);
            });
    } else {
        pdbSelect.prop('disabled', true);
        $('#ligand').prop('disabled', true);
        $('#chain').empty().append('<option value="">--Select Chain--</option>').prop('disabled', true);
    }
}

function handlePDBCodeChange() {
    const pdbCode = $(this).val();
    if (pdbCode) {
        loadLigands();
        fetch(`/check_functional_groups/${pdbCode}`)
            .then(response => response.json())
            .then(data => {
                $('#functional_groups').prop('disabled', !data.has_functional_groups);
                $('#functional_group_label').toggleClass('disabled', !data.has_functional_groups);
            });
    } else {
        $('#ligand').prop('disabled', true);
        $('#chain').empty().append('<option value="">--Select Chain--</option>').prop('disabled', true);
    }
}

function loadLigands() {
    const pdbCode = $('#pdb_code').val();
    const ligandSelect = $('#ligand');

    ligandSelect.empty().append('<option value="">--Select Ligand--</option>');

    if (pdbCode) {
        fetch(`/get_ligands/${pdbCode}`)
            .then(response => response.json())
            .then(data => {
                const occurrenceLigands = [];
                const uniqueLigands = [];
                const seenLigands = new Set();
                (data.ligands || []).forEach(ligand => {
                    const code = String((ligand && ligand.ligand) || '').trim().toUpperCase();
                    if (!code) return;
                    const normalized = Object.assign({}, ligand, { ligand: code });
                    occurrenceLigands.push(normalized);
                    if (!seenLigands.has(code)) {
                        seenLigands.add(code);
                        uniqueLigands.push(normalized);
                    }
                });

                uniqueLigands.forEach(ligand => {
                    ligandSelect.append(new Option(ligand.ligand, ligand.ligand));
                });

                // Re-initialize select2
                ligandSelect.select2();
                ligandSelect.prop('disabled', false).data('ligands', occurrenceLigands);
            });
    }
}





function handleLigandChange() {
    const selectedLigand = $(this).val();
    const ligands = $(this).data('ligands');
    const chainSelect = $('#chain');
    $('#ligand_instance_id').val('');
    chainSelect.empty().append('<option value="">--Select Chain--</option>');

    if (ligands) {
        ligands.forEach(item => {
            if (item.ligand === selectedLigand) {
                const instanceId = String(item.ligand_instance_id || '').trim();
                const label = instanceId
                    ? `Chain ${item.chain} · residue ${item.ligand_id || '?'} · model ${item.model_id || '1'}`
                    : item.chain;
                const option = new Option(label, item.chain);
                option.dataset.ligandInstanceId = instanceId;
                chainSelect.append(option);
            }
        });
        $('#chain-container').show();
        chainSelect.prop('disabled', false);
    }
}

$('#chain').on('change', function() {
    $('#chain_hidden').val($(this).val());  // Set the hidden input value to the selected chain
    const selectedChain = $(this).val();
    $('#ligand_instance_id').val(this.options[this.selectedIndex]?.dataset.ligandInstanceId || '');
    if (selectedChain) {
        $('#generateLigandImagesButton').prop('disabled', false);
    } else {
        $('#generateLigandImagesButton').prop('disabled', true);
    }
});

$('#chain').on('change', function() {
    $('#chain_hidden').val($(this).val());  // Set the hidden input value to the selected chain
    $('#generateLigandImagesButton').prop('disabled', false); // Enable the button once the chain is selected
});


$('#generateLigandImagesButton').on('click', function() {
    if ($('#ligand').val()) { // Ensure there is a selected ligand
        submitLigandImagesForm();
    } else {
        alert('Please select a ligand to generate images.');
    }
});


function submitLigandImagesForm() {
    $('#virus_hidden').val($('#virus').val());
    $('#pdb_code_hidden').val($('#pdb_code').val());
    $('#ligand_hidden').val($('#ligand').val());
    $('#ligand_instance_id_hidden').val($('#ligand_instance_id').val());
    $('#ligandImagesForm').submit();
}

function checkForPyMOLAnimation() {
    const urlParams = new URLSearchParams(window.location.search);
    const activatePymol = urlParams.get('activate_pymol');
    if (activatePymol === 'true') {
        animatePymolButton();
    }
}

function animatePymolButton() {
    const pymolButton = $('#pymol-integration-button');
    if (pymolButton.length) {
        pymolButton.css({
            transition: "transform 0.5s ease-in-out",
            transform: "scale(2)"
        });

        setTimeout(() => {
            pymolButton.css('transform', 'scale(1)');
        }, 2000);
    }
}

function redirectToComingSoon() {
    window.location.href = "/coming-soon"; // Adjust the path to match your Flask route
}


function animatePymolButton() {
    const pymolButton = document.getElementById('pymol-integration-button');
    if (pymolButton) {
        pymolButton.style.transition = "transform 0.5s ease-in-out";
        pymolButton.style.transform = "scale(2)";

        setTimeout(() => {
            pymolButton.style.transform = "scale(1)";
        }, 2000);
    }
}

window.onload = function () {
    const urlParams = new URLSearchParams(window.location.search);
    const activatePymol = urlParams.get('activate_pymol');

    if (activatePymol === 'true') {
        animatePymolButton();
    }
};









// ################################################################
// ################Scripts For Display Images Page#################
// ################################################################
// ################################################################
// ################################################################


// Update the hidden inputs when the form is submitted
function submitLigandImagesForm() {
    document.getElementById('virus_hidden').value = document.getElementById('virus').value;
    document.getElementById('pdb_code_hidden').value = document.getElementById('pdb_code').value;
    document.getElementById('ligand_hidden').value = document.getElementById('ligand').value;
    document.getElementById('ligand_instance_id_hidden').value = document.getElementById('ligand_instance_id').value;
    showLigandImageGenerationPromo();
    document.getElementById('ligandImagesForm').submit();
}

let ligandImageGenerationPromoTimer = null;
let ligandImageGenerationPromoIndex = 0;

function showLigandImageGenerationPromo() {
    const overlay = document.getElementById('ligand-image-generation-promo');
    if (!overlay) return;

    const cards = Array.from(overlay.querySelectorAll('.ligand-loading-promo-card'));
    const dots = overlay.querySelector('.ligand-loading-promo-dots');
    if (!cards.length) return;

    const renderPromo = (index) => {
        ligandImageGenerationPromoIndex = ((index % cards.length) + cards.length) % cards.length;
        cards.forEach((card, cardIndex) => card.classList.toggle('is-active', cardIndex === ligandImageGenerationPromoIndex));
        if (dots) {
            dots.innerHTML = cards.map((_, cardIndex) =>
                `<span class="${cardIndex === ligandImageGenerationPromoIndex ? 'is-active' : ''}"></span>`
            ).join('');
        }
    };

    const chooseNextPromo = () => {
        if (cards.length === 1) return 0;
        const availableIndexes = cards.map((_, index) => index).filter((index) => index !== ligandImageGenerationPromoIndex);
        return availableIndexes[Math.floor(Math.random() * availableIndexes.length)];
    };

    overlay.hidden = false;
    renderPromo(Math.floor(Math.random() * cards.length));
    if (ligandImageGenerationPromoTimer) window.clearInterval(ligandImageGenerationPromoTimer);
    ligandImageGenerationPromoTimer = window.setInterval(() => renderPromo(chooseNextPromo()), 4200);
}


function displayLigandInteractionDiagram(pdbId, ligandCode) {
    const container = document.getElementById('ligand-interaction-container');
    const loadingModal = document.getElementById("loadingModal");
    const spinnerContainer = document.getElementById("spinner-container");
    const errorMessage = document.getElementById("error-message");

    // Reset the modal and error message, then show the modal
    spinnerContainer.innerHTML = `<div class="spinner"></div>`;
    errorMessage.style.display = "none";
    loadingModal.style.display = "block";

    // Start polling the PoseView API after showing the spinner
    setTimeout(() => {
        fetch('https://proteins.plus/api/poseview_rest', {
            method: 'POST',
            headers: {
                'Accept': 'application/json',
                'Content-Type': 'application/json',
            },
            body: JSON.stringify({
                poseview: {
                    pdbCode: pdbId,
                    ligand: ligandCode
                }
            })
        })
        .then(response => {
            if (!response.ok) {
                return response.json().then(err => { throw new Error(err.message); });
            }
            return response.json();
        })
        .then(data => {
            if (data.status_code === 202 || data.status_code === 200) {
                // Poll the API for the PoseView job status
                pollPoseViewJob(data.location);
            } else {
                // If the API request fails, hide the spinner and show the error message
                showError();
            }
        })
        .catch(error => {
            showError();
        });
    }, 100);  // A small delay to allow the spinner to render
}

function showError() {
    const spinnerContainer = document.getElementById("spinner-container");
    const errorMessage = document.getElementById("error-message");

    // Replace the spinner with the error message
    spinnerContainer.innerHTML = '';
    errorMessage.style.display = "block";
}

function hideSpinner() {
    const loadingModal = document.getElementById("loadingModal");
    // Hide the modal once the job is done
    loadingModal.style.display = "none";
}




function toggleChatbot() {
    const chatbotPopup = document.getElementById("chatbot-popup");
    chatbotPopup.style.display =
        chatbotPopup.style.display === "none" || chatbotPopup.style.display === ""
            ? "block"
            : "none";
}
