(function () {
    const modal = document.getElementById('analysis-builder-modal');
    if (!modal) return;

    const virus = document.getElementById('analysis-virus');
    const pdb = document.getElementById('analysis-pdb-code');
    const ligand = document.getElementById('analysis-ligand');
    const chain = document.getElementById('analysis-chain');
    const chainContainer = document.getElementById('analysis-chain-container');
    const instanceId = document.getElementById('analysis-ligand-instance-id');
    const imageButton = document.getElementById('analysis-generate-images');
    const imageForm = document.getElementById('analysis-ligand-images-form');
    const functionalGroups = document.getElementById('analysis-functional-groups');
    const functionalGroupLabel = document.getElementById('analysis-functional-group-label');
    let ligandOccurrences = [];
    let lastFocusedElement = null;

    if (typeof window.showLigandImageGenerationPromo !== 'function') {
        let promoTimer = null;
        let promoIndex = 0;
        window.showLigandImageGenerationPromo = () => {
            const overlay = document.getElementById('ligand-image-generation-promo');
            if (!overlay) return;
            const cards = Array.from(overlay.querySelectorAll('.ligand-loading-promo-card'));
            const dots = overlay.querySelector('.ligand-loading-promo-dots');
            if (!cards.length) return;
            const render = (index) => {
                promoIndex = ((index % cards.length) + cards.length) % cards.length;
                cards.forEach((card, cardIndex) => card.classList.toggle('is-active', cardIndex === promoIndex));
                if (dots) dots.innerHTML = cards.map((_, cardIndex) => `<span class="${cardIndex === promoIndex ? 'is-active' : ''}"></span>`).join('');
            };
            const nextIndex = () => cards.length < 2 ? 0 : cards.map((_, index) => index).filter((index) => index !== promoIndex)[Math.floor(Math.random() * (cards.length - 1))];
            overlay.hidden = false;
            render(Math.floor(Math.random() * cards.length));
            if (promoTimer) window.clearInterval(promoTimer);
            promoTimer = window.setInterval(() => render(nextIndex()), 4200);
        };
    }

    const resetSelect = (select, label, disabled) => {
        select.replaceChildren(new Option(label, ''));
        select.disabled = disabled;
    };

    const open = () => {
        lastFocusedElement = document.activeElement;
        modal.hidden = false;
        document.body.classList.add('analysis-builder-open');
        window.setTimeout(() => virus.focus({ preventScroll: true }), 50);
    };
    const close = () => {
        modal.hidden = true;
        document.body.classList.remove('analysis-builder-open');
        if (lastFocusedElement && typeof lastFocusedElement.focus === 'function') lastFocusedElement.focus({ preventScroll: true });
    };

    const requestJson = (url) => fetch(url).then((response) => {
        if (!response.ok) throw new Error(`Request failed (${response.status})`);
        return response.json();
    });

    const resolveVirusForPdb = async (pdbCode) => {
        if (!pdbCode) return '';
        const data = await requestJson('/get_viruses');
        for (const candidate of data.viruses || []) {
            const pdbs = await requestJson(`/get_pdb_codes/${encodeURIComponent(candidate)}`);
            if ((pdbs.pdb_codes || []).some((value) => String(value).toUpperCase() === String(pdbCode).toUpperCase())) return candidate;
        }
        return '';
    };

    window.openAnalysisBuilderWithContext = async (context = {}) => {
        open();
        const normalized = {
            virus: String(context.virus || '').trim(),
            pdbCode: String(context.pdbCode || '').trim().toUpperCase(),
            ligand: String(context.ligand || '').trim().toUpperCase(),
            chain: String(context.chain || '').trim(),
            ligandInstanceId: String(context.ligandInstanceId || '').trim()
        };
        try {
            const virusName = normalized.virus || await resolveVirusForPdb(normalized.pdbCode);
            if (!virusName) return;
            virus.value = virusName;
            resetSelect(pdb, '--Select PDB Code--', true);
            const pdbData = await requestJson(`/get_pdb_codes/${encodeURIComponent(virusName)}`);
            (pdbData.pdb_codes || []).forEach((value) => pdb.add(new Option(value, value)));
            pdb.disabled = false;
            if (!normalized.pdbCode || ![...pdb.options].some((option) => option.value.toUpperCase() === normalized.pdbCode)) return;
            pdb.value = [...pdb.options].find((option) => option.value.toUpperCase() === normalized.pdbCode).value;

            const ligandData = await requestJson(`/get_ligands/${encodeURIComponent(pdb.value)}`);
            ligandOccurrences = ligandData.ligands || [];
            resetSelect(ligand, '--Select Ligand--', true);
            [...new Set(ligandOccurrences.map((item) => String(item.ligand || '').trim().toUpperCase()).filter(Boolean))]
                .forEach((value) => ligand.add(new Option(value, value)));
            ligand.disabled = false;
            if (!normalized.ligand || ![...ligand.options].some((option) => option.value === normalized.ligand)) return;
            ligand.value = normalized.ligand;

            resetSelect(chain, '--Select Chain--', true);
            ligandOccurrences.filter((item) => String(item.ligand || '').trim().toUpperCase() === normalized.ligand).forEach((item) => {
                const occurrence = String(item.ligand_instance_id || '').trim();
                const label = occurrence ? `Chain ${item.chain} · residue ${item.ligand_id || '?'} · model ${item.model_id || '1'}` : `Chain ${item.chain}`;
                const option = new Option(label, item.chain);
                option.dataset.ligandInstanceId = occurrence;
                chain.add(option);
            });
            chainContainer.hidden = false;
            chain.disabled = false;
            const requestedOption = [...chain.options].find((option) => option.dataset.ligandInstanceId === normalized.ligandInstanceId)
                || [...chain.options].find((option) => option.value === normalized.chain);
            if (requestedOption) {
                chain.value = requestedOption.value;
                instanceId.value = requestedOption.dataset.ligandInstanceId || normalized.ligandInstanceId;
                imageButton.disabled = false;
            }
        } catch (error) {
            console.warn('[analysis-builder] unable to prefill context', error);
        }
    };

    document.querySelectorAll('[data-analysis-builder-open]').forEach((trigger) => {
        trigger.addEventListener('click', (event) => {
            // Explorer links remain usable without JavaScript, but open the shared
            // in-page builder whenever the popup runtime is available.
            event.preventDefault();
            open();
        });
    });
    document.addEventListener('click', (event) => {
        const trigger = event.target.closest('[data-analysis-context]');
        if (!trigger) return;
        event.preventDefault();
        window.openAnalysisBuilderWithContext({
            virus: trigger.dataset.analysisVirus,
            pdbCode: trigger.dataset.analysisPdbCode,
            ligand: trigger.dataset.analysisLigand,
            chain: trigger.dataset.analysisChain,
            ligandInstanceId: trigger.dataset.analysisLigandInstanceId
        });
    });
    modal.querySelector('.analysis-builder-close').addEventListener('click', close);
    modal.addEventListener('click', (event) => { if (event.target === modal) close(); });
    document.addEventListener('keydown', (event) => { if (event.key === 'Escape' && !modal.hidden) close(); });

    fetch('/get_viruses').then((response) => response.json()).then((data) => {
        (data.viruses || []).forEach((value) => virus.add(new Option(value, value)));
    }).catch(() => {});

    virus.addEventListener('change', () => {
        resetSelect(pdb, '--Select PDB Code--', true);
        resetSelect(ligand, '--Select Ligand--', true);
        resetSelect(chain, '--Select Chain--', true);
        chainContainer.hidden = true;
        imageButton.disabled = true;
        instanceId.value = '';
        if (!virus.value) return;
        fetch(`/get_pdb_codes/${encodeURIComponent(virus.value)}`).then((response) => response.json()).then((data) => {
            (data.pdb_codes || []).forEach((value) => pdb.add(new Option(value, value)));
            pdb.disabled = false;
        }).catch(() => {});
    });

    pdb.addEventListener('change', () => {
        resetSelect(ligand, '--Select Ligand--', true);
        resetSelect(chain, '--Select Chain--', true);
        chainContainer.hidden = true;
        imageButton.disabled = true;
        instanceId.value = '';
        if (!pdb.value) return;
        fetch(`/get_ligands/${encodeURIComponent(pdb.value)}`).then((response) => response.json()).then((data) => {
            ligandOccurrences = data.ligands || [];
            [...new Set(ligandOccurrences.map((item) => String(item.ligand || '').trim().toUpperCase()).filter(Boolean))]
                .forEach((value) => ligand.add(new Option(value, value)));
            ligand.disabled = false;
        }).catch(() => {});
        fetch(`/check_functional_groups/${encodeURIComponent(pdb.value)}`).then((response) => response.json()).then((data) => {
            functionalGroups.disabled = !data.has_functional_groups;
            functionalGroupLabel.classList.toggle('disabled', !data.has_functional_groups);
        }).catch(() => {});
    });

    ligand.addEventListener('change', () => {
        resetSelect(chain, '--Select Chain--', true);
        instanceId.value = '';
        imageButton.disabled = true;
        const selected = ligand.value;
        ligandOccurrences.filter((item) => String(item.ligand || '').trim().toUpperCase() === selected).forEach((item) => {
            const occurrence = String(item.ligand_instance_id || '').trim();
            const label = occurrence ? `Chain ${item.chain} · residue ${item.ligand_id || '?'} · model ${item.model_id || '1'}` : `Chain ${item.chain}`;
            const option = new Option(label, item.chain);
            option.dataset.ligandInstanceId = occurrence;
            chain.add(option);
        });
        chainContainer.hidden = false;
        chain.disabled = false;
    });

    chain.addEventListener('change', () => {
        instanceId.value = chain.options[chain.selectedIndex]?.dataset.ligandInstanceId || '';
        imageButton.disabled = !chain.value;
    });

    imageButton.addEventListener('click', () => {
        if (!virus.value || !pdb.value || !ligand.value || !chain.value) return;
        const fields = imageForm.elements;
        fields.virus.value = virus.value;
        fields.pdb_code.value = pdb.value;
        fields.ligand.value = ligand.value;
        fields.chain.value = chain.value;
        fields.ligand_instance_id.value = instanceId.value;
        if (typeof window.showLigandImageGenerationPromo === 'function') window.showLigandImageGenerationPromo();
        imageForm.submit();
    });
})();
