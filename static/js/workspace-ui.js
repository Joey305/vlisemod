/* Shared Phase-1 research-workspace UI.  Scientific identity stays in data
 * attributes/URL parameters; this helper never reconstructs it from labels. */
(function (global) {
    'use strict';
    const methods = {
        mapping: ['Atom mapping', 'CIF-to-SMILES mapping', 'legacy_mcs_etkdg_uff_cif_v2.5'],
        sasa: ['Solvent exposure', 'Shrake–Rupley SASA, 1.40 Å probe', 'biopython-shrake_rupley-1.40-cif-v2.1'],
        interactions: ['Interactions', 'PDBe Arpeggio, occurrence-resolved contact analysis', 'arpeggio-cif-v2.2'],
        geometry: ['Geometry', 'CIF ligand geometry', 'cif-ligand-geometry-v2.4'],
        functionalGroups: ['Functional groups', 'RDKit SMARTS functional groups', 'rdkit-smarts-functional-groups-v2.3'],
        protacability: ['PROTACability', 'Ligand-centered and target-context triage', 'protacability-cif-v2.8'],
        attachment: ['Attachment sites', 'Atom-specific attachment-site evidence', 'attachment-sites-cif-v2.6']
    };
    const esc = (value) => String(value == null ? '' : value).replace(/[&<>'"]/g, c => ({'&':'&amp;','<':'&lt;','>':'&gt;',"'":'&#39;','"':'&quot;'}[c]));
    const paramsUrl = (path, params) => {
        const search = new URLSearchParams();
        Object.entries(params || {}).forEach(([key, value]) => {
            if (Array.isArray(value)) value.forEach(item => item != null && item !== '' && search.append(key, item));
            else if (value != null && value !== '') search.set(key, value);
        });
        return path + (search.toString() ? '?' + search.toString() : '');
    };
    function occurrenceCard(ctx) {
        const details = [ctx.pdbCode, ctx.chain && `Chain ${ctx.chain}`, ctx.residueId && `Residue ${ctx.residueId}`].filter(Boolean);
        const facts = [];
        if (ctx.mappedAtoms != null) facts.push(`${ctx.mappedAtoms} mapped atom${Number(ctx.mappedAtoms) === 1 ? '' : 's'}`);
        if (ctx.exposedAtoms != null) facts.push(`${ctx.exposedAtoms} solvent-exposed atom${Number(ctx.exposedAtoms) === 1 ? '' : 's'}`);
        const data = {
            'data-ligand-instance-id': ctx.ligandInstanceId,
            'data-ligand-code': ctx.ligand,
            'data-pdb-code': ctx.pdbCode,
            'data-chain': ctx.chain,
            'data-residue-id': ctx.residueId,
            'data-model-id': ctx.modelId
        };
        const attrs = Object.entries(data).filter(([, v]) => v != null && v !== '').map(([k, v]) => `${k}="${esc(v)}"`).join(' ');
        return `<article class="workspace-occurrence-card" ${attrs}><strong>${esc(ctx.ligand || 'Ligand occurrence')}</strong><span>${esc(details.join(' · '))}</span>${facts.length ? `<small>${esc(facts.join(' · '))}</small>` : ''}</article>`;
    }
    function ensure(id) { return document.getElementById(id); }
    function render(id, config) {
        const root = ensure(id); if (!root) return;
        if (!config) { root.innerHTML = ''; root.classList.add('hidden'); return; }
        const context = config.context || {};
        const chips = (context.chips || []).filter(Boolean).map(item => `<span class="workspace-chip">${esc(item)}</span>`).join('');
        const canonicalUrl = config.url || window.location.pathname + window.location.search;
        const summary = config.summary ? `<section class="workspace-summary"><span class="workspace-eyebrow">Evidence summary</span><h2>${esc(config.summary.title || '')}</h2><div class="workspace-metrics">${(config.summary.metrics || []).filter(m => m && m.label).map(m => `<div><strong>${esc(m.value)}</strong><span>${esc(m.label)}</span></div>`).join('')}</div>${config.summary.note ? `<p>${esc(config.summary.note)}</p>` : ''}</section>` : '';
        const actions = (config.actions || []).filter(a => a && a.href).map(a => `<a class="${esc(a.className || 'button-secondary')}" ${a.newTab ? 'target="_blank" rel="noopener noreferrer"' : ''} href="${esc(a.href)}">${esc(a.label)}</a>`).join('');
        const provenance = (config.methods || []).map(key => methods[key]).filter(Boolean).map(item => `<li><strong>${esc(item[0])}</strong><span>${esc(item[1])}</span><small>${esc(item[2])}</small></li>`).join('');
        root.classList.remove('hidden');
        root.innerHTML = `<div class="workspace-context"><div><span class="workspace-eyebrow">Current context</span><div class="workspace-chips">${chips}</div></div><div class="workspace-utility"><button type="button" class="button-secondary workspace-copy-link" data-url="${esc(canonicalUrl)}">Copy link</button>${config.clearHref ? `<a class="button-secondary" href="${esc(config.clearHref)}">Clear</a>` : ''}</div></div>${context.occurrence ? occurrenceCard(context.occurrence) : ''}${summary}${actions ? `<section class="workspace-next"><span class="workspace-eyebrow">What next?</span><div class="action-row">${actions}</div></section>` : ''}${provenance ? `<details class="workspace-provenance"><summary>Analysis details</summary><ul>${provenance}</ul></details>` : ''}`;
    }
    function noData(id, title, message, action) {
        const root = ensure(id); if (!root) return;
        root.classList.remove('hidden');
        root.innerHTML = `<section class="workspace-no-data" role="status"><h2>${esc(title)}</h2><p>${esc(message)}</p>${action && action.href ? `<a class="button-secondary" href="${esc(action.href)}">${esc(action.label)}</a>` : ''}</section>`;
    }
    async function copyText(text, button) {
        try { await navigator.clipboard.writeText(new URL(text, window.location.origin).href); }
        catch (_) { const field = document.createElement('textarea'); field.value = new URL(text, window.location.origin).href; document.body.appendChild(field); field.select(); document.execCommand('copy'); field.remove(); }
        if (button) { const original = button.textContent; button.textContent = 'Link copied'; button.setAttribute('aria-label', 'Link copied'); setTimeout(() => { button.textContent = original; button.removeAttribute('aria-label'); }, 1800); }
    }
    document.addEventListener('click', event => { const button = event.target.closest('.workspace-copy-link'); if (button) copyText(button.dataset.url || window.location.href, button); });
    global.VLiSEMODWorkspace = { render, noData, paramsUrl, occurrenceCard, methods, copyText };
})(window);
