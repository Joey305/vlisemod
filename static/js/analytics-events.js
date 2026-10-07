/* Privacy-safe, fixed V-LiSEMOD workflow taxonomy. Never reads form values. */
(function () {
  const allowed = new Set([
    'analysis_builder','protein_query','ligand_query','ligand_comparison','protacability',
    'ligand_interactions','protacability_search','pymol_session','ligand_images',
    'protein_query_export','protacability_evidence_export','builder_from_landing',
    'builder_from_ligand_query','builder_from_ligand_comparison','builder_from_protacability',
    'builder_from_viral_protac_design'
  ]);
  const newUuid = () => {
    if (crypto.randomUUID) return crypto.randomUUID();
    const bytes = crypto.getRandomValues(new Uint8Array(16)); bytes[6] = (bytes[6] & 15) | 64; bytes[8] = (bytes[8] & 63) | 128;
    return [...bytes].map((b, i) => `${b.toString(16).padStart(2, '0')}${[3,5,7,9].includes(i) ? '-' : ''}`).join('');
  };
  const send = (event_type, feature, failure_stage, handoffId) => {
    if (!allowed.has(feature)) return;
    const payload = { event_type, feature, failure_stage: failure_stage || '', handoff_id: handoffId || '', event_id: newUuid() };
    fetch('/analytics/event', { method: 'POST', headers: {'Content-Type': 'application/json'}, keepalive: true, body: JSON.stringify(payload) }).catch(() => {});
  };
  const originForPage = () => ({
    '/': 'builder_from_landing', '/ligand_query': 'builder_from_ligand_query',
    '/compare_ligands': 'builder_from_ligand_comparison', '/protacability_page': 'builder_from_protacability',
    '/viral-protac-design': 'builder_from_viral_protac_design'
  })[location.pathname] || 'builder_from_landing';
  const openBuilderHandoff = ({url, origin, target = '_blank', open = true}) => {
    if (!allowed.has(origin)) return null;
    const finalUrl = new URL(url, location.href);
    if (!/(^|\.)protacbuilder\.com$/i.test(finalUrl.hostname)) return null;
    const handoffId = newUuid();
    finalUrl.searchParams.set('utm_source', 'vlisemod');
    finalUrl.searchParams.set('utm_medium', 'ecosystem_referral');
    finalUrl.searchParams.set('utm_content', origin);
    finalUrl.searchParams.set('handoff_id', handoffId);
    send('companion_handoff', origin, '', handoffId);
    if (open) window.open(finalUrl.toString(), target, 'noopener,noreferrer');
    return finalUrl.toString();
  };
  const prepareBuilderHandoff = (link) => {
    if (link.dataset.analyticsHandoff) return;
    const url = new URL(link.href, location.href);
    if (!/(^|\.)protacbuilder\.com$/i.test(url.hostname)) return;
    const feature = originForPage();
    link.href = openBuilderHandoff({url: url.toString(), origin: feature, open: false}) || url.toString();
    link.dataset.analyticsHandoff = '1';
  };
  window.vlismodAnalytics = { send, prepareBuilderHandoff, openBuilderHandoff };
  document.addEventListener('click', (event) => {
    const target = event.target.closest('a,button');
    if (!target) return;
    if (target.matches('[data-analysis-builder-open], [data-analysis-context]')) send('workflow_started', 'analysis_builder');
    if (target.closest('#showInteractionButton')) { send('workflow_started', 'ligand_query'); send('analysis_submitted', 'ligand_interactions'); }
    if (target.closest('#showComparisonButton')) { send('workflow_started', 'ligand_comparison'); send('analysis_submitted', 'ligand_comparison'); }
    if (target.closest('#apply-protac-filters')) { send('workflow_started', 'protacability'); send('analysis_submitted', 'protacability_search'); }
    if (target.closest('#exportButton')) send('analysis_submitted', 'protein_query');
    if (target.closest('#export-protac-csv')) send('analysis_submitted', 'protacability_search');
    if (target.tagName === 'A' && (target.matches('[data-builder-handoff], .external-link, .nav-external'))) prepareBuilderHandoff(target);
  }, true);
  document.addEventListener('submit', (event) => {
    const form = event.target;
    if (form.matches('#analysis-builder-form')) send('analysis_submitted', 'pymol_session');
    if (form.matches('#analysis-ligand-images-form')) send('analysis_submitted', 'ligand_images');
  });
}());
