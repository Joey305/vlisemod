/* Privacy-safe, fixed V-LiSEMOD workflow taxonomy. Never reads form values. */
(function () {
  const allowed = new Set([
    'analysis_builder','protein_query','ligand_query','ligand_comparison','protacability',
    'ligand_interactions','protacability_search','pymol_session','ligand_images',
    'protein_query_export','protacability_evidence_export','builder_from_landing',
    'builder_from_ligand_query','builder_from_ligand_comparison','builder_from_protacability',
    'builder_from_viral_protac_design'
  ]);
  const send = (event_type, feature, failure_stage) => {
    if (!allowed.has(feature)) return;
    const payload = { event_type, feature, failure_stage: failure_stage || '', event_id: crypto.randomUUID ? crypto.randomUUID() : `${Date.now()}-${Math.random()}` };
    fetch('/analytics/event', { method: 'POST', headers: {'Content-Type': 'application/json'}, keepalive: true, body: JSON.stringify(payload) }).catch(() => {});
  };
  const originForPage = () => ({
    '/': 'builder_from_landing', '/ligand_query': 'builder_from_ligand_query',
    '/compare_ligands': 'builder_from_ligand_comparison', '/protacability_page': 'builder_from_protacability',
    '/viral-protac-design': 'builder_from_viral_protac_design'
  })[location.pathname] || 'builder_from_landing';
  const prepareBuilderHandoff = (link) => {
    if (link.dataset.analyticsHandoff) return;
    const url = new URL(link.href, location.href);
    if (!/protacbuilder\.com$/i.test(url.hostname)) return;
    const feature = originForPage();
    const handoffId = crypto.randomUUID ? crypto.randomUUID() : `${Date.now()}-${Math.random()}`;
    url.searchParams.set('utm_source', 'vlisemod');
    url.searchParams.set('utm_medium', 'ecosystem_referral');
    url.searchParams.set('utm_content', feature);
    url.searchParams.set('handoff_id', handoffId);
    link.href = url.toString(); link.dataset.analyticsHandoff = '1';
    send('companion_handoff', feature);
  };
  window.vlismodAnalytics = { send, prepareBuilderHandoff };
  document.addEventListener('click', (event) => {
    const target = event.target.closest('a,button');
    if (!target) return;
    if (target.matches('[data-analysis-builder-open], [data-analysis-context]')) send('workflow_started', 'analysis_builder');
    if (target.closest('#showInteractionButton')) { send('workflow_started', 'ligand_query'); send('analysis_submitted', 'ligand_interactions'); }
    if (target.closest('#showComparisonButton')) { send('workflow_started', 'ligand_comparison'); send('analysis_submitted', 'ligand_comparison'); }
    if (target.closest('#apply-protac-filters')) { send('workflow_started', 'protacability'); send('analysis_submitted', 'protacability_search'); }
    if (target.closest('#exportButton')) send('analysis_submitted', 'protein_query');
    if (target.closest('#export-protac-csv')) send('analysis_submitted', 'protacability_search');
    if (target.tagName === 'A') prepareBuilderHandoff(target);
  }, true);
  document.addEventListener('submit', (event) => {
    const form = event.target;
    if (form.matches('#analysis-builder-form')) send('analysis_submitted', 'pymol_session');
    if (form.matches('#analysis-ligand-images-form')) send('analysis_submitted', 'ligand_images');
  });
}());
