/* Sends only fixed workflow labels; no form values or scientific content. */
(function () {
  const send = (event_type, feature, failure_stage) => {
    const payload = { event_type, feature: feature || '', failure_stage: failure_stage || '', event_id: crypto.randomUUID ? crypto.randomUUID() : `${Date.now()}-${Math.random()}` };
    fetch('/analytics/event', { method: 'POST', headers: {'Content-Type': 'application/json'}, keepalive: true, body: JSON.stringify(payload) }).catch(() => {});
  };
  window.vlismodAnalytics = { send };
  document.addEventListener('click', (event) => {
    if (event.target.closest('[data-analysis-builder-open]')) send('workflow_started', 'analysis_builder');
  });
  document.addEventListener('submit', (event) => {
    const form = event.target;
    if (form.matches('#analysis-builder-form')) send('analysis_submitted', 'pymol_session');
    if (form.matches('#analysis-ligand-images-form')) send('analysis_submitted', 'ligand_images');
    if (form.matches('[data-analytics-export]')) send('export_generated', 'data_export');
  });
}());
