(function () {
  const HELPER_VERSION = 'pdb-ligand-default-focus-v4';
  const viewerRegistry = new Map();
  let ligandElementSchemeId = null;
  const WATER_RESNAMES = new Set(['HOH', 'WAT', 'DOD', 'H2O', 'TIP', 'SOL']);

  function normalizePdbCode(value) {
    return String(value || '').trim().toUpperCase();
  }

  function normalizeResno(value) {
    return String(value === null || value === undefined ? '' : value).trim().toUpperCase();
  }

  function normalizeChain(value) {
    return String(value || '').trim().toUpperCase();
  }

  function ligandDebugEnabled() {
    try {
      if (window.VLISEMOD_DEBUG_LIGANDS === true) return true;
      const params = new URLSearchParams(window.location.search);
      if (params.get('debug_ligands') === '1') return true;
      if (window.localStorage && localStorage.getItem('vlisemod_debug_ligands') === '1') return true;
    } catch (e) {}
    return false;
  }

  function ligandDebugLog() {
    if (!ligandDebugEnabled()) return;
    console.log.apply(console, arguments);
  }

  function makeLigandCodeAliases(code) {
    const raw = String(code || '').trim().toUpperCase();
    if (!raw) return [];
    const aliases = [raw];
    const oToZero = raw.replace(/O/g, '0');
    const zeroToO = raw.replace(/0/g, 'O');
    if (oToZero && aliases.indexOf(oToZero) === -1) aliases.push(oToZero);
    if (zeroToO && aliases.indexOf(zeroToO) === -1) aliases.push(zeroToO);
    return aliases;
  }

  function inferElementSymbol(atom) {
    const rawElement = String((atom && atom.element) || '').trim().toUpperCase();
    if (rawElement) {
      if (rawElement === 'CL' || rawElement === 'BR') return rawElement;
      if (rawElement.length > 1) return rawElement.replace(/[^A-Z]/g, '').charAt(0);
      return rawElement;
    }
    const atomName = String((atom && (atom.atomname || atom.atomName || atom.name)) || '').trim().toUpperCase();
    if (!atomName) return '';
    const letters = atomName.replace(/[^A-Z]/g, '');
    if (!letters) return '';
    if (letters.startsWith('CL')) return 'CL';
    if (letters.startsWith('BR')) return 'BR';
    return letters.charAt(0);
  }

  function isWaterResidue(rp, resname) {
    if (typeof rp.isWater === 'function' && rp.isWater()) return true;
    return WATER_RESNAMES.has(String(resname || '').trim().toUpperCase());
  }

  function parseLigandContext(opts) {
    const rawChain = String(opts.ligandChain || opts.chainId || '').trim();
    const rawResid = String(opts.ligandResidueId || opts.ligandResno || '').trim();
    let chain = normalizeChain(rawChain);
    let resno = normalizeResno(rawResid);

    const compactResid = rawResid.replace(/\s+/g, '').toUpperCase();
    const splitDirect = compactResid.match(/^([A-Z])([0-9]+[A-Z]?)$/);
    if (!chain && splitDirect) {
      chain = splitDirect[1];
      resno = splitDirect[2];
    }

    const spacedResid = rawResid.match(/^([A-Za-z])\s+([0-9]+[A-Za-z]?)$/);
    if (!chain && spacedResid) {
      chain = normalizeChain(spacedResid[1]);
      resno = normalizeResno(spacedResid[2]);
    }

    if (chain && resno && resno.startsWith(chain) && /^([A-Z][0-9]+[A-Z]?)$/.test(resno)) {
      resno = resno.slice(1);
    }

    return {
      chain: chain,
      resno: resno,
      raw: {
        ligandChain: rawChain,
        ligandResidueId: String(opts.ligandResidueId || '').trim(),
        ligandResno: String(opts.ligandResno || '').trim(),
        chainId: String(opts.chainId || '').trim()
      }
    };
  }

  function disposeViewer(containerId) {
    const entry = viewerRegistry.get(containerId);
    if (entry && entry.stage) {
      try {
        if (typeof entry.stage.removeAllComponents === 'function') entry.stage.removeAllComponents();
      } catch (e) {}
      try { entry.stage.dispose(); } catch (e) {}
    }
    viewerRegistry.delete(containerId);
  }

  function getLigandElementScheme() {
    if (ligandElementSchemeId) return ligandElementSchemeId;
    if (typeof NGL === 'undefined' || !NGL.ColormakerRegistry) return null;
    try {
      ligandElementSchemeId = NGL.ColormakerRegistry.addScheme(function () {
        this.atomColor = function (atom) {
          const element = inferElementSymbol(atom);
          if (element === 'C') return 0xf4a261;
          if (element === 'O') return 0xff3b30;
          if (element === 'N') return 0x1d4ed8;
          if (element === 'H') return 0xf3f4f6;
          if (element === 'S') return 0xfacc15;
          if (element === 'P') return 0xff69b4;
          if (element === 'F') return 0x9be564;
          if (element === 'CL') return 0x16a34a;
          if (element === 'BR') return 0x15803d;
          if (element === 'I') return 0x166534;
          return 0xf4a261;
        };
      });
    } catch (e) {
      console.warn('[VLNGLViewer] unable to register ligand element color scheme; falling back to uniform color', e);
      ligandElementSchemeId = null;
    }
    return ligandElementSchemeId;
  }

  function getLigandColorValue(opts) {
    const mode = String(opts.ligandColorMode || 'element').trim().toLowerCase();
    if (mode === 'uniform' || mode === 'warhead') {
      return opts.ligandColor || '#f4a261';
    }
    return getLigandElementScheme() || (opts.ligandColor || '#f4a261');
  }

  function residueToSelector(residue) {
    const rn = String(residue.resname || '').trim().toUpperCase();
    const ch = String(residue.chain || '').trim();
    const rs = normalizeResno(residue.resno);
    if (!rn || !rs) return '';
    if (ch) return `(resname ${rn}) and (:${ch}) and (${rs}) and not water`;
    return `(resname ${rn}) and (${rs}) and not water`;
  }

  function collectResidueClasses(component) {
    const proteinAAs = new Set(['ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','ILE','LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL','SEC','PYL','MSE']);
    const heteroResidues = [];
    let waterResidueCount = 0;

    component.structure.eachResidue(function (rp) {
      const resname = String(rp.resname || '').trim().toUpperCase();
      const chain = String(rp.chainname || rp.chainid || rp.chain || '').trim();
      const resno = normalizeResno(rp.resno || rp.residueNumber || '');
      if (!resname || !resno) return;

      const water = isWaterResidue(rp, resname);
      if (water) {
        waterResidueCount += 1;
        return;
      }

      const isProtein = (typeof rp.isProtein === 'function' && rp.isProtein()) || proteinAAs.has(resname);
      const isNucleic = (typeof rp.isNucleic === 'function' && rp.isNucleic());
      if (isProtein || isNucleic) return;

      let atomCount = 0;
      try { rp.eachAtom(function () { atomCount += 1; }); } catch (e) {}
      if (atomCount <= 0) return;

      heteroResidues.push({
        resname: resname,
        chain: chain,
        resno: resno,
        atomCount: atomCount,
        selector: residueToSelector({ resname: resname, chain: chain, resno: resno })
      });
    });

    return { heteroResidues: heteroResidues, waterResidueCount: waterResidueCount };
  }

  function countAtomsForSelection(component, sele) {
    if (!sele) return 0;
    let count = 0;
    try {
      const selection = new NGL.Selection(sele);
      component.structure.eachAtom(function () { count += 1; }, selection);
    } catch (e) {
      console.warn('[NGL ligand debug] selection count failed', sele, e);
    }
    return count;
  }

  function describeRepresentation(repr) {
    let selectionString = '';
    try {
      if (repr && repr.selection && repr.selection.string) selectionString = repr.selection.string;
      else if (repr && typeof repr.getSelection === 'function' && repr.getSelection()) selectionString = repr.getSelection().string || '';
      else if (repr && repr.parameters && repr.parameters.sele) selectionString = repr.parameters.sele;
    } catch (e) {}
    return {
      wrapperType: repr && repr.type,
      actualType: (repr && repr.repr && (repr.repr.type || repr.repr.name || (repr.repr.constructor && repr.repr.constructor.name))) || '',
      sele: selectionString || '',
      params: (repr && repr.parameters) || null
    };
  }

  function validateAndPruneRepresentations(component, allowedRepresentations) {
    const allowed = allowedRepresentations || new Set();
    (component.reprList || []).slice().forEach(function (repr) {
      if (!allowed.has(repr)) {
        console.warn('[VLNGLViewer] removing unexpected representation', describeRepresentation(repr));
        try { component.removeRepresentation(repr); } catch (e) {}
      }
    });
  }

  function resolveDefaultFocusMode(opts, hasLigandComponent) {
    const explicit = String((opts && opts.defaultFocus) || '').trim().toLowerCase();
    if (explicit === 'ligand' || explicit === 'protein' || explicit === 'all') return explicit;
    return hasLigandComponent ? 'ligand' : 'protein';
  }

  function applyDefaultFocus(entry, opts, duration) {
    if (!entry || !entry.stage) return false;
    const focusMode = resolveDefaultFocusMode(opts || entry.options || {}, Boolean(entry.ligandComponent && entry.ligandAdded));
    console.log('[VLNGLViewer] default focus', focusMode);

    if (focusMode === 'all') {
      try {
        entry.stage.autoView(duration || 900);
        return true;
      } catch (e) {
        console.warn('[VLNGLViewer] full-scene focus failed', e);
      }
    }

    if (focusMode === 'protein' && entry.proteinComponent) {
      try {
        entry.proteinComponent.autoView('protein', duration || 900);
        return true;
      } catch (e) {
        console.warn('[VLNGLViewer] protein focus failed; falling back', e);
      }
    }

    if (entry.ligandComponent && entry.ligandAdded) {
      try {
        entry.ligandComponent.autoView(undefined, duration || 900);
        console.log('[VLNGLViewer] centered ligand SDF component', entry.ligandAdded);
        return true;
      } catch (e) {
        console.warn('[VLNGLViewer] ligand focus failed; falling back to full scene', e);
      }
    }

    if (entry.proteinComponent) {
      try {
        entry.proteinComponent.autoView('protein', duration || 900);
        return true;
      } catch (e) {
        console.warn('[VLNGLViewer] fallback protein focus failed', e);
      }
    }

    if (entry.component && entry.component !== entry.proteinComponent) {
      try {
        entry.component.autoView(entry.proteinSele || 'protein', duration || 900);
        return true;
      } catch (e) {
        console.warn('[VLNGLViewer] legacy component focus failed', e);
      }
    }

    try {
      entry.stage.autoView(duration || 900);
      return true;
    } catch (e) {
      console.warn('[VLNGLViewer] final stage focus fallback failed', e);
    }
    return false;
  }

  async function initCleanProteinLigandViewer(containerId, options) {
    const opts = Object.assign({
      pdbCode: '',
      chainId: '',
      ligandResname: '',
      ligandChain: '',
      ligandResidueId: '',
      ligandResno: '',
      structureUrl: '',
      allowBroadLigandFallback: false,
      focusLigand: true,
      retryCifOnLigandMiss: false,
      renderAllHeteroWhenNoExact: true,
      renderAllHeteroOnExactMiss: false,
      proteinColor: '#005030',
      ligandColor: '#f4a261',
      ligandColorMode: 'element',
      defaultFocus: '',
      backgroundColor: '#fbfdfc',
      structureExt: 'pdb',
      useDirectRcsbPdb: false,
      proteinStructureUrl: '',
      proteinStructureExt: 'pdb',
      ligandSdfUrl: '',
      ligandLabelAsymId: ''
    }, options || {});

    const container = document.getElementById(containerId);
    if (!container) return Promise.reject(new Error('Viewer container not found'));
    if (typeof NGL === 'undefined') return Promise.reject(new Error('NGL library is not loaded'));

    const pdbCode = normalizePdbCode(opts.pdbCode);
    const ligandResname = String(opts.ligandResname || '').trim().toUpperCase();
    const ligandAliases = makeLigandCodeAliases(ligandResname);
    const parsedContext = parseLigandContext(opts);

    ligandDebugLog('[VLNGLViewer] incoming ligand options', {
      ligandResname: opts.ligandResname || '',
      ligandChain: opts.ligandChain || '',
      chainId: opts.chainId || '',
      ligandResidueId: opts.ligandResidueId || '',
      ligandResno: opts.ligandResno || ''
    });
    ligandDebugLog('[VLNGLViewer] parsed ligand context', parsedContext);

    if (!pdbCode || !/^[A-Z0-9]{4}$/.test(pdbCode)) {
      if (typeof opts.onStatus === 'function') opts.onStatus('No PDB code available for 3D context.');
      return Promise.reject(new Error('No PDB code available for 3D context.'));
    }

    disposeViewer(containerId);
    container.innerHTML = '';
    if (getComputedStyle(container).display === 'none') container.style.display = 'block';
    if (!container.style.height || container.clientHeight < 80) container.style.height = '420px';

    const ligandSdfUrl = String(opts.ligandSdfUrl || '').trim();
    const stage = new NGL.Stage(containerId, { backgroundColor: opts.backgroundColor || '#fbfdfc' });
    viewerRegistry.set(containerId, {
      stage: stage,
      component: null,
      proteinComponent: null,
      ligandComponent: null,
      options: opts,
      proteinSele: 'protein',
      ligandSele: '',
      defaultFocus: resolveDefaultFocusMode(opts, Boolean(ligandSdfUrl))
    });

    const pdbUrl = `https://files.rcsb.org/view/${pdbCode}.pdb`;
    const pdbDownloadUrl = `https://files.rcsb.org/download/${pdbCode}.pdb`;
    const explicitStructureUrl = String(opts.structureUrl || '').trim();
    const proteinStructureUrl = String(opts.proteinStructureUrl || '').trim();

    function loadStructure(url, extLabel) {
      ligandDebugLog('[NGL ligand debug] load attempt', { url: url, fileType: extLabel });
      return stage.loadFile(url, { ext: 'pdb', defaultRepresentation: false }).then(function (component) {
        ligandDebugLog('[NGL ligand debug] load success', { url: url, fileType: extLabel });
        return { component: component, loadInfo: { url: url, fileType: extLabel, success: true } };
      });
    }

    function loadLigandStructure(url) {
      ligandDebugLog('[NGL ligand debug] ligand load attempt', { url: url, fileType: 'sdf' });
      return stage.loadFile(url, { ext: 'sdf', defaultRepresentation: false }).then(function (component) {
        ligandDebugLog('[NGL ligand debug] ligand load success', { url: url, fileType: 'sdf' });
        return { component: component, loadInfo: { url: url, fileType: 'sdf', success: true } };
      });
    }

    function buildRenderPlan(component) {
      const classes = collectResidueClasses(component);
      const heteroResidues = classes.heteroResidues;
      const waterResidueCount = classes.waterResidueCount;
      const proteinAtomCount = countAtomsForSelection(component, 'protein');

      console.table(heteroResidues.map(function (r) {
        return {
          resname: r.resname,
          chain: r.chain,
          resno: r.resno,
          atomCount: r.atomCount,
          selector: r.selector
        };
      }));

      const exactRequested = Boolean(ligandResname && parsedContext.chain && parsedContext.resno);
      const exactMatches = heteroResidues.filter(function (r) {
        return ligandAliases.indexOf(String(r.resname || '').toUpperCase()) !== -1
          && normalizeChain(r.chain) === parsedContext.chain
          && normalizeResno(r.resno) === parsedContext.resno;
      });

      const resnameMatches = heteroResidues.filter(function (r) {
        return ligandAliases.indexOf(String(r.resname || '').toUpperCase()) !== -1;
      });

      const uniqueResnameContexts = Array.from(new Set(resnameMatches.map(function (r) {
        return `${r.resname}|${normalizeChain(r.chain)}|${normalizeResno(r.resno)}`;
      })));

      let mode = 'protein-only';
      let targetResidues = [];
      let status = 'Structure loaded. Protein cartoon rendered.';

      if (exactRequested) {
        if (exactMatches.length === 1) {
          mode = 'focused-ligand';
          targetResidues = exactMatches;
          status = `Structure loaded. Ligand ${ligandResname} highlighted as ${exactMatches[0].resname} ${exactMatches[0].chain}${exactMatches[0].resno}.`;
        } else if (exactMatches.length > 1) {
          status = `Structure loaded. Multiple exact ligand matches were found for ${ligandResname} ${parsedContext.resno}${parsedContext.chain}; no atom-level ligand rendering applied.`;
        } else if (opts.renderAllHeteroOnExactMiss && heteroResidues.length > 0) {
          mode = 'all-hetero';
          targetResidues = heteroResidues;
          status = `Structure loaded. Exact ligand not found; rendering all non-water hetero ligands (${heteroResidues.length}).`;
        } else {
          status = `Structure loaded. No matching ligand residue was found for ${ligandResname} ${parsedContext.resno}${parsedContext.chain}.`;
        }
      } else {
        if (ligandResname) {
          if (uniqueResnameContexts.length === 1 && resnameMatches.length > 0) {
            mode = 'focused-ligand';
            targetResidues = [resnameMatches[0]];
            status = `Structure loaded. Ligand ${ligandResname} highlighted by unique residue context.`;
          } else if (uniqueResnameContexts.length > 1) {
            status = `Structure loaded. Multiple ligand residues matched ${ligandResname}; select a specific chain/residue context.`;
          } else if (opts.renderAllHeteroWhenNoExact && heteroResidues.length > 0) {
            mode = 'all-hetero';
            targetResidues = heteroResidues;
            status = `Structure loaded. Rendering all non-water hetero ligands (${heteroResidues.length}).`;
          }
        } else if (opts.renderAllHeteroWhenNoExact && heteroResidues.length > 0) {
          mode = 'all-hetero';
          targetResidues = heteroResidues;
          status = `Structure loaded. Rendering all non-water hetero ligands (${heteroResidues.length}).`;
        }
      }

      const exactSelectors = targetResidues.map(function (r) { return r.selector; }).filter(Boolean);
      const ligandSelector = exactSelectors.length ? `(${exactSelectors.join(') or (')})` : '';

      console.log('[VLNGLViewer] protein atom count', proteinAtomCount);
      console.log('[VLNGLViewer] water residue count', waterResidueCount);
      console.log('[VLNGLViewer] hetero non-water residue list', heteroResidues);
      console.log('[VLNGLViewer] final ligand/hetero selector', ligandSelector || null);

      return {
        mode: mode,
        status: status,
        ligandSelector: ligandSelector,
        heteroResidues: heteroResidues
      };
    }

    async function renderFromComponent(componentInfo) {
      const component = componentInfo.component;
      console.log('[VLNGLViewer] helper version', HELPER_VERSION);
      console.log('[VLNGLViewer] loadedSourceUrl', componentInfo.loadInfo ? componentInfo.loadInfo.url : null);
      console.log('[VLNGLViewer] loadedExt', 'pdb');
      if (typeof opts.onSourceLoaded === 'function') {
        try { opts.onSourceLoaded(componentInfo.loadInfo ? componentInfo.loadInfo.url : ''); } catch (e) {}
      }
      component.removeAllRepresentations();

      const allowedRepresentations = new Set();
      const proteinRepr = component.addRepresentation('cartoon', {
        sele: 'protein',
        color: opts.proteinColor || '#005030',
        smoothSheet: true
      });
      allowedRepresentations.add(proteinRepr);

      const plan = buildRenderPlan(component);
      let ligandAdded = false;
      let ligandSele = '';
      let ligandRepr = null;

      if (plan.ligandSelector) {
        const ligandColorScheme = getLigandElementScheme();
        ligandRepr = component.addRepresentation('ball+stick', {
          sele: plan.ligandSelector,
          color: ligandColorScheme || (opts.ligandColor || '#f4a261'),
          multipleBond: true,
          radiusScale: 1.25,
          bondScale: 0.55,
          aspectRatio: 1.8,
          opacity: 1.0
        });
        allowedRepresentations.add(ligandRepr);
        ligandAdded = true;
        ligandSele = plan.ligandSelector;
      }

      validateAndPruneRepresentations(component, allowedRepresentations);
      const reprDebug = (component.reprList || []).map(describeRepresentation);
      console.log('[VLNGLViewer] final component.reprList', reprDebug);

      const focusSele = (opts.focusLigand && ligandAdded) ? ligandSele : 'protein';
      try { component.autoView(focusSele, 900); } catch (e) { component.autoView('protein', 900); }
      try { stage.handleResize(); } catch (e) {}

      if (typeof opts.onStatus === 'function') opts.onStatus(plan.status);

      viewerRegistry.set(containerId, {
        stage: stage,
        component: component,
        options: opts,
        proteinSele: 'protein',
        ligandSele: ligandSele,
        ligandAdded: ligandAdded,
        defaultFocus: resolveDefaultFocusMode(opts, false)
      });

      return {
        entry: viewerRegistry.get(containerId),
        ligandAdded: ligandAdded,
        heteroResidues: plan.heteroResidues,
        loadInfo: componentInfo.loadInfo
      };
    }

    async function renderSeparateProteinAndLigand(proteinComponentInfo, ligandComponentInfo) {
      const proteinComponent = proteinComponentInfo.component;
      const ligandComponent = ligandComponentInfo ? ligandComponentInfo.component : null;
      console.log('[VLNGLViewer] helper version', HELPER_VERSION);
      console.log('[VLNGLViewer] protein source URL', proteinComponentInfo.loadInfo ? proteinComponentInfo.loadInfo.url : null);
      console.log('[VLNGLViewer] ligand SDF URL', ligandComponentInfo && ligandComponentInfo.loadInfo ? ligandComponentInfo.loadInfo.url : null);
      if (typeof opts.onSourceLoaded === 'function') {
        try {
          opts.onSourceLoaded({
            proteinSourceUrl: proteinComponentInfo.loadInfo ? proteinComponentInfo.loadInfo.url : '',
            ligandSdfUrl: ligandComponentInfo && ligandComponentInfo.loadInfo ? ligandComponentInfo.loadInfo.url : ''
          });
        } catch (e) {}
      }

      proteinComponent.removeAllRepresentations();
      const proteinRepr = proteinComponent.addRepresentation('cartoon', {
        sele: 'protein',
        color: opts.proteinColor || '#005030',
        smoothSheet: true
      });

      const proteinAtomCount = countAtomsForSelection(proteinComponent, 'protein');
      console.log('[VLNGLViewer] protein atom count', proteinAtomCount);
      console.log('[VLNGLViewer] final protein representation cartoon', describeRepresentation(proteinRepr));

      let ligandRepr = null;
      let ligandAtomCount = 0;
      let ligandAdded = false;
      const ligandContextLabel = [ligandResname || '', parsedContext.chain || '', parsedContext.resno || ''].filter(Boolean).join(' ').trim();
      if (ligandComponent) {
        ligandComponent.removeAllRepresentations();
        ligandRepr = ligandComponent.addRepresentation('ball+stick', {
          sele: 'all',
          color: getLigandColorValue(opts),
          multipleBond: true,
          radiusScale: 1.25,
          bondScale: 0.55,
          aspectRatio: 1.8,
          opacity: 1.0
        });
        ligandAtomCount = countAtomsForSelection(ligandComponent, 'all');
        ligandAdded = true;
        console.log('[VLNGLViewer] ligand SDF atom count', ligandAtomCount);
        console.log('[VLNGLViewer] final ligand representation ball+stick all', describeRepresentation(ligandRepr));
      }

      const resolvedDefaultFocus = resolveDefaultFocusMode(opts, ligandAdded);
      const entry = {
        stage: stage,
        component: proteinComponent,
        proteinComponent: proteinComponent,
        ligandComponent: ligandComponent,
        options: Object.assign({}, opts, { defaultFocus: resolvedDefaultFocus }),
        proteinSele: 'protein',
        ligandSele: 'all',
        ligandAdded: ligandAdded,
        proteinRepr: proteinRepr,
        ligandRepr: ligandRepr,
        ligandSdfUrl: ligandSdfUrl,
        defaultFocus: resolvedDefaultFocus
      };

      console.log('[VLNGLViewer] components loaded: protein=%s ligand=%s', true, ligandAdded);
      let defaultFocusApplied = false;
      try { stage.handleResize(); } catch (e) {}
      defaultFocusApplied = applyDefaultFocus(entry, entry.options, 900);
      setTimeout(function () {
        try { stage.handleResize(); } catch (e) {}
        applyDefaultFocus(entry, entry.options, 700);
      }, 150);

      const status = ligandAdded && defaultFocusApplied
        ? `Structure loaded. Protein shown as cartoon; ${ligandContextLabel || (ligandResname || 'Ligand')} loaded from ligand instance SDF and centered.`
        : ligandAdded
        ? `Structure loaded. Protein shown as cartoon; ${ligandContextLabel || (ligandResname || 'Ligand')} loaded from ligand instance SDF.`
        : 'Protein loaded, but ligand instance SDF could not be loaded.';
      if (typeof opts.onStatus === 'function') opts.onStatus(status);

      viewerRegistry.set(containerId, entry);
      return viewerRegistry.get(containerId);
    }

    if (ligandSdfUrl) {
      try {
        const proteinUrlToUse = opts.useDirectRcsbPdb === true ? pdbUrl : (proteinStructureUrl || explicitStructureUrl || pdbUrl);
        const proteinInfo = await loadStructure(proteinUrlToUse, opts.proteinStructureExt || 'pdb');
        let ligandInfo = null;
        try {
          ligandInfo = await loadLigandStructure(ligandSdfUrl);
        } catch (ligandErr) {
          console.error('[VLNGLViewer] ligand SDF load failed', ligandErr);
        }
        return await renderSeparateProteinAndLigand(proteinInfo, ligandInfo);
      } catch (separateErr) {
        console.error('[VLNGLViewer] separate protein/ligand rendering failed; falling back to single-file mode', separateErr);
        if (typeof opts.onStatus === 'function') opts.onStatus('Protein or ligand component failed to load; falling back to single-file viewer.');
      }
    }

    let pdbComponentInfo;
    try {
      if (opts.useDirectRcsbPdb === true) {
        pdbComponentInfo = await loadStructure(pdbUrl, 'pdb');
      } else if (explicitStructureUrl) {
        pdbComponentInfo = await loadStructure(explicitStructureUrl, 'pdb');
      } else {
        pdbComponentInfo = await loadStructure(pdbUrl, 'pdb');
      }
      ligandDebugLog('[VLNGLViewer] PDB loaded', pdbCode);
    } catch (pdbLoadErr) {
      console.warn('[VLNGLViewer] primary PDB load failed, trying fallback PDB URL:', pdbDownloadUrl, pdbLoadErr);
      try {
        const fallbackPdbInfo = await loadStructure(pdbDownloadUrl, 'pdb');
        const fallbackRenderResult = await renderFromComponent(fallbackPdbInfo);
        return fallbackRenderResult.entry;
      } catch (loadFallbackErr) {
        console.error('[VLNGLViewer] all structure loading attempts failed for', pdbCode, loadFallbackErr);
        if (typeof opts.onStatus === 'function') opts.onStatus('Unable to load 3D context from RCSB right now.');
        throw loadFallbackErr;
      }
    }

    let pdbRenderResult;
    try {
      pdbRenderResult = await renderFromComponent(pdbComponentInfo);
    } catch (renderErr) {
      console.error('[VLNGLViewer] PDB loaded but render failed', renderErr);
      if (typeof opts.onStatus === 'function') opts.onStatus('Structure loaded, but viewer rendering failed. See console.');
      throw renderErr;
    }

    return pdbRenderResult.entry;
  }

  function resetView(containerId) {
    const entry = viewerRegistry.get(containerId);
    if (!entry || !entry.stage) return;
    if (applyDefaultFocus(entry, entry.options || {}, 900)) return;
    if (!entry.component) return;
    entry.component.autoView(entry.proteinSele, 900);
  }

  function focusLigand(containerId) {
    const entry = viewerRegistry.get(containerId);
    if (!entry) return;
    if (entry.ligandComponent && entry.ligandAdded) {
      try { entry.ligandComponent.autoView(undefined, 900); return; } catch (e) {}
    }
    if (!entry.component) return;
    if (!entry.ligandAdded || !entry.ligandSele) {
      entry.component.autoView(entry.proteinSele, 900);
      return;
    }
    entry.component.autoView(entry.ligandSele, 900);
  }

  function fitAll(containerId) {
    const entry = viewerRegistry.get(containerId);
    if (!entry || !entry.stage) return;
    try { entry.stage.autoView(900); } catch (e) {}
  }

  function toggleSurface(containerId, enabled) {
    const entry = viewerRegistry.get(containerId);
    const targetComponent = entry && (entry.proteinComponent || entry.component);
    if (!entry || !targetComponent) return;
    if (enabled) {
      entry.surfaceRepr = targetComponent.addRepresentation('surface', {
        sele: 'protein',
        opacity: 0.2,
        color: entry.options.proteinColor || '#005030'
      });
    } else if (entry.surfaceRepr) {
      targetComponent.removeRepresentation(entry.surfaceRepr);
      entry.surfaceRepr = null;
    }
  }

  function atomSerialSelection(atomSerials) {
    const serials = (Array.isArray(atomSerials) ? atomSerials : [])
      .map(function (value) { return parseInt(value, 10); })
      .filter(function (value) { return Number.isFinite(value) && value > 0; });
    if (!serials.length) return '';
    return serials.map(function (serial) { return `@${serial}`; }).join(' or ');
  }

  function clearAttachmentHighlights(containerId) {
    const entry = viewerRegistry.get(containerId);
    if (!entry) return;
    const targetComponent = entry.component || entry.proteinComponent;
    if (targetComponent && entry.attachmentRepr) {
      try { targetComponent.removeRepresentation(entry.attachmentRepr); } catch (e) {}
    }
    entry.attachmentRepr = null;
  }

  function highlightAtomSerials(containerId, atomSerials, options) {
    const entry = viewerRegistry.get(containerId);
    const targetComponent = entry && (entry.component || entry.proteinComponent);
    if (!entry || !targetComponent) return false;
    const opts = options || {};
    const sele = atomSerialSelection(atomSerials);
    clearAttachmentHighlights(containerId);
    if (!sele) return false;
    try {
      entry.attachmentRepr = targetComponent.addRepresentation(opts.representation || 'spacefill', {
        sele: sele,
        color: opts.color || '#d94f3d',
        radiusScale: opts.radiusScale || 0.72,
        opacity: opts.opacity || 1.0,
        name: 'attachment-site-highlight'
      });
      targetComponent.autoView(sele, opts.duration || 800);
      return true;
    } catch (e) {
      console.warn('[VLNGLViewer] attachment-site highlight failed', e);
      clearAttachmentHighlights(containerId);
      return false;
    }
  }

  function resizeViewer(containerId) {
    const entry = viewerRegistry.get(containerId);
    if (!entry || !entry.stage) return;
    try { entry.stage.handleResize(); } catch (e) {}
  }

  window.VLNGLViewer = {
    initCleanProteinLigandViewer: initCleanProteinLigandViewer,
    initProteinWithLigandSdfViewer: initCleanProteinLigandViewer,
    disposeViewer: disposeViewer,
    resetView: resetView,
    focusLigand: focusLigand,
    fitAll: fitAll,
    toggleSurface: toggleSurface,
    highlightAtomSerials: highlightAtomSerials,
    clearAttachmentHighlights: clearAttachmentHighlights,
    resizeViewer: resizeViewer,
    ligandDebugEnabled: ligandDebugEnabled,
    makeLigandCodeAliases: makeLigandCodeAliases,
    VERSION: HELPER_VERSION
  };
  console.log('[VLNGLViewer] helper version', HELPER_VERSION);
})();
