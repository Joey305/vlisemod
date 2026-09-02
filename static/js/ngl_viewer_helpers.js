(function () {
  const HELPER_VERSION = 'rcsb-mmcif-direct-v1';
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
      try { entry.stage.setSpin(false); } catch (e) {}
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
      proteinStructureExt: 'cif',
      ligandSdfUrl: '',
      ligandLabelAsymId: '',
      ligandMetadataUrl: '',
      attachmentSerialMap: {}
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

    const directMmcifUrl = `https://files.rcsb.org/download/${pdbCode}.cif`;
    const stage = new NGL.Stage(containerId, { backgroundColor: opts.backgroundColor || '#fbfdfc' });
    viewerRegistry.set(containerId, {
      stage: stage,
	      component: null,
	      proteinComponent: null,
	      ligandComponent: null,
	      options: opts,
	      attachmentSerialMap: normalizeSerialIndexMap(opts.attachmentSerialMap),
	      proteinSele: 'protein',
      ligandSele: '',
      defaultFocus: resolveDefaultFocusMode(opts, Boolean(opts.ligandSdfUrl || opts.ligandLabelAsymId))
    });

    const explicitStructureUrl = String(opts.structureUrl || '').trim();
    const proteinStructureUrl = String(opts.proteinStructureUrl || directMmcifUrl).trim();

    function buildLigandSdfUrl(labelAsymId) {
      const label = String(labelAsymId || '').trim().toUpperCase();
      const residue = encodeURIComponent(String(parsedContext.resno || opts.ligandResidueId || '').trim());
      if (!label || !residue || !ligandResname) return '';
      const pdbLower = pdbCode.toLowerCase();
      return `https://models.rcsb.org/v1/${pdbLower}/ligand?auth_seq_id=${residue}&label_asym_id=${encodeURIComponent(label)}&encoding=sdf&filename=${pdbLower}_${encodeURIComponent(label)}_${encodeURIComponent(ligandResname)}.sdf`;
    }

    function parseLabelAsymIdFromMmcif(text) {
      const ligand = ligandResname;
      const authChain = parsedContext.chain;
      const authResidue = normalizeResno(parsedContext.resno);
      const lines = String(text || '').split(/\r?\n/);
      let headers = [];
      for (let index = 0; index < lines.length; index += 1) {
        const line = lines[index].trim();
        if (line === 'loop_') { headers = []; continue; }
        if (line.startsWith('_atom_site.')) { headers.push(line); continue; }
        if (!headers.length || !line || line.startsWith('_')) continue;
        const parts = line.split(/\s+/);
        if (parts.length < headers.length) continue;
        const values = Object.fromEntries(headers.map((header, headerIndex) => [header, parts[headerIndex]]));
        if (String(values['_atom_site.group_PDB'] || '').toUpperCase() !== 'HETATM') continue;
        const comp = String(values['_atom_site.auth_comp_id'] || values['_atom_site.label_comp_id'] || '').toUpperCase();
        const chain = normalizeChain(values['_atom_site.auth_asym_id']);
        const residue = normalizeResno(values['_atom_site.auth_seq_id']);
        if (comp === ligand && chain === authChain && residue === authResidue) {
          return String(values['_atom_site.label_asym_id'] || '').trim().toUpperCase();
        }
      }
      return '';
    }

    async function resolveLigandSdfUrl() {
      if (String(opts.ligandSdfUrl || '').trim()) return String(opts.ligandSdfUrl).trim();
      let labelAsymId = String(opts.ligandLabelAsymId || '').trim().toUpperCase();
      if (!labelAsymId && String(opts.ligandMetadataUrl || '').trim()) {
        try {
          const response = await fetch(opts.ligandMetadataUrl, { cache: 'no-store' });
          if (response.ok) labelAsymId = String((await response.json()).label_asym_id || '').trim().toUpperCase();
        } catch (error) {
          ligandDebugLog('[VLNGLViewer] local ligand occurrence metadata unavailable', error);
        }
      }
      if (!labelAsymId) {
        try {
          const response = await fetch(directMmcifUrl, { cache: 'force-cache' });
          if (response.ok) labelAsymId = parseLabelAsymIdFromMmcif(await response.text());
        } catch (error) {
          ligandDebugLog('[VLNGLViewer] RCSB ligand metadata lookup failed', error);
        }
      }
      return buildLigandSdfUrl(labelAsymId);
    }

    function mmcifWithoutChemicalComponentLoops(text) {
      const lines = String(text || '').split(/\r?\n/);
      const output = [];
      for (let index = 0; index < lines.length;) {
        if (lines[index].trim() !== 'loop_') {
          output.push(lines[index]);
          index += 1;
          continue;
        }
        let headerEnd = index + 1;
        while (headerEnd < lines.length && lines[headerEnd].trim().startsWith('_')) headerEnd += 1;
        const headers = lines.slice(index + 1, headerEnd);
        const isChemicalComponentLoop = headers.some(function (header) {
          return header.trim().startsWith('_chem_comp.');
        });
        if (!isChemicalComponentLoop) {
          output.push(lines[index]);
          index += 1;
          continue;
        }
        index = headerEnd;
        while (index < lines.length && lines[index].trim() !== '#') index += 1;
        if (index < lines.length) index += 1;
      }
      return output.join('\n');
    }

    function loadStructure(url, extLabel) {
      ligandDebugLog('[NGL ligand debug] load attempt', { url: url, fileType: extLabel });
      const extension = extLabel || 'cif';
      const sourcePromise = extension === 'cif'
        ? fetch(url, { cache: 'force-cache' }).then(function (response) {
            if (!response.ok) throw new Error(`RCSB mmCIF request failed (${response.status})`);
            return response.text();
          }).then(function (text) {
            return new Blob([mmcifWithoutChemicalComponentLoops(text)], { type: 'chemical/x-cif' });
          })
        : Promise.resolve(url);
      return sourcePromise.then(function (source) {
        return stage.loadFile(source, { ext: extension, defaultRepresentation: false });
      }).then(function (component) {
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
      console.log('[VLNGLViewer] loadedExt', componentInfo.loadInfo ? componentInfo.loadInfo.fileType : 'cif');
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
	        attachmentSerialMap: normalizeSerialIndexMap(opts.attachmentSerialMap),
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
	        attachmentSerialMap: normalizeSerialIndexMap(opts.attachmentSerialMap),
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

    const proteinUrlToUse = opts.useDirectRcsbPdb === true ? directMmcifUrl : (proteinStructureUrl || explicitStructureUrl || directMmcifUrl);
    const proteinPromise = loadStructure(proteinUrlToUse, opts.proteinStructureExt || 'cif');
    const ligandSdfUrl = await resolveLigandSdfUrl();
    if (ligandSdfUrl) {
      try {
        const ligandPromise = loadLigandStructure(ligandSdfUrl).catch(function (ligandErr) {
          console.error('[VLNGLViewer] ligand SDF load failed', ligandErr);
          return null;
        });
        const [proteinInfo, ligandInfo] = await Promise.all([proteinPromise, ligandPromise]);
        return await renderSeparateProteinAndLigand(proteinInfo, ligandInfo);
      } catch (separateErr) {
        console.error('[VLNGLViewer] separate protein/ligand rendering failed; falling back to single-file mode', separateErr);
        if (typeof opts.onStatus === 'function') opts.onStatus('Protein or ligand component failed to load; falling back to single-file viewer.');
      }
    }

    let pdbComponentInfo;
    try {
      pdbComponentInfo = await proteinPromise;
      ligandDebugLog('[VLNGLViewer] mmCIF loaded', pdbCode);
    } catch (pdbLoadErr) {
      console.error('[VLNGLViewer] RCSB mmCIF load failed for', pdbCode, pdbLoadErr);
      if (typeof opts.onStatus === 'function') opts.onStatus('The 3D viewer could not load RCSB coordinates. Stored analysis results remain available.');
      throw pdbLoadErr;
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

  function focusLigand(containerId, duration) {
    const entry = viewerRegistry.get(containerId);
    if (!entry) return;
    if (entry.ligandComponent && entry.ligandAdded) {
      try { entry.ligandComponent.autoView(undefined, duration || 900); return; } catch (e) {}
    }
    if (!entry.component) return;
    if (!entry.ligandAdded || !entry.ligandSele) {
      entry.component.autoView(entry.proteinSele, duration || 900);
      return;
    }
    entry.component.autoView(entry.ligandSele, duration || 900);
  }

  function stopSpin(containerId) {
    const entry = viewerRegistry.get(containerId);
    if (!entry || !entry.stage) return false;
    try { entry.stage.setSpin(false); } catch (e) {}
    entry.isSpinning = false;
    const container = document.getElementById(containerId);
    if (container) container.dataset.ligandSpin = 'off';
    return true;
  }

  function setLigandSpin(containerId, enabled, options) {
    const entry = viewerRegistry.get(containerId);
    if (!entry || !entry.stage || typeof entry.stage.setSpin !== 'function') return false;
    const shouldSpin = enabled !== false;
    const container = document.getElementById(containerId);
    if (!shouldSpin) return stopSpin(containerId);
    if (!entry.ligandAdded) return false;
    const opts = options || {};
    try {
      focusLigand(containerId, opts.focusDuration || 450);
      entry.stage.setSpin(true);
      entry.isSpinning = true;
      if (container) container.dataset.ligandSpin = 'on';
      return true;
    } catch (e) {
      console.warn('[VLNGLViewer] ligand spin failed', e);
      stopSpin(containerId);
      return false;
    }
  }

  function toggleLigandSpin(containerId, enabled, options) {
    const entry = viewerRegistry.get(containerId);
    if (!entry) return false;
    const nextState = typeof enabled === 'boolean' ? enabled : !entry.isSpinning;
    return setLigandSpin(containerId, nextState, options);
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

  function cleanPositiveIntegers(values) {
    return (Array.isArray(values) ? values : [])
      .map(function (value) { return parseInt(value, 10); })
      .filter(function (value) { return Number.isFinite(value) && value > 0; });
  }

  function cleanNonNegativeIntegers(values) {
    return (Array.isArray(values) ? values : [])
      .map(function (value) { return parseInt(value, 10); })
      .filter(function (value) { return Number.isFinite(value) && value >= 0; });
  }

  function atomIndexSelection(atomIndices) {
    const indices = cleanNonNegativeIntegers(atomIndices);
    if (!indices.length) return '';
    return indices.map(function (index) { return `@${index}`; }).join(' or ');
  }

  function atomSerialSelection(atomSerials) {
    return atomIndexSelection(cleanPositiveIntegers(atomSerials));
  }

  function normalizeSerialIndexMap(value) {
    const map = {};
    const entries = value && typeof value === 'object' ? Object.entries(value) : [];
    entries.forEach(function (entry) {
      const serial = parseInt(entry[0], 10);
      const atomIndex = parseInt(entry[1], 10);
      if (Number.isFinite(serial) && serial > 0 && Number.isFinite(atomIndex) && atomIndex >= 0) {
        map[String(serial)] = atomIndex;
      }
    });
    return map;
  }

  function ligandIndicesFromSerials(entry, atomSerials) {
    const serialMap = normalizeSerialIndexMap(entry && entry.attachmentSerialMap);
    const seen = new Set();
    return cleanPositiveIntegers(atomSerials).map(function (serial) {
      const mapped = serialMap[String(serial)];
      return Number.isFinite(mapped) ? mapped : null;
    }).filter(function (atomIndex) {
      if (!Number.isFinite(atomIndex) || atomIndex < 0 || seen.has(atomIndex)) return false;
      seen.add(atomIndex);
      return true;
    });
  }

  function removeRepresentationFromAnyComponent(entry, repr) {
    if (!entry || !repr) return;
    if (Array.isArray(repr)) {
      repr.forEach(function (item) { removeRepresentationFromAnyComponent(entry, item); });
      return;
    }
    const components = [entry.ligandComponent, entry.component, entry.proteinComponent];
    const seen = new Set();
    components.forEach(function (component) {
      if (!component || seen.has(component)) return;
      seen.add(component);
      try { component.removeRepresentation(repr); } catch (e) {}
    });
  }

  function clearAttachmentHighlights(containerId) {
    const entry = viewerRegistry.get(containerId);
    if (!entry) return;
    removeRepresentationFromAnyComponent(entry, entry.attachmentRepr);
    removeRepresentationFromAnyComponent(entry, entry.attachmentSurfaceRepr);
    removeRepresentationFromAnyComponent(entry, entry.attachmentReprs);
    entry.attachmentRepr = null;
    entry.attachmentSurfaceRepr = null;
    entry.attachmentReprs = [];
    entry.lastAttachmentHighlight = null;
    const container = document.getElementById(containerId);
    if (container) {
      delete container.dataset.attachmentHighlight;
    }
  }

  function setAttachmentSerialMap(containerId, serialMap) {
    const entry = viewerRegistry.get(containerId);
    if (!entry) return false;
    entry.attachmentSerialMap = normalizeSerialIndexMap(serialMap);
    const container = document.getElementById(containerId);
    if (container) {
      container.dataset.attachmentSerialMapSize = String(Object.keys(entry.attachmentSerialMap).length);
    }
    return true;
  }

  function mapAttachmentAtomsToLigandIndices(containerId, attachmentAtoms) {
    const entry = viewerRegistry.get(containerId);
    if (!entry || !entry.ligandComponent || !entry.ligandComponent.structure) return {};
    const ligandAtoms = [];
    try {
      entry.ligandComponent.structure.eachAtom(function (atom) {
        ligandAtoms.push({ index: atom.index, x: atom.x, y: atom.y, z: atom.z });
      });
    } catch (error) {
      console.warn('[VLNGLViewer] unable to read loaded ligand coordinates', error);
      return {};
    }
    const serialMap = {};
    (Array.isArray(attachmentAtoms) ? attachmentAtoms : []).forEach(function (attachmentAtom) {
      const serial = Number(attachmentAtom && attachmentAtom.pdb_atom_serial);
      const x = Number(attachmentAtom && attachmentAtom.x);
      const y = Number(attachmentAtom && attachmentAtom.y);
      const z = Number(attachmentAtom && attachmentAtom.z);
      if (!Number.isFinite(serial) || !Number.isFinite(x) || !Number.isFinite(y) || !Number.isFinite(z)) return;
      let closest = null;
      let closestDistanceSquared = Infinity;
      ligandAtoms.forEach(function (ligandAtom) {
        const distanceSquared = (ligandAtom.x - x) ** 2 + (ligandAtom.y - y) ** 2 + (ligandAtom.z - z) ** 2;
        if (distanceSquared < closestDistanceSquared) {
          closest = ligandAtom;
          closestDistanceSquared = distanceSquared;
        }
      });
      if (closest && closestDistanceSquared <= 0.01) serialMap[String(serial)] = closest.index;
    });
    setAttachmentSerialMap(containerId, serialMap);
    return serialMap;
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

  function highlightAttachmentRegionSets(containerId, regionSets, options) {
    const entry = viewerRegistry.get(containerId);
    const fallbackComponent = entry && (entry.component || entry.proteinComponent);
    if (!entry || (!entry.ligandComponent && !fallbackComponent)) return false;
    const opts = options || {};
    const regions = (Array.isArray(regionSets) ? regionSets : [])
      .map(function (region, index) {
        return {
          regionId: String(region && region.regionId || `region-${index + 1}`),
          color: String(region && region.color || opts.candidateColor || opts.color || '#d94f3d'),
          candidateSerials: cleanPositiveIntegers((region && region.candidateSerials) || []),
          surfaceSerials: cleanPositiveIntegers((region && region.surfaceSerials) || [])
        };
      })
      .filter(function (region) { return region.candidateSerials.length || region.surfaceSerials.length; });
    clearAttachmentHighlights(containerId);
    if (!regions.length) return false;

    const hasLigandSerialMap = Object.keys(entry.attachmentSerialMap || {}).length > 0;
    const hasLigandSelection = entry.ligandComponent && entry.ligandAdded && regions.some(function (region) {
      const candidateSele = atomIndexSelection(ligandIndicesFromSerials(entry, region.candidateSerials));
      const surfaceSele = atomIndexSelection(ligandIndicesFromSerials(entry, region.surfaceSerials));
      return (candidateSele && countAtomsForSelection(entry.ligandComponent, candidateSele) > 0) ||
        (surfaceSele && countAtomsForSelection(entry.ligandComponent, surfaceSele) > 0);
    });
    const useLigandComponent = entry.ligandComponent && entry.ligandAdded && (hasLigandSelection || hasLigandSerialMap);

    const targetComponent = useLigandComponent ? entry.ligandComponent : fallbackComponent;
    const highlightRecords = [];
    let focusSele = '';
    try {
      entry.attachmentReprs = [];
      regions.forEach(function (region) {
        const candidateSerialSet = new Set(region.candidateSerials);
        const surfaceSerials = region.surfaceSerials.filter(function (serial) { return !candidateSerialSet.has(serial); });
        const ligandCandidateIndices = ligandIndicesFromSerials(entry, region.candidateSerials);
        const ligandSurfaceIndices = ligandIndicesFromSerials(entry, surfaceSerials)
          .filter(function (atomIndex) { return ligandCandidateIndices.indexOf(atomIndex) === -1; });
        const candidateSele = useLigandComponent ? atomIndexSelection(ligandCandidateIndices) : atomSerialSelection(region.candidateSerials);
        const surfaceSele = useLigandComponent ? atomIndexSelection(ligandSurfaceIndices) : atomSerialSelection(surfaceSerials);
        if (surfaceSele) {
          entry.attachmentReprs.push(targetComponent.addRepresentation('spacefill', {
            sele: surfaceSele,
            color: region.color,
            radiusScale: opts.surfaceRadiusScale || 0.38,
            opacity: opts.surfaceOpacity || 0.30,
            transparent: true,
            name: `attachment-site-surface-atoms-${region.regionId}`
          }));
        }
        if (candidateSele) {
          entry.attachmentReprs.push(targetComponent.addRepresentation(opts.representation || 'spacefill', {
            sele: candidateSele,
            color: region.color,
            radiusScale: opts.candidateRadiusScale || opts.radiusScale || 0.50,
            opacity: opts.candidateOpacity || opts.opacity || 0.48,
            transparent: true,
            name: `attachment-site-candidate-atoms-${region.regionId}`
          }));
        }
        if (!focusSele && (candidateSele || surfaceSele)) focusSele = candidateSele || surfaceSele;
        highlightRecords.push({
          regionId: region.regionId,
          color: region.color,
          candidateSerials: region.candidateSerials,
          surfaceSerials: surfaceSerials,
          candidateLigandIndices: ligandCandidateIndices,
          surfaceLigandIndices: ligandSurfaceIndices
        });
      });
      if (!entry.attachmentReprs.length || !focusSele) return false;
      targetComponent.autoView(focusSele, opts.duration || 800);
      entry.lastAttachmentHighlight = {
        mode: useLigandComponent ? 'ligand-sdf-index' : 'component-atom-index',
        regions: highlightRecords
      };
      const container = document.getElementById(containerId);
      if (container) {
        container.dataset.attachmentHighlight = JSON.stringify(entry.lastAttachmentHighlight);
      }
      console.log('[VLNGLViewer] attachment highlight', JSON.stringify(entry.lastAttachmentHighlight));
      return true;
    } catch (e) {
      console.warn('[VLNGLViewer] attachment-site highlight failed', e);
      clearAttachmentHighlights(containerId);
      return false;
    }
  }

  function highlightAttachmentSerialSets(containerId, serialSets, options) {
    const opts = options || {};
    return highlightAttachmentRegionSets(containerId, [{
      regionId: opts.regionId || 'attachment-site',
      color: opts.candidateColor || opts.color || '#d94f3d',
      candidateSerials: Array.isArray(serialSets && serialSets.candidateSerials) ? serialSets.candidateSerials : [],
      surfaceSerials: Array.isArray(serialSets && serialSets.surfaceSerials) ? serialSets.surfaceSerials : []
    }], opts);
  }

  function resizeViewer(containerId) {
    const entry = viewerRegistry.get(containerId);
    if (!entry || !entry.stage) return;
    try { entry.stage.handleResize(); } catch (e) {}
  }

  function debugViewerState(containerId) {
    const entry = viewerRegistry.get(containerId);
    if (!entry) return null;
    return {
      helperVersion: HELPER_VERSION,
      hasProteinComponent: Boolean(entry.proteinComponent || entry.component),
      hasLigandComponent: Boolean(entry.ligandComponent),
      ligandAdded: Boolean(entry.ligandAdded),
      attachmentSerialMapSize: Object.keys(entry.attachmentSerialMap || {}).length,
      lastAttachmentHighlight: entry.lastAttachmentHighlight || null
    };
  }

  window.VLNGLViewer = {
    initCleanProteinLigandViewer: initCleanProteinLigandViewer,
    initProteinWithLigandSdfViewer: initCleanProteinLigandViewer,
    disposeViewer: disposeViewer,
    resetView: resetView,
    focusLigand: focusLigand,
    fitAll: fitAll,
	    stopSpin: stopSpin,
	    setLigandSpin: setLigandSpin,
	    toggleLigandSpin: toggleLigandSpin,
	    toggleSurface: toggleSurface,
	    highlightAtomSerials: highlightAtomSerials,
	    highlightAttachmentSerialSets: highlightAttachmentSerialSets,
	    highlightAttachmentRegionSets: highlightAttachmentRegionSets,
    clearAttachmentHighlights: clearAttachmentHighlights,
    setAttachmentSerialMap: setAttachmentSerialMap,
    mapAttachmentAtomsToLigandIndices: mapAttachmentAtomsToLigandIndices,
    resizeViewer: resizeViewer,
	    ligandDebugEnabled: ligandDebugEnabled,
	    makeLigandCodeAliases: makeLigandCodeAliases,
	    debugViewerState: debugViewerState,
	    VERSION: HELPER_VERSION
  };
  console.log('[VLNGLViewer] helper version', HELPER_VERSION);
})();
