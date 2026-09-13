/* Offline, DOM-free taxonomy-flow model plus the small SVG viewer adapter. */
(function (root) {
  "use strict";

  const RANKS = ["D", "K", "P", "C", "O", "F", "G", "S"];
  const RANK_INDEX = Object.fromEntries(RANKS.map((rank, index) => [rank, index]));
  const MAX_SAFE_INTEGER = 9007199254740991;
  const ENTRY_ID = "synthetic:entry:classified";
  const NODE_KINDS = {entry: 0, taxon: 1, residual: 2, residual_carry: 3};
  const LINK_KINDS = {biological: 0, other_hidden: 1, assigned_above: 2, carry: 3};
  const SUBTYPE_ORDER = {other_hidden: 0, assigned_above: 1};
  const COUNT_MODEL = "classified clade-read flow with persistent explicit residual lanes";

  function fail(message) { throw new Error(message); }

  function safeInteger(value, label) {
    if (!Number.isSafeInteger(value) || value < 0 || value > MAX_SAFE_INTEGER) {
      fail(label + " must be a non-negative JavaScript-safe integer");
    }
    return value;
  }

  function compareBytes(left, right) {
    const a = new TextEncoder().encode(left);
    const b = new TextEncoder().encode(right);
    const length = Math.min(a.length, b.length);
    for (let i = 0; i < length; i += 1) {
      if (a[i] !== b[i]) return a[i] - b[i];
    }
    return a.length - b.length;
  }

  function rotr(value, bits) {
    return (value >>> bits) | (value << (32 - bits));
  }

  /* Small synchronous SHA-256 implementation; required for deterministic IDs
     in a browser without a network or asynchronous crypto dependency. */
  const SHA_K = [
    0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5,
    0x3956c25b, 0x59f111f1, 0x923f82a4, 0xab1c5ed5,
    0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3,
    0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174,
    0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc,
    0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da,
    0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7,
    0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967,
    0x27b70a85, 0x2e1b2138, 0x4d2c6dfc, 0x53380d13,
    0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85,
    0xa2bfe8a1, 0xa81a664b, 0xc24b8b70, 0xc76c51a3,
    0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070,
    0x19a4c116, 0x1e376c08, 0x2748774c, 0x34b0bcb5,
    0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3,
    0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208,
    0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2,
  ];

  function sha256(text) {
    const source = new TextEncoder().encode(text);
    const bitLength = source.length * 8;
    const paddedLength = (((source.length + 9 + 63) >> 6) << 6);
    const bytes = new Uint8Array(paddedLength);
    bytes.set(source);
    bytes[source.length] = 0x80;
    const view = new DataView(bytes.buffer);
    view.setUint32(paddedLength - 4, bitLength >>> 0);
    view.setUint32(paddedLength - 8, Math.floor(bitLength / 4294967296));
    let h0 = 0x6a09e667;
    let h1 = 0xbb67ae85;
    let h2 = 0x3c6ef372;
    let h3 = 0xa54ff53a;
    let h4 = 0x510e527f;
    let h5 = 0x9b05688c;
    let h6 = 0x1f83d9ab;
    let h7 = 0x5be0cd19;
    for (let offset = 0; offset < bytes.length; offset += 64) {
      const w = new Uint32Array(64);
      for (let i = 0; i < 16; i += 1) w[i] = view.getUint32(offset + i * 4);
      for (let i = 16; i < 64; i += 1) {
        const s0 = rotr(w[i - 15], 7) ^ rotr(w[i - 15], 18) ^ (w[i - 15] >>> 3);
        const s1 = rotr(w[i - 2], 17) ^ rotr(w[i - 2], 19) ^ (w[i - 2] >>> 10);
        w[i] = (w[i - 16] + s0 + w[i - 7] + s1) >>> 0;
      }
      let a = h0, b = h1, c = h2, d = h3;
      let e = h4, f = h5, g = h6, h = h7;
      for (let i = 0; i < 64; i += 1) {
        const s1 = rotr(e, 6) ^ rotr(e, 11) ^ rotr(e, 25);
        const ch = (e & f) ^ (~e & g);
        const temp1 = (h + s1 + ch + SHA_K[i] + w[i]) >>> 0;
        const s0 = rotr(a, 2) ^ rotr(a, 13) ^ rotr(a, 22);
        const maj = (a & b) ^ (a & c) ^ (b & c);
        const temp2 = (s0 + maj) >>> 0;
        h = g; g = f; f = e; e = (d + temp1) >>> 0;
        d = c; c = b; b = a; a = (temp1 + temp2) >>> 0;
      }
      h0 = (h0 + a) >>> 0; h1 = (h1 + b) >>> 0;
      h2 = (h2 + c) >>> 0; h3 = (h3 + d) >>> 0;
      h4 = (h4 + e) >>> 0; h5 = (h5 + f) >>> 0;
      h6 = (h6 + g) >>> 0; h7 = (h7 + h) >>> 0;
    }
    return [h0, h1, h2, h3, h4, h5, h6, h7]
      .map(value => value.toString(16).padStart(8, "0")).join("");
  }

  function nodeId(rank, path) { return "taxon:" + rank + ":" + sha256(path); }
  function laneId(subtype, origin, rank) { return sha256([subtype, origin, rank].join("\0")); }
  function residualId(lane, rank) { return "synthetic:residual:" + lane + ":" + rank; }
  function linkId(source, target, kind) { return "link:" + sha256([source, target, kind].join("\0")); }

  function parseRanks(ranks) {
    const values = Array.isArray(ranks) ? ranks.slice() : String(ranks).split(",");
    if (values.length < 2 || values.length > RANKS.length) fail("invalid rank selection");
    if (new Set(values).size !== values.length) fail("rank selection contains duplicates");
    if (values.some(rank => !RANK_INDEX.hasOwnProperty(rank))) fail("invalid rank selection");
    for (let i = 1; i < values.length; i += 1) {
      if (RANK_INDEX[values[i]] <= RANK_INDEX[values[i - 1]]) fail("rank selection is not canonical");
    }
    return values;
  }

  function sourceFromPayload(payload) {
    if (!payload || typeof payload !== "object") fail("Sankey payload must be an object");
    const totals = payload.totals || {};
    const total = safeInteger(totals.total, "totals.total");
    const classified = safeInteger(totals.classified, "totals.classified");
    const unclassified = safeInteger(totals.unclassified, "totals.unclassified");
    if (total !== classified + unclassified) fail("payload totals do not conserve");
    const nodes = Array.isArray(payload.source_nodes) ? payload.source_nodes : [];
    const byPath = new Map();
    const byId = new Map();
    const byRank = Object.fromEntries(RANKS.map(rank => [rank, []]));
    for (const raw of nodes) {
      if (!raw || typeof raw !== "object") fail("source node must be an object");
      const rank = raw.rank;
      if (!RANK_INDEX.hasOwnProperty(rank)) fail("source node has invalid rank");
      safeInteger(raw.clade, "source node clade");
      safeInteger(raw.direct, "source node direct");
      safeInteger(raw.source_order_index, "source node source_order_index");
      if (raw.direct > raw.clade || typeof raw.path !== "string") fail("invalid source node counts/path");
      const node = Object.assign({}, raw, {id: raw.id || nodeId(rank, raw.path)});
      if (byPath.has(node.path) || byId.has(node.id)) fail("duplicate source node");
      byPath.set(node.path, node); byId.set(node.id, node); byRank[rank].push(node);
    }
    for (const node of nodes) {
      let parent = node.parent_id ? byId.get(node.parent_id) : null;
      const parentPath = parent ? parent.path : "";
      if (parent && !node.path.startsWith(parent.path + ";")) fail("source path/parent mismatch");
      if (!parent && RANK_INDEX[node.rank] !== 0) fail("missing source ancestor");
      node.parent_path = parentPath;
      if (parent) {
        if (RANK_INDEX[parent.rank] + 1 !== RANK_INDEX[node.rank]) fail("non-adjacent source ancestor");
      }
    }
    const children = new Map();
    for (const node of nodes) {
      const key = node.parent_path || "";
      if (!children.has(key)) children.set(key, []);
      children.get(key).push(node);
    }
    for (const node of nodes) {
      const childSum = (children.get(node.path) || []).reduce((sum, child) => sum + child.clade, 0);
      if (node.clade !== node.direct + childSum) fail("source clade arithmetic failed");
    }
    if (nodes.length && byRank.D.reduce((sum, node) => sum + node.clade, 0) !== classified) {
      fail("source Domain clades do not equal classified reads");
    }
    if (nodes.reduce((sum, node) => sum + node.direct, 0) !== classified) {
      fail("source direct counts do not equal classified reads");
    }
    return {total, classified, unclassified, nodes, byPath, byId, byRank};
  }

  function selection(source, ranks, maxN) {
    if (!Number.isInteger(maxN) || maxN < 1 || maxN > 100) fail("invalid max N");
    const retained = Object.fromEntries(ranks.map(rank => [rank, new Set()]));
    const reasons = new Map();
    for (const rank of ranks) {
      const candidates = source.byRank[rank].filter(node => node.clade > 0).slice();
      candidates.sort((a, b) => b.clade - a.clade || compareBytes(a.path, b.path) ||
        a.source_order_index - b.source_order_index || compareBytes(a.id, b.id));
      for (const node of candidates.slice(0, maxN)) {
        retained[rank].add(node.path); reasons.set(node.path, "top_n");
      }
    }
    for (let later = 1; later < ranks.length; later += 1) {
      for (const path of Array.from(retained[ranks[later]])) {
        let current = source.byPath.get(path);
        for (const earlierRank of ranks.slice(0, later)) {
          while (current.rank !== earlierRank) {
            if (!current.parent_id) fail("retained node has missing ancestor");
            current = source.byId.get(current.parent_id);
          }
          retained[earlierRank].add(current.path);
          if (!reasons.has(current.path)) reasons.set(current.path, "ancestor_closure");
          current = source.byPath.get(path);
        }
      }
    }
    return {retained, reasons};
  }

  function nodeRecord(fields) {
    return {
      id: fields.id, kind: fields.kind, subtype: fields.subtype,
      biological: fields.biological, carried: fields.carried, lane_id: fields.lane_id,
      taxid: fields.taxid, rank: fields.rank, name: fields.name, path: fields.path,
      parent_id: fields.parent_id, clade: fields.clade, direct: fields.direct,
      status: fields.status, selection_reason: fields.selection_reason,
      source_order_index: fields.source_order_index, visual_order: null,
      origin_source_id: fields.origin_source_id, origin_target_rank: fields.origin_target_rank,
      column_rank: fields.column_rank, value: fields.value,
      member_count: fields.member_count, member_paths_sha256: fields.member_paths_sha256,
    };
  }

  function bioNode(sourceNode, rank, reason) {
    return nodeRecord({
      id: sourceNode.id, kind: "taxon", subtype: null, biological: true, carried: false,
      lane_id: null, taxid: sourceNode.taxid, rank: sourceNode.rank,
      name: sourceNode.name, path: sourceNode.path, parent_id: sourceNode.parent_id,
      clade: sourceNode.clade, direct: sourceNode.direct, status: sourceNode.status,
      selection_reason: reason, source_order_index: sourceNode.source_order_index,
      origin_source_id: null, origin_target_rank: null, column_rank: rank,
      value: sourceNode.clade, member_count: null, member_paths_sha256: null,
    });
  }

  function residualNode(subtype, origin, originIndex, originRank, columnRank, value, carried,
                       memberCount, memberDigest) {
    const lane = laneId(subtype, origin, originRank);
    return nodeRecord({
      id: residualId(lane, columnRank), kind: carried ? "residual_carry" : "residual",
      subtype, biological: false, carried, lane_id: lane, taxid: null, rank: null,
      name: subtype === "other_hidden" ? "Other " + columnRank : "Assigned above " + columnRank,
      path: null, parent_id: null, clade: null, direct: null, status: null,
      selection_reason: carried ? "carry" : "residual", source_order_index: originIndex,
      origin_source_id: origin, origin_target_rank: originRank, column_rank: columnRank,
      value, member_count: memberCount, member_paths_sha256: memberDigest,
    });
  }

  function link(source, target, kind, value, transition) {
    return {id: linkId(source, target, kind), source, target, kind, value,
      transition_index: transition};
  }

  function digestMembers(paths) {
    return sha256(paths.slice().sort(compareBytes).join("\n"));
  }

  function descendant(node, ancestorPath) {
    return node.path === ancestorPath || node.path.startsWith(ancestorPath + ";");
  }

  function assignVisual(nodes, source, ranks) {
    for (const rank of ranks) {
      const biological = nodes.filter(node => node.kind === "taxon" && node.column_rank === rank)
        .sort((a, b) => a.source_order_index - b.source_order_index);
      const indexes = new Map(biological.map((node, index) => [node.id, index]));
      const pathOrder = biological.map(node => node.path).sort(compareBytes);
      const items = biological.map(node => ({key: [indexes.get(node.id), 1, 0, "", node.id], node}));
      for (const node of nodes.filter(item => item.column_rank === rank &&
          (item.kind === "residual" || item.kind === "residual_carry"))) {
        let slot = biological.length;
        if (node.origin_source_id !== ENTRY_ID) {
          const origin = source.byId.get(node.origin_source_id);
          if (!origin) fail("residual origin is unknown");
          const descendants = biological.filter(item => descendant(item, origin.path))
            .map(item => indexes.get(item.id));
          if (descendants.length) slot = Math.max(...descendants) + 1;
          else {
            slot = pathOrder.filter(path => compareBytes(path, origin.path) <= 0).length;
          }
        }
        items.push({key: [slot, 0, SUBTYPE_ORDER[node.subtype], node.lane_id, node.id], node});
      }
      items.sort((a, b) => {
        for (let i = 0; i < a.key.length; i += 1) {
          const left = a.key[i], right = b.key[i];
          if (typeof left === "string") { const result = compareBytes(left, right); if (result) return result; }
          else if (left !== right) return left - right;
        }
        return 0;
      });
      items.forEach((item, index) => { item.node.visual_order = index; });
    }
    const entry = nodes.find(node => node.kind === "entry");
    if (entry) entry.visual_order = 0;
  }

  function canonicalNodes(nodes) {
    return nodes.slice().sort((a, b) => {
      const ca = a.kind === "entry" ? -1 : RANK_INDEX[a.column_rank];
      const cb = b.kind === "entry" ? -1 : RANK_INDEX[b.column_rank];
      return ca - cb || NODE_KINDS[a.kind] - NODE_KINDS[b.kind] ||
        (a.source_order_index == null ? -1 : a.source_order_index) -
        (b.source_order_index == null ? -1 : b.source_order_index) || compareBytes(a.id, b.id);
    });
  }

  function canonicalLinks(links, nodes) {
    const order = new Map(nodes.map((node, index) => [node.id, index]));
    const seen = new Set();
    for (const edge of links) {
      const key = edge.source + "\0" + edge.target + "\0" + edge.kind;
      if (seen.has(key)) fail("duplicate link");
      seen.add(key);
      safeInteger(edge.value, "link value");
      if (!order.has(edge.source) || !order.has(edge.target)) fail("unknown link endpoint");
    }
    return links.slice().sort((a, b) => a.transition_index - b.transition_index ||
      order.get(a.source) - order.get(b.source) || order.get(a.target) - order.get(b.target) ||
      LINK_KINDS[a.kind] - LINK_KINDS[b.kind] || compareBytes(a.id, b.id));
  }

  function assertConservation(source, ranks, nodes, links) {
    const classified = source.classified;
    if (classified === 0) {
      if (nodes.length || links.length) fail("zero-classified view is not empty");
      return {classified: 0,
        column_totals: ranks.map(rank => ({rank, value: 0})),
        transition_totals: ranks.map((rank, index) => ({from: index ? ranks[index - 1] : "ENTRY", to: rank, value: 0})),
        rightmost_flow: 0};
    }
    const transitions = ranks.map((rank, index) => {
      const value = links.filter(edge => edge.transition_index === index)
        .reduce((sum, edge) => sum + edge.value, 0);
      if (value !== classified) fail("transition conservation failed");
      return {from: index ? ranks[index - 1] : "ENTRY", to: rank, value};
    });
    const columns = ranks.map(rank => {
      const value = nodes.filter(node => node.column_rank === rank)
        .reduce((sum, node) => sum + node.value, 0);
      if (value !== classified) fail("column conservation failed");
      return {rank, value};
    });
    const outgoing = new Map();
    const incoming = new Map();
    for (const edge of links) {
      if (!outgoing.has(edge.source)) outgoing.set(edge.source, []);
      if (!incoming.has(edge.target)) incoming.set(edge.target, []);
      outgoing.get(edge.source).push(edge); incoming.get(edge.target).push(edge);
    }
    for (const node of nodes) {
      const edges = outgoing.get(node.id) || [];
      if (node.kind === "taxon" && node.column_rank !== ranks[ranks.length - 1]) {
        if (edges.reduce((sum, edge) => sum + edge.value, 0) !== node.clade) fail("source conservation failed");
      }
      if ((node.kind === "residual" || node.kind === "residual_carry") &&
          node.column_rank !== ranks[ranks.length - 1] &&
          (edges.length !== 1 || edges[0].value !== node.value)) fail("carry conservation failed");
      if (node.kind === "taxon" && node.column_rank !== ranks[0]) {
        const biological = (incoming.get(node.id) || []).filter(edge => edge.kind === "biological")
          .reduce((sum, edge) => sum + edge.value, 0);
        if (biological !== node.clade) fail("target conservation failed");
      }
    }
    return {classified, column_totals: columns, transition_totals: transitions,
      rightmost_flow: columns[columns.length - 1].value};
  }

  function buildView(payload, ranks, maxN) {
    const source = sourceFromPayload(payload);
    const selectedRanks = parseRanks(ranks);
    const selected = selection(source, selectedRanks, maxN);
    if (source.classified === 0) return {nodes: [], links: [],
      conservation: assertConservation(source, selectedRanks, [], [])};
    const nodes = [];
    const links = [];
    const entry = nodeRecord({
      id: ENTRY_ID, kind: "entry", subtype: null, biological: false, carried: false,
      lane_id: null, taxid: null, rank: null, name: "Classified reads", path: null,
      parent_id: null, clade: null, direct: null, status: null, selection_reason: "entry",
      source_order_index: null, origin_source_id: null, origin_target_rank: null,
      column_rank: "ENTRY", value: source.classified, member_count: null, member_paths_sha256: null,
    });
    nodes.push(entry);
    const bioByRank = new Map();
    let residuals = [];
    const first = selectedRanks[0];
    const firstBio = Array.from(selected.retained[first]).map(path =>
      bioNode(source.byPath.get(path), first, selected.reasons.get(path)))
      .sort((a, b) => a.source_order_index - b.source_order_index);
    bioByRank.set(first, firstBio); nodes.push(...firstBio);
    const visibleFirst = selected.retained[first];
    const allFirst = source.byRank[first];
    const hiddenFirst = allFirst.filter(node => !visibleFirst.has(node.path));
    const hiddenMass = hiddenFirst.reduce((sum, node) => sum + node.clade, 0);
    const above = source.classified - allFirst.reduce((sum, node) => sum + node.clade, 0);
    if (above < 0) fail("first rank accounts for too many reads");
    if (hiddenMass > 0) {
      const residual = residualNode("other_hidden", ENTRY_ID, null, first, first, hiddenMass,
        false, hiddenFirst.length, digestMembers(hiddenFirst.map(node => node.path)));
      residuals.push(residual); nodes.push(residual);
      links.push(link(ENTRY_ID, residual.id, "other_hidden", hiddenMass, 0));
    }
    if (above > 0) {
      const residual = residualNode("assigned_above", ENTRY_ID, null, first, first, above,
        false, null, null);
      residuals.push(residual); nodes.push(residual);
      links.push(link(ENTRY_ID, residual.id, "assigned_above", above, 0));
    }
    for (const node of firstBio) links.push(link(ENTRY_ID, node.id, "biological", node.value, 0));
    for (let transition = 1; transition < selectedRanks.length; transition += 1) {
      const sourceRank = selectedRanks[transition - 1];
      const targetRank = selectedRanks[transition];
      const targetBio = Array.from(selected.retained[targetRank]).map(path =>
        bioNode(source.byPath.get(path), targetRank, selected.reasons.get(path)))
        .sort((a, b) => a.source_order_index - b.source_order_index);
      bioByRank.set(targetRank, targetBio); nodes.push(...targetBio);
      const targetByPath = new Map(targetBio.map(node => [node.path, node]));
      const nextResiduals = [];
      for (const sourceNode of bioByRank.get(sourceRank)) {
        const allTargets = source.byRank[targetRank].filter(node => descendant(node, sourceNode.path));
        const visible = allTargets.filter(node => targetByPath.has(node.path));
        const hidden = allTargets.filter(node => !targetByPath.has(node.path));
        const allMass = allTargets.reduce((sum, node) => sum + node.clade, 0);
        const assigned = sourceNode.clade - allMass;
        if (assigned < 0) fail("negative assigned-above mass");
        for (const target of visible) links.push(link(sourceNode.id, targetByPath.get(target.path).id,
          "biological", target.clade, transition));
        const hiddenValue = hidden.reduce((sum, node) => sum + node.clade, 0);
        if (hiddenValue > 0) {
          const residual = residualNode("other_hidden", sourceNode.id, sourceNode.source_order_index,
            targetRank, targetRank, hiddenValue, false, hidden.length,
            digestMembers(hidden.map(node => node.path)));
          nextResiduals.push(residual); nodes.push(residual);
          links.push(link(sourceNode.id, residual.id, "other_hidden", hiddenValue, transition));
        }
        if (assigned > 0) {
          const residual = residualNode("assigned_above", sourceNode.id, sourceNode.source_order_index,
            targetRank, targetRank, assigned, false, null, null);
          nextResiduals.push(residual); nodes.push(residual);
          links.push(link(sourceNode.id, residual.id, "assigned_above", assigned, transition));
        }
      }
      for (const prior of residuals) {
        const carry = residualNode(prior.subtype, prior.origin_source_id, prior.source_order_index,
          prior.origin_target_rank, targetRank, prior.value, true, prior.member_count,
          prior.member_paths_sha256);
        nextResiduals.push(carry); nodes.push(carry);
        links.push(link(prior.id, carry.id, "carry", prior.value, transition));
      }
      residuals = nextResiduals;
    }
    assignVisual(nodes, source, selectedRanks);
    const orderedNodes = canonicalNodes(nodes);
    const orderedLinks = canonicalLinks(links, orderedNodes);
    const conservation = assertConservation(source, selectedRanks, orderedNodes, orderedLinks);
    return {nodes: orderedNodes, links: orderedLinks, conservation};
  }

  function decodePayload(encoded) {
    const binary = atob(encoded.trim());
    const bytes = Uint8Array.from(binary, character => character.charCodeAt(0));
    return JSON.parse(new TextDecoder("utf-8").decode(bytes));
  }

  function initBrowser() {
    const payloadElement = document.getElementById("wf16s-payload");
    if (!payloadElement) return;
    const payload = decodePayload(payloadElement.textContent);
    const defaultRanks = payload.defaults.ranks.slice();
    const defaultN = payload.defaults.max_taxa_per_rank;
    const rankControl = document.getElementById("ranks");
    const maxControl = document.getElementById("max-n");
    const status = document.getElementById("status");
    const svg = document.getElementById("sankey");
    const search = document.getElementById("search");
    const percent = document.getElementById("percent");
    const pin = {id: null};
    RANKS.forEach(rank => {
      const option = document.createElement("option"); option.value = rank; option.textContent = rank;
      option.selected = defaultRanks.indexOf(rank) >= 0; rankControl.appendChild(option);
    });
    maxControl.value = String(defaultN);
    const rawBytes = Uint8Array.from(atob(payloadElement.textContent.trim()), character => character.charCodeAt(0));
    const jsonUrl = URL.createObjectURL(new Blob([rawBytes], {type: "application/json"}));
    document.getElementById("download-json").href = jsonUrl;

    function selectedRanks() { return Array.from(rankControl.selectedOptions).map(option => option.value); }
    function svgElement(name) {
      const namespace = "ht" + "tp://www.w3.org/2000/svg";
      return document.createElementNS(namespace, name);
    }
    function draw(view, verified) {
      while (svg.firstChild) svg.removeChild(svg.firstChild);
      const ranks = ["ENTRY"].concat(selectedRanks());
      const columns = new Map(ranks.map(rank => [rank, []]));
      view.nodes.forEach(node => { if (columns.has(node.column_rank)) columns.get(node.column_rank).push(node); });
      const positions = new Map();
      const maxValue = Math.max(1, payload.totals.classified);
      const width = 1200, height = 640, columnWidth = width / (ranks.length + 1);
      ranks.forEach((rank, column) => {
        const items = columns.get(rank).slice().sort((a, b) => a.visual_order - b.visual_order);
        const gap = Math.min(18, (height - 40) / Math.max(1, items.length * 2));
        let cursor = 20;
        items.forEach(node => {
          const nodeHeight = Math.max(8, Math.min(70, 8 + 70 * node.value / maxValue));
          const x = 20 + column * columnWidth;
          const y = cursor;
          positions.set(node.id, {x, y, width: 13, height: nodeHeight});
          cursor += nodeHeight + gap;
        });
      });
      view.links.forEach(edge => {
        const source = positions.get(edge.source), target = positions.get(edge.target);
        if (!source || !target) return;
        const path = svgElement("path");
        const x1 = source.x + source.width, x2 = target.x;
        const y1 = source.y + source.height / 2, y2 = target.y + target.height / 2;
        const mid = (x1 + x2) / 2;
        path.setAttribute("d", `M ${x1} ${y1} C ${mid} ${y1}, ${mid} ${y2}, ${x2} ${y2}`);
        path.setAttribute("class", "link");
        path.setAttribute("stroke", edge.kind === "biological" ? "#4682b4" : "#9e9e9e");
        path.setAttribute("stroke-width", String(Math.max(1, 18 * edge.value / maxValue)));
        path.setAttribute("aria-label", `${edge.kind} ${edge.value}`);
        svg.appendChild(path);
      });
      const needle = search.value.toLocaleLowerCase();
      view.nodes.forEach(node => {
        const pos = positions.get(node.id); if (!pos) return;
        const group = svgElement("g"); group.setAttribute("class", `node ${node.kind} ${node.status || ""}`);
        group.setAttribute("tabindex", "0"); group.setAttribute("role", "button");
        group.dataset.nodeId = node.id;
        const isMatch = !needle || String(node.name || "").toLocaleLowerCase().includes(needle) ||
          String(node.path || "").toLocaleLowerCase().includes(needle);
        if (!isMatch) group.setAttribute("opacity", "0.25");
        if (pin.id && pin.id !== node.id) {
          const related = view.links.some(edge => (edge.source === pin.id && edge.target === node.id) ||
            (edge.target === pin.id && edge.source === node.id));
          if (!related) group.setAttribute("opacity", "0.2");
        }
        const rect = svgElement("rect"); rect.setAttribute("x", pos.x); rect.setAttribute("y", pos.y);
        rect.setAttribute("width", pos.width); rect.setAttribute("height", pos.height);
        group.appendChild(rect);
        const label = svgElement("text"); label.setAttribute("x", pos.x + 18); label.setAttribute("y", pos.y + pos.height / 2);
        const amount = percent.checked ? `${(100 * node.value / Math.max(1, payload.totals.classified)).toFixed(2)}%` : String(node.value);
        label.textContent = `${node.name || node.kind} (${amount})`;
        group.appendChild(label);
        group.addEventListener("click", () => { pin.id = pin.id === node.id ? null : node.id; draw(view, verified); });
        group.addEventListener("keydown", event => { if (event.key === "Enter" || event.key === " ") { event.preventDefault(); pin.id = pin.id === node.id ? null : node.id; draw(view, verified); } });
        svg.appendChild(group);
      });
      status.textContent = verified ? "Verified default view" : "Interactive transient view";
    }

    function render(verified) {
      let ranks = selectedRanks();
      if (ranks.length < 2) { ranks = defaultRanks.slice(); rankControl.querySelectorAll("option").forEach(option => { option.selected = ranks.includes(option.value); }); }
      let value = Number(maxControl.value);
      if (!Number.isSafeInteger(value) || value < 1 || value > 100) value = defaultN;
      try { draw(buildView(payload, ranks, value), verified); }
      catch (error) { status.textContent = `Invalid transient view: ${error.message}`; }
    }
    rankControl.addEventListener("change", () => render(false));
    maxControl.addEventListener("change", () => render(false));
    search.addEventListener("input", () => render(false));
    percent.addEventListener("change", () => render(false));
    document.getElementById("reset").addEventListener("click", () => {
      rankControl.querySelectorAll("option").forEach(option => { option.selected = defaultRanks.includes(option.value); });
      maxControl.value = String(defaultN); search.value = ""; percent.checked = false; pin.id = null; render(true);
    });
    document.getElementById("download-svg").addEventListener("click", () => {
      const blob = new Blob([new XMLSerializer().serializeToString(svg)], {type: "image/svg+xml"});
      const url = URL.createObjectURL(blob); const link = document.createElement("a");
      link.href = url; link.download = `${payload.sample_id}.sankey.svg`; link.click(); URL.revokeObjectURL(url);
    });
    render(true);
  }

  const api = {buildView, parseRanks, sourceFromPayload, sha256, MAX_SAFE_INTEGER, COUNT_MODEL};
  if (typeof module !== "undefined" && module.exports) module.exports = api;
  else {
    root.WF16STaxonomySankey = api;
    if (document.readyState === "loading") document.addEventListener("DOMContentLoaded", initBrowser);
    else initBrowser();
  }
}(typeof globalThis === "undefined" ? this : globalThis));
