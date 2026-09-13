import assert from "node:assert/strict";
import {createRequire} from "node:module";

const require = createRequire(import.meta.url);
const client = require("../analysis/utils/taxonomy_sankey_client.js");

function sourceNode(rank, name, parentPath, clade, direct, index) {
  const path = parentPath ? `${parentPath};${name}` : name;
  return {
    id: `taxon:${rank}:${client.sha256(path)}`,
    taxid: String(100 + index),
    rank,
    name,
    path,
    parent_id: parentPath ? `taxon:${parentPath === "Alpha" ? "D" : parentPath === "Beta" ? "D" : "K"}:${client.sha256(parentPath)}` : null,
    clade,
    direct,
    status: "resolved",
    source_order_index: index,
  };
}

const nodes = [];
const add = (rank, name, parent, clade, direct) => {
  const node = sourceNode(rank, name, parent, clade, direct, nodes.length);
  nodes.push(node);
  return node.path;
};
add("D", "Alpha", "", 60, 0);
add("K", "A1", "Alpha", 40, 0);
add("P", "A1a", "Alpha;A1", 20, 20);
add("P", "A1b", "Alpha;A1", 20, 20);
add("K", "A2", "Alpha", 20, 0);
add("P", "A2a", "Alpha;A2", 20, 20);
add("D", "Beta", "", 25, 0);
add("K", "B1", "Beta", 25, 0);
add("P", "B1a", "Beta;B1", 25, 25);
add("D", "Épsilon", "", 15, 15);

const payload = {
  totals: {total: 100, classified: 100, unclassified: 0},
  source_nodes: nodes,
};

assert.equal(client.sha256("abc"), "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad");
for (const ranks of [["D", "K"], ["D", "P"], ["D", "K", "P"], ["K", "P"]]) {
  for (const maxN of [1, 2, 3]) {
    const view = client.buildView(payload, ranks, maxN);
    assert.equal(view.conservation.rightmost_flow, 100);
    assert.deepEqual(view.conservation.column_totals,
      ranks.map(rank => ({rank, value: 100})));
    assert.ok(view.nodes.every(node => Number.isSafeInteger(node.value) && node.value > 0));
    assert.ok(view.links.every(link => link.value > 0));
    const ids = new Set(view.nodes.map(node => node.id));
    assert.equal(ids.size, view.nodes.length);
    for (const link of view.links) {
      assert.ok(ids.has(link.source));
      assert.ok(ids.has(link.target));
    }
  }
}

const zero = client.buildView({totals: {total: 10, classified: 0, unclassified: 10}, source_nodes: []}, ["D", "K"], 10);
assert.deepEqual(zero.nodes, []);
assert.deepEqual(zero.links, []);
assert.throws(() => client.buildView(payload, ["D", "K"], 0), /invalid max N/);
assert.throws(() => client.buildView(payload, ["K", "D"], 2), /canonical/);

console.log("taxonomy Sankey client tests passed");
