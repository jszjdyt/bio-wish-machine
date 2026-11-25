/**
 * System-prompt guidance for the virtual-tumor tools. Kept short: it tells the
 * agent when the tools apply and what the numbers mean, and leaves the science
 * to the model's own knowledge.
 * @module dsh-nsclc-vt/prompt
 */

/**
 * @param {{modelName:string, nodeCount:number, phenotypes:string[], backgrounds:string[]}} info
 * @returns {string}
 */
export function buildGuidance(info) {
    return [
        '## NSCLC virtual tumor (nsclc_vt plugin)',
        '',
        `An in-process qualitative-network simulation of NSCLC adenocarcinoma signalling — the "virtual tumor" of Clarke et al., Cancer Research 2026 (doi:10.1158/0008-5472.CAN-25-2476), running the BioModelAnalyzer VMCAI stability engine. The loaded model is "${info.modelName}" with ${info.nodeCount} signalling nodes; the readouts are ${info.phenotypes.join(' and ')}.`,
        '',
        'Reach for these tools when the user asks what a genetic background or a drug target does to NSCLC proliferation/death, which drug targets combine well, or how a gene is wired into the signalling network:',
        '- `nsclc_model` — discover node names, gene symbols and the bundled cell-line backgrounds. Start here if unsure of exact names.',
        '- `nsclc_node` — one node\'s target function, curated activators/inhibitors with citations, and its wiring.',
        '- `nsclc_predict` — steady state for one genotype plus interventions; reports phenotypes and which nodes moved most.',
        '- `nsclc_screen` — mono and pairwise perturbation screens with the manuscript\'s additivity score; use `candidates` to bound pair screens.',
        '',
        'Conventions, so results are read correctly:',
        '- A "knockout"/"drug" is the node forced to level 0; higher levels mean more activity. `interventions` entries are `NODE` (knockout) or `NODE:level`.',
        '- `Difference = mean(Proliferation) - mean(Death)`; lower is therapeutically better.',
        '- `additive = Difference(pair) - min(Difference of the two singles)`; **negative means the pair beats its better single component (synergy)**, positive means antagonism. This is the manuscript\'s own definition.',
        '- Levels are integer intervals, not concentrations. `uncertainty = (hi - lo)/range` is 0 only when the analysis pinned the node to one value; a non-zero value means the network did not prove a unique steady state for that node.',
        '- Every result is a model prediction on a qualitative abstraction, useful for ranking hypotheses and designing experiments. Say so when reporting; do not present it as clinical guidance.',
    ].join('\n');
}
