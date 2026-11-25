/**
 * Worker entry for large screens. Receives a chunk of jobs and the model
 * coordinates, rebuilds the analyst in-thread (model parsing is ~10 ms) and
 * returns plain result objects, which are structured-clone friendly.
 * @module dsh-nsclc-vt/screen-worker
 */
import { parentPort, workerData } from 'node:worker_threads';
import { createAnalyst } from './analysis.js';
import { runJobsInline } from './screen.js';

try {
    const analyst = createAnalyst({
        modelPath: workerData.modelPath,
        phenotypes: workerData.phenotypes,
    });
    parentPort.postMessage({ ok: true, results: runJobsInline(analyst, workerData.jobs) });
} catch (error) {
    parentPort.postMessage({ ok: false, error: error instanceof Error ? error.message : String(error) });
}
