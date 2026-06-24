#!/usr/bin/env node
/**
 * Augment existing demo JSON with ari, sampleSilhouette, stabilitySD
 * (fallback when full R export is unavailable).
 */
import { readFileSync, writeFileSync } from 'fs'
import { join, dirname } from 'path'
import { fileURLToPath } from 'url'

const __dirname = dirname(fileURLToPath(import.meta.url))
const dataDir = join(__dirname, '../ui/src/data')

function euclidean(a, b) {
  let s = 0
  for (let i = 0; i < a.length; i++) s += (a[i] - b[i]) ** 2
  return Math.sqrt(s)
}

function adjustedRandIndex(labels1, labels2) {
  const n = labels1.length
  const map1 = new Map()
  const map2 = new Map()
  for (let i = 0; i < n; i++) {
    map1.set(labels1[i], (map1.get(labels1[i]) ?? 0) + 1)
    map2.set(labels2[i], (map2.get(labels2[i]) ?? 0) + 1)
  }
  let sumComb = 0
  for (const c of map1.values()) sumComb += (c * (c - 1)) / 2
  for (const c of map2.values()) sumComb += (c * (c - 1)) / 2

  const contingency = new Map()
  for (let i = 0; i < n; i++) {
    const key = `${labels1[i]}|${labels2[i]}`
    contingency.set(key, (contingency.get(key) ?? 0) + 1)
  }
  let index = 0
  for (const nij of contingency.values()) index += (nij * (nij - 1)) / 2

  const expected =
    sumComb * [...map2.values()].reduce((s, c) => s + (c * (c - 1)) / 2, 0) /
    ((n * (n - 1)) / 2)
  const maxIndex =
    (sumComb + [...map2.values()].reduce((s, c) => s + (c * (c - 1)) / 2, 0)) / 2
  if (maxIndex === expected) return 1
  return (index - expected) / (maxIndex - expected)
}

function silhouetteWidths(samples, clusters, metricKeys) {
  const n = samples.length
  const vectors = samples.map((s) => metricKeys.map((m) => s.metrics[m] ?? 0))

  const byCluster = new Map()
  for (let i = 0; i < n; i++) {
    const c = clusters[i]
    if (!byCluster.has(c)) byCluster.set(c, [])
    byCluster.get(c).push(i)
  }

  const widths = new Array(n)
  for (let i = 0; i < n; i++) {
    const ci = clusters[i]
    const same = byCluster.get(ci).filter((j) => j !== i)
    let a = 0
    if (same.length === 0) {
      a = 0
    } else {
      a = same.reduce((s, j) => s + euclidean(vectors[i], vectors[j]), 0) / same.length
    }
    let b = Infinity
    for (const [c, members] of byCluster) {
      if (c === ci) continue
      const meanDist =
        members.reduce((s, j) => s + euclidean(vectors[i], vectors[j]), 0) / members.length
      b = Math.min(b, meanDist)
    }
    if (!isFinite(b)) b = 0
    const m = Math.max(a, b)
    widths[i] = m === 0 ? 0 : (b - a) / m
  }
  return widths.map((w) => Math.round(w * 1e6) / 1e6)
}

function std(values) {
  const m = values.reduce((a, b) => a + b, 0) / values.length
  const v = values.reduce((s, x) => s + (x - m) ** 2, 0) / values.length
  return Math.round(Math.sqrt(v) * 1e6) / 1e6
}

function augment(demo) {
  const metricKeys = demo.meta.metrics
  const typeToInt = new Map()
  let ti = 0
  const trueLabels = demo.samples.map((s) => {
    if (!typeToInt.has(s.cancerType)) typeToInt.set(s.cancerType, ti++)
    return typeToInt.get(s.cancerType)
  })

  const stabilitySD = {}
  const sampleSilhouette = {}

  for (const k of Object.keys(demo.clusters)) {
    const clusters = demo.clusters[k]
    demo.kSummary[k].ari = Math.round(adjustedRandIndex(trueLabels, clusters) * 1e6) / 1e6
    sampleSilhouette[k] = silhouetteWidths(demo.samples, clusters, metricKeys)

    const stabVals = Object.values(demo.stability[k] ?? {}).filter((v) => typeof v === 'number')
    stabilitySD[k] = stabVals.length > 0 ? std(stabVals) : 0
  }

  demo.stabilitySD = stabilitySD
  demo.sampleSilhouette = sampleSilhouette
  return demo
}

for (const file of ['nci60.json', 'golub.json']) {
  const path = join(dataDir, file)
  const demo = JSON.parse(readFileSync(path, 'utf8'))
  augment(demo)
  writeFileSync(path, JSON.stringify(demo))
  console.log('Augmented', file)
}
