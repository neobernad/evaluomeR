// Shape of ui/src/data/*.json, written by tools/export-demo-data.R. Keep the two in sync.

export interface KSummaryEntry {
  avgStability: number
  avgSilhouette: number
  composite: number
  ari?: number
}

export type ByK<T> = Record<string, T>

export interface DemoData {
  meta: {
    dataset: string
    nSamples: number
    nMetrics: number
    cancerTypes: string[]
    bs: number
  }
  samples: { id: string; cancerType: string }[]
  optimalK: number
  kSummary: ByK<KSummaryEntry>
  clusters: ByK<number[]>
  pca: {
    coords: { pc1: number; pc2: number }[]
    varExplained: number[]
  }
}
