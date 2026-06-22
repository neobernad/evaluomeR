export interface PCACoord {
  pc1: number
  pc2: number
}

export interface Sample {
  id: string
  cancerType: string
  metrics: Record<string, number>
}

export interface OptimalKRow {
  metric: string
  stabilityMaxK: number
  qualityMaxK: number
  globalOptimalK: number
}

export interface KSummaryEntry {
  avgStability: number
  avgSilhouette: number
  composite: number
}

export type MetricMap = Record<string, number>
export type ByK<T> = Record<string, T>

export interface DemoData {
  meta: {
    dataset: string
    nSamples: number
    nCancerTypes: number
    cancerTypes: string[]
    metrics: string[]
    chartMetrics: string[]
    kRange: [number, number]
    bs: number
    allMetrics?: boolean
    previewMetrics?: string[]
  }
  samples: Sample[]
  stability: ByK<MetricMap>
  quality: {
    silhouette: ByK<MetricMap>
    ch: ByK<MetricMap>
  }
  optimalK: number
  optimalKPerMetric: Record<string, number>
  optimalKDetail: OptimalKRow[]
  kSummary: Record<string, KSummaryEntry>
  clusters: ByK<number[]>
  pca: {
    coords: PCACoord[]
    varExplained: [number, number]
  }
}
