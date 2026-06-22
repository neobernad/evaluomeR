import type { DemoData, Sample } from '@/types/demo'
import { getCancerColor } from '@/lib/cancerColors'
import { computeEllipsePoints } from '@/lib/ellipse'

const CLUSTER_COLORS = [
  '#3b82f6', '#a78bfa', '#22c55e', '#f59e0b',
  '#f43f5e', '#06b6d4', '#fb923c', '#84cc16',
]

export function buildPCAClusterOption(params: {
  title?: string
  titleColor?: string
  samples: Sample[]
  clusters: DemoData['clusters']
  k: number
  pca: DemoData['pca']
  height?: number
  compact?: boolean
}) {
  const { title, titleColor, samples, clusters, k, pca, compact = false } = params
  const [var1, var2] = pca.varExplained
  const xLabel = `PC1 (${(var1 * 100).toFixed(1)}%)`
  const yLabel = `PC2 (${(var2 * 100).toFixed(1)}%)`

  const clusterAssignments = clusters[String(k)] ?? []
  const clusterIds = [...new Set(clusterAssignments)].sort((a, b) => a - b)
  const cancerTypes = [...new Set(samples.map((s) => s.cancerType))].sort()

  const ellipseSeries = clusterIds.map((cid, idx) => {
    const points: [number, number][] = samples
      .map((_, i) => {
        if (clusterAssignments[i] !== cid) return null
        const coord = pca.coords[i]
        return [coord.pc1, coord.pc2] as [number, number]
      })
      .filter((p): p is [number, number] => p !== null)

    const color = CLUSTER_COLORS[idx % CLUSTER_COLORS.length]
    const ellipse = computeEllipsePoints(points)

    return {
      name: `Cluster ${cid}`,
      type: 'line',
      data: ellipse,
      showSymbol: false,
      smooth: true,
      lineStyle: { color, width: 2, opacity: 0.55 },
      areaStyle: { color, opacity: 0.1 },
      silent: true,
      z: 1,
      legendHoverLink: false,
    }
  })

  const scatterSeries = cancerTypes.map((ct) => ({
    name: ct,
    type: 'scatter',
    symbolSize: compact ? 9 : 11,
    itemStyle: { color: getCancerColor(ct) },
    z: 3,
    data: samples
      .map((s, i) => ({
        s,
        cluster: clusterAssignments[i] ?? 0,
        coord: pca.coords[i],
      }))
      .filter(({ s }) => s.cancerType === ct)
      .map(({ s, cluster, coord }) => ({
        name: s.id,
        value: [coord.pc1, coord.pc2],
        cluster,
      })),
  }))

  return {
    backgroundColor: 'transparent',
    title: title
      ? {
          text: title,
          left: 'center',
          textStyle: { color: titleColor ?? '#e2e8f0', fontSize: 13, fontWeight: 700 },
        }
      : undefined,
    tooltip: {
      trigger: 'item',
      formatter: (p: { data: { name: string; cluster: number }; seriesName: string }) =>
        `<b>${p.data.name}</b><br/>Tissue: ${p.seriesName}<br/>Cluster: ${p.data.cluster}`,
    },
    legend: {
      data: cancerTypes,
      textStyle: { color: '#94a3b8', fontSize: compact ? 10 : 11 },
      top: title ? 24 : 0,
      type: 'scroll',
    },
    grid: {
      left: '10%',
      right: '4%',
      bottom: '12%',
      top: title ? '30%' : compact ? '22%' : '18%',
    },
    xAxis: {
      type: 'value',
      name: xLabel,
      nameTextStyle: { color: '#94a3b8', fontSize: compact ? 10 : 12 },
      axisLabel: { color: '#94a3b8', fontSize: compact ? 9 : 11 },
      splitLine: { lineStyle: { color: '#1e293b' } },
    },
    yAxis: {
      type: 'value',
      name: yLabel,
      nameTextStyle: { color: '#94a3b8', fontSize: compact ? 10 : 12 },
      axisLabel: { color: '#94a3b8', fontSize: compact ? 9 : 11 },
      splitLine: { lineStyle: { color: '#1e293b' } },
    },
    series: [...ellipseSeries, ...scatterSeries],
  }
}
