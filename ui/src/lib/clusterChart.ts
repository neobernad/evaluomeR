import type { DemoData, Sample } from '@/types/demo'
import { getCancerColor } from '@/lib/cancerColors'
import { computeEllipsePoints } from '@/lib/ellipse'

const CLUSTER_COLORS = [
  '#3b82f6', '#a78bfa', '#22c55e', '#f59e0b',
  '#f43f5e', '#06b6d4', '#fb923c', '#84cc16',
]

function computeClusterSphere(points: [number, number, number][], sigmaScale = 1.5) {
  if (points.length === 0) return null

  const n = points.length
  const cx = points.reduce((sum, point) => sum + point[0], 0) / n
  const cy = points.reduce((sum, point) => sum + point[1], 0) / n
  const cz = points.reduce((sum, point) => sum + point[2], 0) / n

  const meanRadius =
    points.reduce((sum, point) => {
      const dx = point[0] - cx
      const dy = point[1] - cy
      const dz = point[2] - cz
      return sum + Math.sqrt(dx * dx + dy * dy + dz * dz)
    }, 0) / n

  return {
    cx,
    cy,
    cz,
    radius: Math.max(meanRadius * sigmaScale, 0.08),
  }
}

function buildFibonacciSphere(
  sphere: { cx: number; cy: number; cz: number; radius: number },
  color: string,
  nPoints = 400,
) {
  const { cx, cy, cz, radius: r } = sphere
  const phi = Math.PI * (3 - Math.sqrt(5))
  const pts = Array.from({ length: nPoints }, (_, i) => {
    const y = 1 - (i / (nPoints - 1)) * 2
    const ry = Math.sqrt(Math.max(0, 1 - y * y))
    const theta = phi * i
    return {
      value: [
        cx + r * Math.cos(theta) * ry,
        cy + r * y,
        cz + r * Math.sin(theta) * ry,
      ] as [number, number, number],
    }
  })

  return [
    {
      type: 'scatter3D',
      coordinateSystem: 'cartesian3D',
      data: pts,
      symbolSize: 3,
      itemStyle: { color, opacity: 0.45 },
      silent: true,
    },
  ]
}

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

export function buildPCACluster3DOption(params: {
  samples: Sample[]
  clusters: DemoData['clusters']
  k: number
  pca: DemoData['pca']
}) {
  const { samples, clusters, k, pca } = params
  const [var1, var2, var3] = pca.varExplained
  const xLabel = `PC1 (${(var1 * 100).toFixed(1)}%)`
  const yLabel = `PC2 (${(var2 * 100).toFixed(1)}%)`
  const zLabel = `PC3 (${((var3 ?? 0) * 100).toFixed(1)}%)`

  const clusterAssignments = clusters[String(k)] ?? []
  const clusterIds = [...new Set(clusterAssignments)].sort((a, b) => a - b)
  const cancerTypes = [...new Set(samples.map((s) => s.cancerType))].sort()

  const axisStyle = {
    nameTextStyle: { color: '#94a3b8', fontSize: 11 },
    axisLabel: { color: '#94a3b8', fontSize: 10 },
    axisLine: { lineStyle: { color: '#475569' } },
    splitLine: { lineStyle: { color: '#1e293b' } },
  }

  const sphereSeries = clusterIds.flatMap((cid, idx) => {
    const color = CLUSTER_COLORS[idx % CLUSTER_COLORS.length]
    const points: [number, number, number][] = samples
      .map((_, i) => {
        if (clusterAssignments[i] !== cid) return null
        const coord = pca.coords[i]
        if (coord.pc3 === undefined) return null
        return [coord.pc1, coord.pc2, coord.pc3] as [number, number, number]
      })
      .filter((point): point is [number, number, number] => point !== null)

    const sphere = computeClusterSphere(points)
    return sphere ? buildFibonacciSphere(sphere, color) : []
  })

  const tissueScatterSeries = cancerTypes.map((ct) => ({
    name: ct,
    type: 'scatter3D',
    data: samples
      .map((sample, i) => {
        if (sample.cancerType !== ct) return null
        const coord = pca.coords[i]
        if (coord.pc3 === undefined) return null
        return {
          name: sample.id,
          value: [coord.pc1, coord.pc2, coord.pc3] as [number, number, number],
          tissue: sample.cancerType,
          cluster: clusterAssignments[i],
        }
      })
      .filter((point): point is NonNullable<typeof point> => point !== null),
    symbolSize: 9,
    itemStyle: { color: getCancerColor(ct), opacity: 0.95 },
    z: 3,
  }))

  return {
    backgroundColor: 'transparent',
    tooltip: {
      trigger: 'item',
      formatter: (params: {
        data?: { name?: string; tissue?: string; cluster?: number }
        seriesName?: string
      }) => {
        const data = params.data
        if (!data?.name) return ''
        return `<b>${data.name}</b><br/>Tissue: ${data.tissue ?? '—'}<br/>Cluster: ${data.cluster ?? params.seriesName}`
      },
    },
    legend: {
      data: cancerTypes,
      textStyle: { color: '#94a3b8', fontSize: 11 },
      top: 0,
      type: 'scroll',
    },
    xAxis3D: { type: 'value', name: xLabel, ...axisStyle },
    yAxis3D: { type: 'value', name: yLabel, ...axisStyle },
    zAxis3D: { type: 'value', name: zLabel, ...axisStyle },
    grid3D: {
      boxWidth: 120,
      boxDepth: 120,
      boxHeight: 120,
      viewControl: {
        projection: 'perspective',
        distance: 220,
        alpha: 80,
        beta: 0,
      },
      light: {
        main: { intensity: 1.1, shadow: false },
        ambient: { intensity: 0.45 },
      },
    },
    series: [...sphereSeries, ...tissueScatterSeries],
  }
}
