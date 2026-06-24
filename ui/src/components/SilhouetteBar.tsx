import { useMemo } from 'react'
import ReactECharts from 'echarts-for-react'
import type { Sample } from '@/types/demo'

const CLUSTER_COLORS = [
  '#3b82f6', '#a78bfa', '#22c55e', '#f59e0b',
  '#f43f5e', '#06b6d4', '#fb923c', '#84cc16',
]

interface SilhouetteBarProps {
  samples: Sample[]
  clusters: number[]
  silhouetteWidths: number[]
}

export function SilhouetteBar({ samples, clusters, silhouetteWidths }: SilhouetteBarProps) {
  const option = useMemo(() => {
    const entries = samples.map((s, i) => ({
      id: s.id,
      cluster: clusters[i] ?? 0,
      width: silhouetteWidths[i] ?? 0,
    }))

    entries.sort((a, b) => {
      if (a.cluster !== b.cluster) return a.cluster - b.cluster
      return b.width - a.width
    })

    const labels = entries.map((e) => e.id)
    const data = entries.map((e) => ({
      value: e.width,
      itemStyle: {
        color: e.width < 0 ? '#ef4444' : CLUSTER_COLORS[e.cluster % CLUSTER_COLORS.length],
      },
    }))

    const minW = Math.min(...entries.map((e) => e.width), 0)
    const maxW = Math.max(...entries.map((e) => e.width), 1)

    return {
      backgroundColor: 'transparent',
      tooltip: {
        trigger: 'axis',
        axisPointer: { type: 'shadow' },
        formatter: (params: { name: string; value: number }[]) => {
          const p = params[0]
          if (!p) return ''
          const entry = entries.find((e) => e.id === p.name)
          return `<b>${p.name}</b><br/>Cluster: ${entry?.cluster ?? '—'}<br/>Silhouette: ${p.value.toFixed(3)}`
        },
      },
      grid: { left: '18%', right: '6%', bottom: '8%', top: '4%' },
      xAxis: {
        type: 'value',
        min: Math.min(minW - 0.05, -0.2),
        max: Math.max(maxW + 0.05, 1),
        axisLabel: { color: '#94a3b8', fontSize: 10 },
        splitLine: { lineStyle: { color: '#1e293b' } },
      },
      yAxis: {
        type: 'category',
        data: labels,
        axisLabel: { color: '#94a3b8', fontSize: 8 },
        inverse: true,
      },
      series: [
        {
          type: 'bar',
          data,
          barMaxWidth: 12,
        },
      ],
    }
  }, [samples, clusters, silhouetteWidths])

  return (
    <ReactECharts option={option} style={{ height: Math.max(280, samples.length * 14) }} opts={{ renderer: 'canvas' }} notMerge />
  )
}
