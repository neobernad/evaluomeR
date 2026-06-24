import { useMemo } from 'react'
import ReactECharts from 'echarts-for-react'
import type { Sample } from '@/types/demo'
import { formatMetricLabel } from '@/lib/metricLabels'

const CLUSTER_COLORS = [
  '#3b82f6', '#a78bfa', '#22c55e', '#f59e0b',
  '#f43f5e', '#06b6d4', '#fb923c', '#84cc16',
]

interface ClusterRadarProps {
  samples: Sample[]
  clusters: number[]
  chartMetrics: string[]
}

export function ClusterRadar({ samples, clusters, chartMetrics }: ClusterRadarProps) {
  const option = useMemo(() => {
    const clusterIds = [...new Set(clusters)].sort((a, b) => a - b)

    const ranges = chartMetrics.map((m) => {
      const vals = samples.map((s) => s.metrics[m] ?? 0)
      return { min: Math.min(...vals), max: Math.max(...vals) }
    })

    const normalize = (metric: string, value: number) => {
      const idx = chartMetrics.indexOf(metric)
      const { min, max } = ranges[idx]
      if (max === min) return 0.5
      return (value - min) / (max - min)
    }

    const indicators = chartMetrics.map((m) => ({
      name: formatMetricLabel(m),
      max: 1,
    }))

    const seriesData = clusterIds.map((cid, idx) => {
      const members = samples.filter((_, i) => clusters[i] === cid)
      const means = chartMetrics.map((m) => {
        const avg = members.reduce((s, sample) => s + (sample.metrics[m] ?? 0), 0) / members.length
        return normalize(m, avg)
      })
      return {
        name: `Cluster ${cid}`,
        value: means,
        lineStyle: { color: CLUSTER_COLORS[idx % CLUSTER_COLORS.length] },
        areaStyle: { color: CLUSTER_COLORS[idx % CLUSTER_COLORS.length], opacity: 0.15 },
        itemStyle: { color: CLUSTER_COLORS[idx % CLUSTER_COLORS.length] },
      }
    })

    return {
      backgroundColor: 'transparent',
      tooltip: { trigger: 'item' },
      legend: {
        data: clusterIds.map((cid) => `Cluster ${cid}`),
        textStyle: { color: '#94a3b8', fontSize: 11 },
        top: 0,
      },
      radar: {
        indicator: indicators,
        center: ['50%', '55%'],
        radius: '62%',
        axisName: { color: '#94a3b8', fontSize: 9 },
        splitLine: { lineStyle: { color: '#1e293b' } },
        splitArea: { areaStyle: { color: ['rgba(15,23,42,0.3)', 'rgba(15,23,42,0.5)'] } },
        axisLine: { lineStyle: { color: '#334155' } },
      },
      series: [
        {
          type: 'radar',
          data: seriesData,
        },
      ],
    }
  }, [samples, clusters, chartMetrics])

  return (
    <div>
      <p className="mb-2 text-xs text-slate-500">
        Normalized mean metric profile per cluster — higher values indicate stronger presence of that metric.
      </p>
      <ReactECharts option={option} style={{ height: 400 }} opts={{ renderer: 'canvas' }} notMerge />
    </div>
  )
}
