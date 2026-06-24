import { useMemo } from 'react'
import ReactECharts from 'echarts-for-react'
import type { KSummaryEntry } from '@/types/demo'

interface KParetoChartProps {
  kSummary: Record<string, KSummaryEntry>
  optimalK: number
  currentK: number
}

export function KParetoChart({ kSummary, optimalK, currentK }: KParetoChartProps) {
  const kValues = Object.keys(kSummary).sort((a, b) => Number(a) - Number(b))

  const option = useMemo(() => {
    const points = kValues.map((kv) => {
      const entry = kSummary[kv]
      const kNum = Number(kv)
      const isOptimal = kNum === optimalK
      const isCurrent = kNum === currentK
      return {
        name: `k = ${kv}`,
        value: [entry.avgStability, entry.avgSilhouette],
        symbolSize: isOptimal ? 18 : isCurrent ? 14 : 10,
        itemStyle: {
          color: isOptimal ? '#22c55e' : isCurrent ? '#3b82f6' : '#94a3b8',
          borderColor: isOptimal ? '#4ade80' : isCurrent ? '#60a5fa' : '#64748b',
          borderWidth: isOptimal || isCurrent ? 2 : 0,
        },
        label: {
          show: true,
          formatter: `k=${kv}`,
          position: 'top',
          color: '#94a3b8',
          fontSize: 10,
        },
      }
    })

    const maxVal = Math.max(
      ...kValues.flatMap((kv) => [kSummary[kv].avgStability, kSummary[kv].avgSilhouette]),
      1,
    )

    return {
      backgroundColor: 'transparent',
      tooltip: {
        trigger: 'item',
        formatter: (p: { name: string; value: [number, number] }) =>
          `<b>${p.name}</b><br/>Stability: ${p.value[0].toFixed(3)}<br/>Silhouette: ${p.value[1].toFixed(3)}`,
      },
      grid: { left: '12%', right: '8%', bottom: '12%', top: '12%' },
      xAxis: {
        type: 'value',
        name: 'Avg stability',
        min: 0,
        max: maxVal * 1.05,
        nameTextStyle: { color: '#94a3b8' },
        axisLabel: { color: '#94a3b8' },
        splitLine: { lineStyle: { color: '#1e293b' } },
      },
      yAxis: {
        type: 'value',
        name: 'Avg silhouette',
        min: 0,
        max: maxVal * 1.05,
        nameTextStyle: { color: '#94a3b8' },
        axisLabel: { color: '#94a3b8' },
        splitLine: { lineStyle: { color: '#1e293b' } },
      },
      series: [
        {
          type: 'line',
          data: [
            [0, 0],
            [maxVal, maxVal],
          ],
          symbol: 'none',
          lineStyle: { color: '#334155', type: 'dashed', width: 1 },
          silent: true,
          z: 0,
        },
        {
          name: 'k values',
          type: 'scatter',
          data: points,
          z: 2,
        },
      ],
    }
  }, [kSummary, kValues, optimalK, currentK])

  return (
    <div>
      <p className="mb-2 text-xs text-slate-500">
        Stability vs silhouette per k — optimal k should sit closest to the top-right.
      </p>
      <ReactECharts option={option} style={{ height: 300 }} opts={{ renderer: 'canvas' }} notMerge />
    </div>
  )
}
