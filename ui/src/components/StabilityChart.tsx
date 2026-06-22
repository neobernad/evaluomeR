import { useMemo } from 'react'
import ReactECharts from 'echarts-for-react'
import { motion } from 'framer-motion'
import type { DemoData } from '@/types/demo'
import { formatMetricLabel } from '@/lib/metricLabels'

const ZONES = [
  { name: 'Unstable',      min: 0,    max: 0.6,  color: 'rgba(239,68,68,0.10)' },
  { name: 'Doubtful',      min: 0.6,  max: 0.75, color: 'rgba(234,179,8,0.10)' },
  { name: 'Stable',        min: 0.75, max: 0.85, color: 'rgba(34,197,94,0.10)' },
  { name: 'Highly Stable', min: 0.85, max: 1.0,  color: 'rgba(59,130,246,0.10)' },
]

const LINE_COLORS = [
  '#3b82f6', '#a78bfa', '#22c55e', '#f59e0b',
  '#f43f5e', '#06b6d4', '#fb923c', '#84cc16',
]

interface StabilityChartProps {
  stability: DemoData['stability']
  k: number
  optimalK: number
  chartMetrics: string[]
}

export function StabilityChart({ stability, k, optimalK, chartMetrics }: StabilityChartProps) {
  const kValues = Object.keys(stability).sort((a, b) => Number(a) - Number(b))
  const isOptimal = k === optimalK

  const option = useMemo(
    () => ({
      backgroundColor: 'transparent',
      tooltip: { trigger: 'axis' },
      legend: {
        data: chartMetrics.map(formatMetricLabel),
        textStyle: { color: '#94a3b8' },
        top: 0,
        type: 'scroll',
      },
      grid: { left: '8%', right: '4%', bottom: '12%', top: '16%' },
      xAxis: {
        type: 'category',
        data: kValues,
        name: 'k',
        nameTextStyle: { color: '#94a3b8' },
        axisLabel: { color: '#94a3b8' },
      },
      yAxis: {
        type: 'value',
        min: 0,
        max: 1,
        name: 'Jaccard stability',
        nameTextStyle: { color: '#94a3b8' },
        axisLabel: { color: '#94a3b8' },
        splitLine: { lineStyle: { color: '#1e293b' } },
      },
      series: chartMetrics.map((metric, idx) => ({
        name: formatMetricLabel(metric),
        type: 'line',
        smooth: true,
        symbol: 'circle',
        symbolSize: 8,
        data: kValues.map((kv) => stability[kv]?.[metric] ?? null),
        lineStyle: { width: 2.5 },
        itemStyle: { color: LINE_COLORS[idx % LINE_COLORS.length] },
        markArea: idx === 0
          ? {
              silent: true,
              data: [
                ...ZONES.map((z) => [
                  { yAxis: z.min, itemStyle: { color: z.color } },
                  { yAxis: z.max },
                ]),
                [
                  {
                    xAxis: String(optimalK),
                    itemStyle: { color: 'rgba(34,197,94,0.12)' },
                    label: {
                      show: true,
                      formatter: `Optimal k = ${optimalK}`,
                      color: '#22c55e',
                      fontSize: 10,
                      position: 'insideTop',
                    },
                  },
                  { xAxis: String(optimalK) },
                ],
              ],
            }
          : undefined,
        markLine:
          idx === 0 && k !== optimalK
            ? {
                silent: true,
                symbol: 'none',
                lineStyle: { color: '#cbd5e1', type: 'dashed', width: 1.5 },
                data: [{ xAxis: String(k) }],
                label: {
                  formatter: `k = ${k}`,
                  color: '#cbd5e1',
                  fontSize: 10,
                },
              }
            : undefined,
      })),
    }),
    [stability, chartMetrics, kValues, k, optimalK],
  )

  return (
    <div>
      {!isOptimal && (
        <p className="mb-2 rounded-lg border border-amber-500/20 bg-amber-500/5 px-3 py-2 text-xs text-amber-400">
          Viewing k = {k}. Green band marks the recommended k = {optimalK}.
        </p>
      )}
      <motion.div
        key={k}
        initial={{ opacity: 0.4 }}
        animate={{ opacity: 1 }}
        transition={{ duration: 0.35 }}
      >
        <ReactECharts option={option} style={{ height: 380 }} opts={{ renderer: 'canvas' }} />
      </motion.div>
    </div>
  )
}
