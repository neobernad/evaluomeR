import { useMemo } from 'react'
import ReactECharts from 'echarts-for-react'
import { motion } from 'framer-motion'
import { Shuffle } from 'lucide-react'
import { ChartCaption } from '@/components/ChartCaption'
import type { DemoData } from '@/types/demo'
import { formatMetricLabel } from '@/lib/metricLabels'
import { buildAnimatedLines, buildGraphicKLabels } from '@/lib/kChartMarks'

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
              data: ZONES.map((z) => [
                {
                  yAxis: z.min,
                  itemStyle: { color: z.color },
                  label: {
                    show: true,
                    formatter: z.name,
                    position: 'insideTopLeft',
                    color: '#94a3b8',
                    fontSize: 9,
                    opacity: 0.8,
                  },
                },
                { yAxis: z.max },
              ]),
            }
          : undefined,
        ...(idx === 0 ? buildAnimatedLines(k, optimalK, 0, 1) : {}),
      })),
      ...buildGraphicKLabels(k, optimalK, kValues, { left: 8, right: 4, top: 16 }),
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
        <ReactECharts option={option} style={{ height: 380 }} opts={{ renderer: 'canvas' }} notMerge />
        <ChartCaption
          icon={Shuffle}
          text="Each line is the mean Jaccard similarity of cluster assignments across 100 bootstrap resamples — how consistently a metric groups samples under k."
          highlights={[
            { label: '≥ 0.85 Highly stable', color: 'blue' },
            { label: '0.75 – 0.85 Stable', color: 'emerald' },
            { label: '0.60 – 0.75 Doubtful', color: 'amber' },
            { label: '< 0.60 Unstable', color: 'red' },
          ]}
        />
      </motion.div>
    </div>
  )
}
