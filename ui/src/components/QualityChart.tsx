import { useMemo } from 'react'
import ReactECharts from 'echarts-for-react'
import { motion } from 'framer-motion'
import { OptimalKBadge } from '@/components/OptimalKBadge'
import type { DemoData } from '@/types/demo'
import { formatMetricLabel } from '@/lib/metricLabels'

const LINE_COLORS = [
  '#22c55e', '#38bdf8', '#a78bfa', '#f59e0b',
  '#f43f5e', '#06b6d4', '#fb923c', '#84cc16',
]

function buildKMarks(selectedK: number, optimalK: number) {
  const marks: Record<string, unknown> = {}

  marks.markArea = {
    silent: true,
    data: [
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

  if (selectedK !== optimalK) {
    marks.markLine = {
      silent: true,
      symbol: 'none',
      lineStyle: { color: '#cbd5e1', type: 'dashed', width: 1.5 },
      data: [{ xAxis: String(selectedK) }],
      label: {
        formatter: `k = ${selectedK}`,
        color: '#cbd5e1',
        fontSize: 10,
      },
    }
  }

  return marks
}

interface QualityChartProps {
  quality: DemoData['quality']
  k: number
  optimalK: number
  chartMetrics: string[]
}

function buildLineOption(
  title: string,
  kValues: string[],
  chartMetrics: string[],
  dataByK: DemoData['quality']['silhouette'],
  yMin: number,
  yMax: number,
  selectedK: number,
  optimalK: number,
) {
  const kMarks = buildKMarks(selectedK, optimalK)

  return {
    backgroundColor: 'transparent',
    title: {
      text: title,
      left: 'center',
      textStyle: { color: '#e2e8f0', fontSize: 14, fontWeight: 600 },
    },
    tooltip: { trigger: 'axis' },
    legend: { data: chartMetrics.map(formatMetricLabel), textStyle: { color: '#94a3b8' }, top: 28, type: 'scroll' },
    grid: { left: '10%', right: '6%', bottom: '12%', top: '24%' },
    xAxis: {
      type: 'category',
      data: kValues,
      axisLabel: { color: '#94a3b8' },
    },
    yAxis: {
      type: 'value',
      min: yMin,
      max: yMax,
      axisLabel: { color: '#94a3b8' },
      splitLine: { lineStyle: { color: '#1e293b' } },
    },
    series: chartMetrics.map((metric, idx) => ({
      name: formatMetricLabel(metric),
      type: 'line',
      smooth: true,
      symbol: 'circle',
      symbolSize: (_: number, params: { dataIndex: number }) =>
        kValues[params.dataIndex] === String(selectedK) ? 14 : 8,
      data: kValues.map((kv) => dataByK[kv]?.[metric] ?? null),
      lineStyle: { width: 2.5 },
      itemStyle: { color: LINE_COLORS[idx % LINE_COLORS.length] },
      ...(idx === 0 ? kMarks : {}),
    })),
  }
}

export function QualityChart({ quality, k, optimalK, chartMetrics }: QualityChartProps) {
  const kValues = Object.keys(quality.silhouette).sort((a, b) => Number(a) - Number(b))

  const chMax = useMemo(() => {
    const vals = kValues.flatMap((kv) =>
      chartMetrics.map((m) => quality.ch[kv]?.[m] ?? 0),
    )
    return Math.ceil(Math.max(...vals) * 1.1)
  }, [quality.ch, kValues, chartMetrics])

  const silOption = useMemo(
    () => buildLineOption('Silhouette width', kValues, chartMetrics, quality.silhouette, -0.2, 1, k, optimalK),
    [quality.silhouette, kValues, chartMetrics, k, optimalK],
  )

  const chOption = useMemo(
    () => buildLineOption('Calinski–Harabasz index', kValues, chartMetrics, quality.ch, 0, chMax, k, optimalK),
    [quality.ch, kValues, chartMetrics, k, optimalK, chMax],
  )

  return (
    <motion.div
      key={k}
      className="space-y-6"
      initial={{ opacity: 0.4 }}
      animate={{ opacity: 1 }}
      transition={{ duration: 0.35 }}
    >
      <OptimalKBadge optimalK={optimalK} currentK={k} />
      <div className="grid gap-6 lg:grid-cols-2">
        <ReactECharts option={silOption} style={{ height: 340 }} opts={{ renderer: 'canvas' }} />
        <ReactECharts option={chOption} style={{ height: 340 }} opts={{ renderer: 'canvas' }} />
      </div>
    </motion.div>
  )
}
