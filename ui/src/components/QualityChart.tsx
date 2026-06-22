import { useMemo } from 'react'
import ReactECharts from 'echarts-for-react'
import { motion } from 'framer-motion'
import { BarChart2, Layers } from 'lucide-react'
import { ChartCaption } from '@/components/ChartCaption'
import { OptimalKBadge } from '@/components/OptimalKBadge'
import type { DemoData } from '@/types/demo'
import { formatMetricLabel } from '@/lib/metricLabels'
import { buildAnimatedLines, buildGraphicKLabels } from '@/lib/kChartMarks'

const LINE_COLORS = [
  '#22c55e', '#38bdf8', '#a78bfa', '#f59e0b',
  '#f43f5e', '#06b6d4', '#fb923c', '#84cc16',
]

interface Band {
  yMin: number
  yMax: number
  color: string
  label: string
}

const SILHOUETTE_BANDS: Band[] = [
  { yMin: 0.7, yMax: 1.0, color: 'rgba(34,197,94,0.09)', label: 'Strong structure' },
  { yMin: 0.5, yMax: 0.7, color: 'rgba(6,182,212,0.09)', label: 'Reasonable' },
  { yMin: 0.25, yMax: 0.5, color: 'rgba(245,158,11,0.09)', label: 'Weak' },
  { yMin: 0.0, yMax: 0.25, color: 'rgba(249,115,22,0.09)', label: 'No structure' },
  { yMin: -0.2, yMax: 0.0, color: 'rgba(239,68,68,0.09)', label: 'Misclassified' },
]

function chBands(chMax: number): Band[] {
  return [
    { yMin: chMax * 0.75, yMax: chMax, color: 'rgba(34,197,94,0.09)', label: 'Excellent' },
    { yMin: chMax * 0.5, yMax: chMax * 0.75, color: 'rgba(6,182,212,0.09)', label: 'Good' },
    { yMin: chMax * 0.25, yMax: chMax * 0.5, color: 'rgba(245,158,11,0.09)', label: 'Fair' },
    { yMin: 0, yMax: chMax * 0.25, color: 'rgba(249,115,22,0.09)', label: 'Poor' },
  ]
}

function buildBandSeries(bands: Band[]) {
  return {
    name: '_bands',
    type: 'line',
    data: [],
    silent: true,
    legendHoverLink: false,
    z: 0,
    markArea: {
      silent: true,
      data: bands.map((b) => [
        {
          yAxis: b.yMin,
          itemStyle: { color: b.color },
          label: {
            show: true,
            formatter: b.label,
            position: 'insideTopLeft',
            color: '#94a3b8',
            fontSize: 9,
            opacity: 0.8,
          },
        },
        { yAxis: b.yMax },
      ]),
    },
  }
}

function buildKMarks(selectedK: number, optimalK: number, yMin: number, yMax: number) {
  return {
    markArea: {
      silent: true,
      data: [
        [
          { xAxis: String(optimalK), itemStyle: { color: 'rgba(34,197,94,0.08)' } },
          { xAxis: String(optimalK) },
        ],
      ],
    },
    ...buildAnimatedLines(selectedK, optimalK, yMin, yMax),
  }
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
  bands?: Band[],
) {
  const kMarks = buildKMarks(selectedK, optimalK, yMin, yMax)

  const lineSeries = chartMetrics.map((metric, idx) => ({
    name: formatMetricLabel(metric),
    type: 'line',
    smooth: true,
    symbol: 'circle',
    symbolSize: (_: number, params: { dataIndex: number }) =>
      kValues[params.dataIndex] === String(selectedK) ? 14 : 8,
    data: kValues.map((kv) => dataByK[kv]?.[metric] ?? null),
    lineStyle: { width: 2.5 },
    itemStyle: { color: LINE_COLORS[idx % LINE_COLORS.length] },
    z: 2,
    ...(idx === 0 ? kMarks : {}),
  }))

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
    series: [
      ...(bands ? [buildBandSeries(bands)] : []),
      ...lineSeries,
    ],
    ...buildGraphicKLabels(selectedK, optimalK, kValues, { left: 10, right: 6, top: 24 }),
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
    () =>
      buildLineOption(
        'Silhouette width',
        kValues,
        chartMetrics,
        quality.silhouette,
        -0.2,
        1,
        k,
        optimalK,
        SILHOUETTE_BANDS,
      ),
    [quality.silhouette, kValues, chartMetrics, k, optimalK],
  )

  const chOption = useMemo(
    () =>
      buildLineOption(
        'Calinski–Harabasz index',
        kValues,
        chartMetrics,
        quality.ch,
        0,
        chMax,
        k,
        optimalK,
        chBands(chMax),
      ),
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
        <div>
          <ReactECharts option={silOption} style={{ height: 340 }} opts={{ renderer: 'canvas' }} notMerge />
          <ChartCaption
            icon={Layers}
            text="Silhouette width measures how much closer a sample is to its own cluster than to the nearest neighbour cluster. Values near +1 indicate well-separated clusters."
            highlights={[
              { label: '> 0.70 Strong', color: 'emerald' },
              { label: '0.50 – 0.70 Reasonable', color: 'cyan' },
              { label: '< 0.25 Weak', color: 'amber' },
            ]}
          />
        </div>
        <div>
          <ReactECharts option={chOption} style={{ height: 340 }} opts={{ renderer: 'canvas' }} notMerge />
          <ChartCaption
            icon={BarChart2}
            text="Calinski–Harabasz is the ratio of between-cluster to within-cluster dispersion — higher means tighter, better-separated groups. No fixed scale; compare across k values."
            highlights={[
              { label: 'Higher = better', color: 'blue' },
              { label: 'No upper bound', color: 'slate' },
            ]}
          />
        </div>
      </div>
    </motion.div>
  )
}
