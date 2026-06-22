import { useMemo } from 'react'
import ReactECharts from 'echarts-for-react'
import type { KSummaryEntry } from '@/types/demo'
import { buildOptimalAnimatedLine, buildGraphicOptimalLabel } from '@/lib/kChartMarks'

interface KScoreChartProps {
  kSummary: Record<string, KSummaryEntry>
  optimalK: number
  currentK: number
  dataset?: string
}

export function KScoreChart({ kSummary, optimalK, currentK, dataset }: KScoreChartProps) {
  const kValues = Object.keys(kSummary).sort((a, b) => Number(a) - Number(b))
  const optLabel = `k = ${optimalK}`
  const currentLabel = `k = ${currentK}`
  const isAtOptimal = currentK === optimalK

  const option = useMemo(() => {
    const categories = kValues.map((kv) => `k = ${kv}`)
    const stabilityData = kValues.map((k) => ({
      value: kSummary[k].avgStability,
      itemStyle: {
        color:
          Number(k) === optimalK
            ? '#3b82f6'
            : Number(k) === currentK
              ? '#60a5fa'
              : 'rgba(59,130,246,0.35)',
      },
    }))
    const silhouetteData = kValues.map((k) => ({
      value: kSummary[k].avgSilhouette,
      itemStyle: {
        color:
          Number(k) === optimalK
            ? '#22c55e'
            : Number(k) === currentK
              ? '#4ade80'
              : 'rgba(34,197,94,0.35)',
      },
    }))
    const compositeData = kValues.map((k) => ({
      value: kSummary[k].composite,
      symbolSize: Number(k) === currentK ? 14 : 8,
      itemStyle: { color: '#f59e0b' },
    }))

    const markAreaData = (isAtOptimal
      ? [
          [
            { xAxis: optLabel, itemStyle: { color: 'rgba(34,197,94,0.15)' } },
            { xAxis: optLabel },
          ],
        ]
      : [
          [
            { xAxis: currentLabel, itemStyle: { color: 'rgba(59,130,246,0.15)' } },
            { xAxis: currentLabel },
          ],
          [
            { xAxis: optLabel, itemStyle: { color: 'rgba(34,197,94,0.08)' } },
            { xAxis: optLabel },
          ],
        ]) as Array<[{ xAxis: string; itemStyle?: { color: string } }, { xAxis: string }]>

    return {
      backgroundColor: 'transparent',
      tooltip: {
        trigger: 'axis',
        axisPointer: { type: 'shadow' },
        formatter: (params: { name: string; seriesName: string; value: number }[]) => {
          const k = params[0]?.name
          const lines = params.map(
            (p) => `${p.seriesName}: <b>${p.value.toFixed(3)}</b>`,
          )
          return [`<b>${k}</b>`, ...lines].join('<br/>')
        },
      },
      legend: {
        data: ['Avg stability', 'Avg silhouette', 'Composite score'],
        textStyle: { color: '#94a3b8' },
        top: 4,
      },
      grid: { left: '8%', right: '6%', bottom: '14%', top: '18%' },
      xAxis: {
        type: 'category',
        data: categories,
        axisLabel: { color: '#94a3b8', fontWeight: 'bold' },
        axisLine: { lineStyle: { color: '#334155' } },
      },
      yAxis: [
        {
          type: 'value',
          min: 0,
          max: 1,
          name: 'Score',
          nameTextStyle: { color: '#94a3b8' },
          axisLabel: { color: '#94a3b8' },
          splitLine: { lineStyle: { color: '#1e293b' } },
        },
      ],
      series: [
        {
          name: 'Avg stability',
          type: 'bar',
          barMaxWidth: 32,
          data: stabilityData,
          markArea: {
            silent: true,
            data: markAreaData,
          },
        },
        {
          name: 'Avg silhouette',
          type: 'bar',
          barMaxWidth: 32,
          data: silhouetteData,
        },
        {
          name: 'Composite score',
          type: 'line',
          smooth: true,
          symbol: 'circle',
          lineStyle: { color: '#f59e0b', width: 3 },
          itemStyle: { color: '#f59e0b' },
          data: compositeData,
          ...buildOptimalAnimatedLine(optLabel, 0, 1),
        },
      ],
      ...buildGraphicOptimalLabel(
        optLabel,
        categories,
        { left: 8, right: 6, top: 18 },
        `✦ best  k = ${optimalK}`,
      ),
    }
  }, [kSummary, kValues, optimalK, currentK, optLabel, currentLabel, isAtOptimal])

  return (
    <div>
      <ReactECharts option={option} style={{ height: 320 }} opts={{ renderer: 'canvas' }} notMerge />
      {dataset === 'nci60_k8' && (
        <p className="mt-1 text-center text-xs text-slate-500">
          k = 2 is excluded — too few groups for {'>'}9 cancer tissue types
        </p>
      )}
    </div>
  )
}
