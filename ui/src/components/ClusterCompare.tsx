import { useMemo } from 'react'
import { motion } from 'framer-motion'
import ReactECharts from 'echarts-for-react'
import type { DemoData, Sample } from '@/types/demo'
import { buildPCAClusterOption } from '@/lib/clusterChart'

interface ClusterCompareProps {
  samples: Sample[]
  clusters: DemoData['clusters']
  userK: number
  optimalK: number
  pca: DemoData['pca']
}

export function ClusterCompare({
  samples,
  clusters,
  userK,
  optimalK,
  pca,
}: ClusterCompareProps) {
  const isOptimal = userK === optimalK

  const userOption = useMemo(
    () =>
      buildPCAClusterOption({
        title: `Your k = ${userK}`,
        titleColor: isOptimal ? '#22c55e' : '#f59e0b',
        samples,
        clusters,
        k: userK,
        pca,
        compact: true,
      }),
    [samples, clusters, userK, pca, isOptimal],
  )

  const optOption = useMemo(
    () =>
      buildPCAClusterOption({
        title: `Recommended k = ${optimalK}`,
        titleColor: '#22c55e',
        samples,
        clusters,
        k: optimalK,
        pca,
        compact: true,
      }),
    [samples, clusters, optimalK, pca],
  )

  return (
    <div className="space-y-3">
      <div className="flex items-center gap-2 text-sm">
        <span className="h-2 w-2 rounded-full bg-blue-500" />
        <span className="text-slate-400">
          PCA projection · color = cancer tissue · ellipses = cluster groups
        </span>
      </div>
      <div className="grid gap-4 lg:grid-cols-2">
        <motion.div
          key={userK}
          className={[
            'rounded-xl border p-1',
            isOptimal ? 'border-emerald-500/40 bg-emerald-950/10' : 'border-slate-700/60',
          ].join(' ')}
          initial={{ opacity: 0.5 }}
          animate={{ opacity: 1 }}
          transition={{ duration: 0.3 }}
        >
          <ReactECharts option={userOption} style={{ height: 340 }} opts={{ renderer: 'canvas' }} notMerge />
        </motion.div>

        <div className="rounded-xl border border-emerald-500/30 bg-emerald-950/10 p-1">
          <ReactECharts option={optOption} style={{ height: 340 }} opts={{ renderer: 'canvas' }} notMerge />
        </div>
      </div>
      {!isOptimal && (
        <p className="text-center text-xs text-amber-400/80">
          At k = {userK}, cancer types mix across groups. At k = {optimalK} they separate more cleanly.
        </p>
      )}
    </div>
  )
}
