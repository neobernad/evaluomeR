import { useMemo } from 'react'
import ReactECharts from 'echarts-for-react'
import { motion } from 'framer-motion'
import type { DemoData, Sample } from '@/types/demo'
import { buildPCAClusterOption } from '@/lib/clusterChart'

interface ClusterScatterProps {
  samples: Sample[]
  clusters: DemoData['clusters']
  k: number
  pca: DemoData['pca']
  optimalK: number
}

export function ClusterScatter({ samples, clusters, k, pca, optimalK }: ClusterScatterProps) {
  const isOptimal = k === optimalK

  const option = useMemo(
    () => buildPCAClusterOption({ samples, clusters, k, pca }),
    [samples, clusters, k, pca],
  )

  return (
    <div>
      <p className="mb-2 text-xs text-slate-500">
        PCA projection of all metrics · ellipses wrap k-means cluster assignments
      </p>
      {!isOptimal && (
        <p className="mb-2 rounded-lg border border-amber-500/20 bg-amber-500/5 px-3 py-2 text-xs text-amber-400">
          k = {k} is not the recommended k. Cluster boundaries may not align well with cancer tissue types.
        </p>
      )}
      <motion.div
        key={k}
        initial={{ opacity: 0.4 }}
        animate={{ opacity: 1 }}
        transition={{ duration: 0.35 }}
      >
        <ReactECharts option={option} style={{ height: 420 }} opts={{ renderer: 'canvas' }} />
      </motion.div>
    </div>
  )
}
