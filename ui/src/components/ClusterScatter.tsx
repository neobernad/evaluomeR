import { useMemo, useState } from 'react'
import ReactECharts from 'echarts-for-react'
import { motion } from 'framer-motion'
import type { DemoData, Sample } from '@/types/demo'
import { buildPCACluster3DOption, buildPCAClusterOption } from '@/lib/clusterChart'

interface ClusterScatterProps {
  samples: Sample[]
  clusters: DemoData['clusters']
  k: number
  pca: DemoData['pca']
  optimalK: number
}

type ClusterView = '2d' | '3d'

export function ClusterScatter({ samples, clusters, k, pca, optimalK }: ClusterScatterProps) {
  const isOptimal = k === optimalK
  const has3D = pca.coords.length > 0 && pca.coords[0].pc3 !== undefined
  const [view, setView] = useState<ClusterView>('2d')

  const activeView: ClusterView = has3D ? view : '2d'

  const option = useMemo(() => {
    if (activeView === '3d') {
      return buildPCACluster3DOption({ samples, clusters, k, pca })
    }
    return buildPCAClusterOption({ samples, clusters, k, pca })
  }, [activeView, samples, clusters, k, pca])

  const caption =
    activeView === '3d'
      ? 'PCA projection (PC1/PC2/PC3) · spheres wrap k-means cluster assignments'
      : 'PCA projection of all metrics · ellipses wrap k-means cluster assignments'

  return (
    <div>
      <div className="mb-3 flex flex-wrap items-center justify-between gap-3">
        <p className="text-xs text-slate-500">{caption}</p>
        {has3D && (
          <div className="inline-flex rounded-lg border border-slate-800 bg-slate-900/60 p-1">
            {(['2d', '3d'] as const).map((mode) => {
              const isActive = activeView === mode
              return (
                <button
                  key={mode}
                  type="button"
                  onClick={() => setView(mode)}
                  className={[
                    'rounded-md px-3 py-1 text-xs font-semibold uppercase tracking-wide transition-all',
                    isActive
                      ? 'bg-slate-800 text-white shadow-sm'
                      : 'text-slate-500 hover:text-slate-300',
                  ].join(' ')}
                >
                  {mode}
                </button>
              )
            })}
          </div>
        )}
      </div>
      {!isOptimal && (
        <p className="mb-2 rounded-lg border border-amber-500/20 bg-amber-500/5 px-3 py-2 text-xs text-amber-400">
          k = {k} is not the recommended k. Cluster boundaries may not align well with cancer tissue types.
        </p>
      )}
      <motion.div
        key={`${k}-${activeView}`}
        initial={{ opacity: 0.4 }}
        animate={{ opacity: 1 }}
        transition={{ duration: 0.35 }}
      >
        <ReactECharts option={option} style={{ height: 420 }} opts={{ renderer: 'canvas' }} notMerge />
      </motion.div>
    </div>
  )
}
