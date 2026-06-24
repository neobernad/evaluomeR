import { useEffect, useMemo, useState } from 'react'
import ReactECharts from 'echarts-for-react'
import { motion } from 'framer-motion'
import { ClusterRadar } from '@/components/ClusterRadar'
import type { DemoData, Sample } from '@/types/demo'
import { buildPCACluster3DOption, buildPCAClusterOption } from '@/lib/clusterChart'

interface ClusterScatterProps {
  samples: Sample[]
  clusters: DemoData['clusters']
  k: number
  pca: DemoData['pca']
  optimalK: number
  chartMetrics: string[]
}

type ClusterView = '2d' | '3d' | 'profiles'

export function ClusterScatter({
  samples,
  clusters,
  k,
  pca,
  optimalK,
  chartMetrics,
}: ClusterScatterProps) {
  const isOptimal = k === optimalK
  const has3D = pca.coords.length > 0 && pca.coords[0].pc3 !== undefined
  const [view, setView] = useState<ClusterView>('2d')
  const [glReady, setGlReady] = useState(false)

  const activeView: ClusterView = view === '3d' && !has3D ? '2d' : view

  useEffect(() => {
    if (activeView === '3d') {
      import('echarts-gl').then(() => setGlReady(true))
    }
  }, [activeView])

  const clusterAssignments = clusters[String(k)] ?? []

  const option = useMemo(() => {
    if (activeView === 'profiles') return null
    if (activeView === '3d') {
      if (!glReady) return null
      return buildPCACluster3DOption({ samples, clusters, k, pca })
    }
    return buildPCAClusterOption({ samples, clusters, k, pca })
  }, [activeView, samples, clusters, k, pca, glReady])

  const pcCount = activeView === '3d' ? 3 : 2
  const totalVar = pca.varExplained
    .slice(0, pcCount)
    .reduce((a, b) => a + b, 0)

  const caption =
    activeView === 'profiles'
      ? 'Normalized mean metric profile per cluster'
      : activeView === '3d'
        ? 'PCA projection (PC1/PC2/PC3) · spheres wrap k-means cluster assignments'
        : 'PCA projection of all metrics · ellipses wrap k-means cluster assignments'

  const viewModes: ClusterView[] = has3D ? ['2d', '3d', 'profiles'] : ['2d', 'profiles']

  return (
    <div>
      <div className="mb-3 flex flex-wrap items-center justify-between gap-3">
        <div className="space-y-1">
          <p className="text-xs text-slate-500">{caption}</p>
          {activeView !== 'profiles' && (
            <span className="inline-flex rounded-full border border-slate-800 bg-slate-900/60 px-2 py-0.5 text-[10px] text-slate-400">
              PC1{pcCount > 2 ? '+PC2+PC3' : '+PC2'} explain {(totalVar * 100).toFixed(1)}% of variance
            </span>
          )}
        </div>
        <div className="inline-flex rounded-lg border border-slate-800 bg-slate-900/60 p-1">
          {viewModes.map((mode) => {
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
      </div>
      {!isOptimal && activeView !== 'profiles' && (
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
        {activeView === 'profiles' ? (
          <ClusterRadar
            samples={samples}
            clusters={clusterAssignments}
            chartMetrics={chartMetrics}
          />
        ) : activeView === '3d' && !glReady ? (
          <div className="flex h-[420px] items-center justify-center text-sm text-slate-500">
            Loading 3D view…
          </div>
        ) : option ? (
          <ReactECharts option={option} style={{ height: 420 }} opts={{ renderer: 'canvas' }} notMerge />
        ) : null}
      </motion.div>
    </div>
  )
}
