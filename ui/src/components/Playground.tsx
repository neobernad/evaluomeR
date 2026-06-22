import { useState } from 'react'
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from '@/components/ui/card'
import { Tabs, TabsContent, TabsList, TabsTrigger } from '@/components/ui/tabs'
import { DatasetSelector, type DatasetKey, type DatasetOption } from '@/components/DatasetSelector'
import { KSlider } from '@/components/KSlider'
import { OptimalKPanel } from '@/components/OptimalKPanel'
import { KScoreChart } from '@/components/KScoreChart'
import { ClusterCompare } from '@/components/ClusterCompare'
import { OptimalKTable } from '@/components/OptimalKTable'
import { StabilityChart } from '@/components/StabilityChart'
import { QualityChart } from '@/components/QualityChart'
import { ClusterScatter } from '@/components/ClusterScatter'
import type { DemoData } from '@/types/demo'

interface PlaygroundProps {
  datasets: DatasetOption[]
  defaultDataset?: DatasetKey
}

export function Playground({ datasets, defaultDataset = 'nci60' }: PlaygroundProps) {
  const [datasetKey, setDatasetKey] = useState<DatasetKey>(defaultDataset)
  const active = datasets.find((d) => d.key === datasetKey) ?? datasets[0]
  const data: DemoData = active.data
  const [k, setK] = useState<number>(data.optimalK)

  const handleDatasetChange = (key: DatasetKey) => {
    const next = datasets.find((d) => d.key === key)
    if (!next) return
    setDatasetKey(key)
    setK(next.data.optimalK)
  }

  const typeLabel = datasetKey === 'golub' ? 'leukemia classes' : 'cancer tissues'
  const subtitle =
    datasetKey === 'golub'
      ? `Explore precomputed evaluomeR analyses on ${data.meta.nSamples} Golub leukemia samples.`
      : `Explore precomputed evaluomeR analyses on ${data.meta.nSamples} NCI-60 cancer cell lines.`

  return (
    <section id="playground" className="mx-auto max-w-5xl px-6 py-20">
      <div className="mb-8 text-center">
        <h2 className="text-3xl font-bold text-white">Interactive Playground</h2>
        <p className="mt-2 text-slate-400">{subtitle}</p>
      </div>

      <div className="mb-6">
        <DatasetSelector
          datasets={datasets}
          selected={datasetKey}
          onSelect={handleDatasetChange}
        />
      </div>

      <div className="mb-8">
        <KSlider
          k={k}
          kMin={data.meta.kRange[0]}
          kMax={data.meta.kRange[1]}
          onChange={setK}
          optimalK={data.optimalK}
          kSummary={data.kSummary}
        />
      </div>

      <Card>
        <CardHeader>
          <CardTitle>Analysis views</CardTitle>
          <CardDescription>
            Optimal k, stability, quality indices, and cluster assignments for k = {k}.
          </CardDescription>
        </CardHeader>
        <CardContent>
          <Tabs defaultValue="optimal-k">
            <TabsList className="w-full flex-wrap">
              <TabsTrigger value="optimal-k">Optimal k</TabsTrigger>
              <TabsTrigger value="stability">Stability</TabsTrigger>
              <TabsTrigger value="quality">Quality</TabsTrigger>
              <TabsTrigger value="clusters">Clusters</TabsTrigger>
            </TabsList>

            <TabsContent value="optimal-k" className="space-y-6 pt-4">
              <OptimalKPanel
                optimalK={data.optimalK}
                nCancerTypes={data.meta.nCancerTypes}
                kSummary={data.kSummary}
                currentK={k}
                typeLabel={typeLabel}
                onJumpToOptimal={() => setK(data.optimalK)}
              />
              <KScoreChart
                kSummary={data.kSummary}
                optimalK={data.optimalK}
                dataset={data.meta.dataset}
              />
              <ClusterCompare
                samples={data.samples}
                clusters={data.clusters}
                userK={k}
                optimalK={data.optimalK}
                pca={data.pca}
              />
              <OptimalKTable
                optimalKDetail={data.optimalKDetail}
                optimalK={data.optimalK}
              />
            </TabsContent>

            <TabsContent value="stability">
              <StabilityChart
                stability={data.stability}
                k={k}
                optimalK={data.optimalK}
                chartMetrics={data.meta.chartMetrics}
              />
            </TabsContent>

            <TabsContent value="quality">
              <QualityChart
                quality={data.quality}
                k={k}
                optimalK={data.optimalK}
                chartMetrics={data.meta.chartMetrics}
              />
            </TabsContent>

            <TabsContent value="clusters">
              <ClusterScatter
                samples={data.samples}
                clusters={data.clusters}
                k={k}
                pca={data.pca}
                optimalK={data.optimalK}
              />
            </TabsContent>
          </Tabs>
        </CardContent>
      </Card>
    </section>
  )
}
