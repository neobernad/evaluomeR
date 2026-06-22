import { Card, CardContent, CardDescription, CardHeader, CardTitle } from '@/components/ui/card'
import { Badge } from '@/components/ui/badge'
import { getCancerColor } from '@/lib/cancerColors'
import type { DemoData } from '@/types/demo'

interface DatasetCardProps {
  data: DemoData
}

export function DatasetCard({ data }: DatasetCardProps) {
  // Show top 3 chart metrics in the table to keep it readable
  const displayMetrics = (data.meta.previewMetrics ?? data.meta.chartMetrics).slice(0, 3)
  const isGolub = data.meta.dataset === 'golub'
  const title = isGolub ? 'Golub leukemia dataset' : 'NCI-60 cancer cell lines'
  const typeLabel = isGolub ? 'leukemia classes' : 'tissue types'

  return (
    <Card>
      <CardHeader>
        <div className="flex flex-wrap items-start justify-between gap-3">
          <div>
            <CardTitle>{title}</CardTitle>
            <CardDescription className="mt-1">
              {data.meta.nSamples} samples · {data.meta.nCancerTypes} {typeLabel} ·{' '}
              {data.meta.metrics.length} metrics · bootstrap = {data.meta.bs}
            </CardDescription>
          </div>
          <div className="flex flex-wrap gap-1.5">
            {data.meta.cancerTypes.map((ct) => (
              <Badge
                key={ct}
                variant="outline"
                className="text-[10px]"
                style={{
                  borderColor: getCancerColor(ct) + '60',
                  color: getCancerColor(ct),
                }}
              >
                {ct}
              </Badge>
            ))}
          </div>
        </div>
      </CardHeader>
      <CardContent>
        <div className="overflow-x-auto rounded-lg border border-slate-800">
          <table className="w-full text-left text-sm">
            <thead className="bg-slate-900 text-slate-400">
              <tr>
                <th className="px-4 py-3 font-medium">Sample</th>
                <th className="px-4 py-3 font-medium">{isGolub ? 'Class' : 'Cancer type'}</th>
                {displayMetrics.map((m) => (
                  <th key={m} className="px-4 py-3 font-medium font-mono text-xs">{m}</th>
                ))}
              </tr>
            </thead>
            <tbody>
              {data.samples.slice(0, 20).map((sample, i) => (
                <tr
                  key={sample.id}
                  className={i % 2 === 0 ? 'bg-slate-950/40' : 'bg-slate-900/20'}
                >
                  <td className="px-4 py-2.5 font-mono text-xs text-slate-300">{sample.id}</td>
                  <td className="px-4 py-2.5">
                    <span
                      className="inline-flex items-center gap-1.5 text-xs font-medium"
                      style={{ color: getCancerColor(sample.cancerType) }}
                    >
                      <span
                        className="h-1.5 w-1.5 rounded-full"
                        style={{ backgroundColor: getCancerColor(sample.cancerType) }}
                      />
                      {sample.cancerType}
                    </span>
                  </td>
                  {displayMetrics.map((m) => (
                    <td key={m} className="px-4 py-2.5 text-xs text-slate-400">
                      {(sample.metrics[m] ?? 0).toFixed(2)}
                    </td>
                  ))}
                </tr>
              ))}
            </tbody>
          </table>
        </div>
        {data.samples.length > 20 && (
          <p className="mt-2 text-center text-xs text-slate-600">
            Showing 20 of {data.samples.length} samples
          </p>
        )}
      </CardContent>
    </Card>
  )
}
