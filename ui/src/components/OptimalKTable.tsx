import { useState } from 'react'
import { ChevronDown, ChevronRight } from 'lucide-react'
import type { OptimalKRow } from '@/types/demo'
import { formatMetricLabel } from '@/lib/metricLabels'

interface OptimalKTableProps {
  optimalKDetail: OptimalKRow[]
  optimalK: number
}

export function OptimalKTable({ optimalKDetail, optimalK }: OptimalKTableProps) {
  const [open, setOpen] = useState(false)

  return (
    <div className="rounded-lg border border-slate-800">
      <button
        onClick={() => setOpen((v) => !v)}
        className="flex w-full items-center justify-between px-4 py-3 text-sm font-medium text-slate-300 hover:bg-slate-900/40 transition-colors"
      >
        <span>Per-metric breakdown</span>
        {open ? (
          <ChevronDown className="h-4 w-4 text-slate-500" />
        ) : (
          <ChevronRight className="h-4 w-4 text-slate-500" />
        )}
      </button>

      {open && (
        <div className="overflow-x-auto border-t border-slate-800">
          <table className="w-full text-left text-xs">
            <thead className="bg-slate-900 text-slate-400">
              <tr>
                <th className="px-4 py-2.5 font-medium">Metric</th>
                <th className="px-4 py-2.5 font-medium">Best stability k</th>
                <th className="px-4 py-2.5 font-medium">Best quality k</th>
                <th className="px-4 py-2.5 font-medium">Global optimal k</th>
              </tr>
            </thead>
            <tbody>
              {optimalKDetail.map((row, i) => (
                <tr
                  key={row.metric}
                  className={[
                    i % 2 === 0 ? 'bg-slate-950/40' : 'bg-slate-900/20',
                    row.globalOptimalK === optimalK ? 'text-emerald-400' : 'text-slate-300',
                  ].join(' ')}
                >
                  <td className="px-4 py-2 font-mono">{formatMetricLabel(row.metric)}</td>
                  <td className="px-4 py-2 text-center">{row.stabilityMaxK}</td>
                  <td className="px-4 py-2 text-center">{row.qualityMaxK}</td>
                  <td className="px-4 py-2 text-center font-semibold">{row.globalOptimalK}</td>
                </tr>
              ))}
            </tbody>
          </table>
        </div>
      )}
    </div>
  )
}
