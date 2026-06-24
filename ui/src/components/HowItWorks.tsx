import { useState } from 'react'
import { ChevronDown, ChevronRight, Layers, Shuffle, Zap } from 'lucide-react'

interface HowItWorksProps {
  onNavigateTab?: (tab: string) => void
}

export function HowItWorks({ onNavigateTab }: HowItWorksProps) {
  const [open, setOpen] = useState(false)

  const steps = [
    {
      icon: Shuffle,
      title: 'Bootstrap stability',
      body: (
        <>
          Resample the dataset 100× and cluster each replicate. Jaccard similarity of assignments
          yields a stability index via <code className="text-blue-300">stabilityRange()</code>.
        </>
      ),
      tab: 'stability',
    },
    {
      icon: Layers,
      title: 'Cluster quality',
      body: (
        <>
          For each k, compute silhouette width and Calinski–Harabasz via{' '}
          <code className="text-emerald-300">qualityRange()</code> — higher means tighter,
          better-separated groups.
        </>
      ),
      tab: 'quality',
    },
    {
      icon: Zap,
      title: 'Optimal k selection',
      body: (
        <>
          Composite score = stability × silhouette.{' '}
          <code className="text-amber-300">getOptimalKValue()</code> picks the k that maximises
          this combined index.
        </>
      ),
      tab: 'optimal-k',
    },
  ]

  return (
    <div className="rounded-lg border border-slate-800">
      <button
        type="button"
        onClick={() => setOpen((v) => !v)}
        className="flex w-full items-center justify-between px-4 py-3 text-sm font-medium text-slate-300 transition-colors hover:bg-slate-900/40"
      >
        <span>How evaluomeR decides optimal k</span>
        {open ? (
          <ChevronDown className="h-4 w-4 text-slate-500" />
        ) : (
          <ChevronRight className="h-4 w-4 text-slate-500" />
        )}
      </button>
      {open && (
        <div className="grid gap-3 border-t border-slate-800 p-4 sm:grid-cols-3">
          {steps.map((step, i) => {
            const Icon = step.icon
            return (
              <div
                key={step.title}
                className="rounded-lg border border-slate-800 bg-slate-900/40 p-4"
              >
                <div className="mb-2 flex items-center gap-2">
                  <span className="flex h-6 w-6 items-center justify-center rounded-full bg-slate-800 text-xs font-bold text-slate-400">
                    {i + 1}
                  </span>
                  <Icon className="h-4 w-4 text-slate-500" />
                  <span className="text-sm font-medium text-slate-200">{step.title}</span>
                </div>
                <p className="text-xs leading-relaxed text-slate-400">{step.body}</p>
                {onNavigateTab && (
                  <button
                    type="button"
                    onClick={() => onNavigateTab(step.tab)}
                    className="mt-2 text-xs text-blue-400 hover:underline"
                  >
                    View {step.tab === 'optimal-k' ? 'optimal k' : step.tab} tab →
                  </button>
                )}
              </div>
            )
          })}
        </div>
      )}
    </div>
  )
}
