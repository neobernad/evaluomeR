import { motion } from 'framer-motion'
import { Check, Zap } from 'lucide-react'
import { Badge } from '@/components/ui/badge'
import type { KSummaryEntry } from '@/types/demo'

interface OptimalKPanelProps {
  optimalK: number
  nCancerTypes: number
  kSummary: Record<string, KSummaryEntry>
  onJumpToOptimal: () => void
  currentK: number
  typeLabel?: string
}

export function OptimalKPanel({
  optimalK,
  nCancerTypes,
  kSummary,
  onJumpToOptimal,
  currentK,
  typeLabel = 'known classes',
}: OptimalKPanelProps) {
  const isAtOptimal = currentK === optimalK
  const optEntry = kSummary[String(optimalK)]

  return (
    <motion.div
      className="rounded-xl border border-slate-700/60 bg-slate-900/50 p-5 space-y-4"
      initial={{ opacity: 0, y: 12 }}
      animate={{ opacity: 1, y: 0 }}
      transition={{ duration: 0.4 }}
    >
      <div className="flex flex-wrap items-start justify-between gap-4">
        <div className="space-y-2">
          <p className="text-sm text-slate-300 max-w-xl">
            Choosing <em>k</em> shapes every cluster — evaluomeR finds the value where metrics are
            both{' '}
            <span className="text-blue-400 font-medium">stable across bootstrap resamples</span> and{' '}
            <span className="text-emerald-400 font-medium">meaningful in cluster quality</span>.
          </p>
          <div className="flex flex-wrap items-center gap-2">
            <Badge variant="success" className="gap-1">
              Recommended k = {optimalK}
            </Badge>
            <Badge variant="outline" className="text-slate-400">
              {nCancerTypes} {typeLabel}
            </Badge>
            {optEntry && (
              <Badge variant="outline" className="text-amber-400 border-amber-500/30">
                composite score {optEntry.composite.toFixed(3)}
              </Badge>
            )}
          </div>
        </div>

        <button
          type="button"
          onClick={isAtOptimal ? undefined : onJumpToOptimal}
          disabled={isAtOptimal}
          className={[
            'inline-flex shrink-0 items-center gap-2 rounded-lg border px-4 py-2 text-sm font-medium transition',
            isAtOptimal
              ? 'cursor-default border-emerald-500/20 bg-emerald-950/10 text-emerald-600'
              : 'cursor-pointer border-emerald-500/40 bg-emerald-600/20 text-emerald-300 hover:bg-emerald-600/30',
          ].join(' ')}
        >
          {isAtOptimal ? <Check className="h-4 w-4" /> : <Zap className="h-4 w-4" />}
          {isAtOptimal ? 'At optimal k' : 'Jump to optimal k'}
        </button>
      </div>
    </motion.div>
  )
}
