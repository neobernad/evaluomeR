import { AnimatePresence, motion } from 'framer-motion'
import { SlidersHorizontal } from 'lucide-react'
import type { KSummaryEntry } from '@/types/demo'

interface KSliderProps {
  k: number
  kMin: number
  kMax: number
  onChange: (k: number) => void
  optimalK: number
  kSummary: Record<string, KSummaryEntry>
  compact?: boolean
}

function scoreBarColor(ratio: number, isOptimal: boolean): string {
  if (isOptimal) return '#22c55e'
  if (ratio >= 0.9) return '#3b82f6'
  if (ratio >= 0.75) return '#60a5fa'
  if (ratio >= 0.6) return '#f59e0b'
  return '#ef4444'
}

export function KSlider({ k, kMin, kMax, onChange, optimalK, kSummary, compact = false }: KSliderProps) {
  const kRange = Array.from({ length: kMax - kMin + 1 }, (_, i) => kMin + i)
  const maxComposite = Math.max(...Object.values(kSummary).map((e) => e.composite))
  const smGridCols: Record<number, string> = {
    3: 'sm:grid-cols-3',
    4: 'sm:grid-cols-4',
    5: 'sm:grid-cols-5',
    6: 'sm:grid-cols-6',
  }
  const selectedEntry = kSummary[String(k)]
  const selectedRatio = selectedEntry ? selectedEntry.composite / maxComposite : 0
  const isSelectedOptimal = k === optimalK

  return (
    <motion.div
      layout
      transition={{ duration: 0.35, ease: [0.4, 0, 0.2, 1] }}
      className={
        compact
          ? 'border-b border-slate-800 bg-slate-950/90 px-6 pb-4 pt-5 backdrop-blur-sm'
          : 'rounded-xl border border-slate-800 bg-slate-900/50 p-6'
      }
    >
      <AnimatePresence mode="wait" initial={false}>
        {compact ? (
          <motion.div
            key="compact"
            initial={{ opacity: 0, y: -6 }}
            animate={{ opacity: 1, y: 0 }}
            exit={{ opacity: 0, y: -6 }}
            transition={{ duration: 0.2 }}
            className="flex items-center gap-3"
          >
            <div className="flex shrink-0 items-center gap-2">
              <SlidersHorizontal className="h-4 w-4 text-slate-500" aria-hidden />
              <span className="text-sm font-medium text-slate-300">Clusters (k)</span>
            </div>
            <div className="mx-1 h-5 w-px shrink-0 bg-slate-700" aria-hidden />
            <div className="flex flex-1 flex-wrap items-center justify-center gap-2.5">
              {kRange.map((kv) => {
                const isSelected = kv === k
                const isOptimal = kv === optimalK
                return (
                  <motion.button
                    key={kv}
                    type="button"
                    onClick={() => onChange(kv)}
                    whileHover={!isSelected ? { scale: 1.08 } : undefined}
                    whileTap={{ scale: 0.95 }}
                    className={[
                      'relative flex h-9 w-10 items-center justify-center rounded-lg border text-sm font-bold transition-all',
                      isSelected
                        ? isOptimal
                          ? 'border-emerald-500/60 bg-emerald-950/40 text-emerald-300'
                          : 'border-blue-500/50 bg-blue-950/30 text-blue-300'
                        : 'border-slate-700/60 bg-slate-900/40 text-slate-400 opacity-60 hover:opacity-90',
                    ].join(' ')}
                  >
                    {isOptimal && (
                      <span className="absolute -top-1.5 -right-1.5 flex h-4 w-4 items-center justify-center rounded-full bg-gradient-to-br from-emerald-500 to-teal-400 text-[8px] font-semibold text-white shadow-[0_0_6px_rgba(16,185,129,0.4)]">
                        ✦
                      </span>
                    )}
                    {kv}
                  </motion.button>
                )
              })}
            </div>
            {selectedEntry && (
              <span className="flex shrink-0 items-center gap-1.5">
                <span
                  className="h-2 w-2 rounded-full"
                  style={{ backgroundColor: scoreBarColor(selectedRatio, isSelectedOptimal) }}
                />
                <span className="font-mono text-xs text-slate-400">
                  {selectedEntry.composite.toFixed(3)}
                </span>
              </span>
            )}
          </motion.div>
        ) : (
          <motion.div
            key="full"
            initial={{ opacity: 0, y: 6 }}
            animate={{ opacity: 1, y: 0 }}
            exit={{ opacity: 0, y: 6 }}
            transition={{ duration: 0.2 }}
          >
            <div className="mb-5">
              <p className="text-sm font-medium text-slate-300">Number of clusters (k)</p>
              <p className="text-xs text-slate-500">
                Select k to compare stability × quality scores — computed by evaluomeR
              </p>
            </div>

            <div
              className={[
                'mt-4 grid grid-cols-3 gap-3 pt-4',
                smGridCols[kRange.length] ?? 'sm:grid-cols-6',
              ].join(' ')}
            >
              {kRange.map((kv) => {
                const entry = kSummary[String(kv)]
                const isSelected = kv === k
                const isOptimal = kv === optimalK
                const ratio = entry ? entry.composite / maxComposite : 0
                const barColor = scoreBarColor(ratio, isOptimal)

                return (
                  <motion.button
                    key={kv}
                    type="button"
                    onClick={() => onChange(kv)}
                    whileHover={!isSelected ? { scale: 1.04, y: -2 } : undefined}
                    whileTap={{ scale: 0.97 }}
                    animate={{
                      y: isSelected ? -4 : 0,
                      boxShadow: isSelected
                        ? isOptimal
                          ? '0 8px 24px rgba(34,197,94,0.25)'
                          : '0 8px 24px rgba(59,130,246,0.2)'
                        : '0 0 0 rgba(0,0,0,0)',
                    }}
                    transition={{ type: 'spring', stiffness: 400, damping: 25 }}
                    className={[
                      'relative flex flex-col items-center rounded-xl border px-2 py-4 text-center transition-all',
                      isSelected
                        ? isOptimal
                          ? 'border-emerald-500/60 bg-emerald-950/30 opacity-100'
                          : 'border-blue-500/50 bg-blue-950/20 opacity-100'
                        : 'border-slate-700/60 bg-slate-900/40 opacity-60 hover:opacity-80',
                    ].join(' ')}
                  >
                    {isOptimal && (
                      <span className="absolute -top-3.5 left-1/2 -translate-x-1/2 whitespace-nowrap rounded-full bg-gradient-to-r from-emerald-500 to-teal-400 px-2 py-0.5 text-[9px] font-semibold tracking-wide text-white shadow-[0_0_8px_rgba(16,185,129,0.45)]">
                        ✦ optimal k
                      </span>
                    )}
                    <span
                      className={[
                        'font-mono text-xl font-bold leading-none sm:text-3xl',
                        isSelected
                          ? isOptimal
                            ? 'text-emerald-300'
                            : 'text-blue-300'
                          : 'text-slate-400',
                      ].join(' ')}
                    >
                      {kv}
                    </span>

                    <div className="mt-3 h-1.5 w-full overflow-hidden rounded-full bg-slate-800">
                      <motion.div
                        className="h-full rounded-full"
                        style={{ backgroundColor: barColor }}
                        initial={{ width: 0 }}
                        animate={{ width: `${ratio * 100}%` }}
                        transition={{ duration: 0.5, delay: 0.1 }}
                      />
                    </div>

                    {entry && (
                      <span className="mt-2 font-mono text-[10px] text-slate-500">
                        {entry.composite.toFixed(3)}
                      </span>
                    )}
                  </motion.button>
                )
              })}
            </div>

            <p className="mt-4 text-center text-xs text-slate-600">
              Composite = stability × silhouette (all-metrics index)
            </p>
          </motion.div>
        )}
      </AnimatePresence>
    </motion.div>
  )
}
