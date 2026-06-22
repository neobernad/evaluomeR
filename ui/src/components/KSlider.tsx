import { motion } from 'framer-motion'
import type { KSummaryEntry } from '@/types/demo'

interface KSliderProps {
  k: number
  kMin: number
  kMax: number
  onChange: (k: number) => void
  optimalK: number
  kSummary: Record<string, KSummaryEntry>
}

function scoreBarColor(ratio: number, isOptimal: boolean): string {
  if (isOptimal) return '#22c55e'
  if (ratio >= 0.9) return '#3b82f6'
  if (ratio >= 0.75) return '#60a5fa'
  if (ratio >= 0.6) return '#f59e0b'
  return '#ef4444'
}

export function KSlider({ k, kMin, kMax, onChange, optimalK, kSummary }: KSliderProps) {
  const kRange = Array.from({ length: kMax - kMin + 1 }, (_, i) => kMin + i)
  const maxComposite = Math.max(...Object.values(kSummary).map((e) => e.composite))

  return (
    <div className="rounded-xl border border-slate-800 bg-slate-900/50 p-6">
      <div className="mb-5">
        <p className="text-sm font-medium text-slate-300">Number of clusters (k)</p>
        <p className="text-xs text-slate-500">
          Select k to compare stability × quality scores — computed by evaluomeR
        </p>
      </div>

      <div
        className="grid gap-3"
        style={{ gridTemplateColumns: `repeat(${kRange.length}, minmax(0, 1fr))` }}
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
                <span className="absolute -top-2 rounded-full bg-emerald-600 px-1.5 py-0.5 text-[9px] font-bold text-white">
                  ★ optimal
                </span>
              )}
              <span
                className={[
                  'font-mono text-3xl font-bold leading-none',
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
    </div>
  )
}
