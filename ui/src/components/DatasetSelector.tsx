import { motion } from 'framer-motion'
import type { DemoData } from '@/types/demo'

export type DatasetKey = 'nci60' | 'golub'

export interface DatasetOption {
  key: DatasetKey
  label: string
  description: string
  data: DemoData
}

interface DatasetSelectorProps {
  datasets: DatasetOption[]
  selected: DatasetKey
  onSelect: (key: DatasetKey) => void
}

function Chip({
  label,
  value,
  accent,
}: {
  label: string
  value: string
  accent?: 'blue' | 'violet'
}) {
  const bg =
    accent === 'blue'
      ? 'bg-blue-950/60 ring-blue-800/40'
      : accent === 'violet'
        ? 'bg-violet-950/60 ring-violet-800/40'
        : 'bg-slate-800/60 ring-slate-700/40'
  const val =
    accent === 'blue'
      ? 'text-blue-300'
      : accent === 'violet'
        ? 'text-violet-300'
        : 'text-slate-300'

  return (
    <span
      className={`inline-flex items-center gap-1 rounded-full px-2 py-0.5 text-[10px] font-mono ring-1 ${bg}`}
    >
      <span className="text-slate-500">{label}</span>
      <span className={val}>{value}</span>
    </span>
  )
}

export function DatasetSelector({ datasets, selected, onSelect }: DatasetSelectorProps) {
  return (
    <div className="grid gap-3 sm:grid-cols-2">
      {datasets.map((ds) => {
        const isSelected = ds.key === selected
        const [kMin, kMax] = ds.data.meta.kRange

        return (
          <motion.button
            key={ds.key}
            type="button"
            onClick={() => onSelect(ds.key)}
            whileHover={!isSelected ? { scale: 1.02 } : undefined}
            whileTap={{ scale: 0.98 }}
            className={[
              'rounded-xl border p-4 text-left transition-all',
              isSelected
                ? 'border-blue-500/60 bg-blue-950/30 shadow-lg shadow-blue-900/20'
                : 'border-slate-800 bg-slate-900/30 opacity-60 hover:opacity-80',
            ].join(' ')}
          >
            <p className={['text-base font-semibold', isSelected ? 'text-white' : 'text-slate-400'].join(' ')}>
              {ds.label}
            </p>
            <div className="mt-2 flex flex-wrap gap-1.5">
              <Chip label="samples" value={String(ds.data.meta.nSamples)} />
              <Chip
                label={ds.key === 'golub' ? 'classes' : 'types'}
                value={String(ds.data.meta.nCancerTypes)}
              />
              <Chip label="k" value={`${kMin}–${kMax}`} accent="blue" />
              <Chip label="bs" value={String(ds.data.meta.bs)} accent="violet" />
            </div>
            <p className="mt-2 text-[11px] text-slate-500">{ds.description}</p>
          </motion.button>
        )
      })}
    </div>
  )
}
