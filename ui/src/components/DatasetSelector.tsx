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

export function DatasetSelector({ datasets, selected, onSelect }: DatasetSelectorProps) {
  return (
    <div className="grid gap-3 sm:grid-cols-2">
      {datasets.map((ds) => {
        const isSelected = ds.key === selected
        const [kMin, kMax] = ds.data.meta.kRange
        const typeLabel =
          ds.key === 'golub'
            ? `${ds.data.meta.nCancerTypes} leukemia classes`
            : `${ds.data.meta.nCancerTypes} tissue types`

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
            <p className="mt-1 text-xs text-slate-500">
              {ds.data.meta.nSamples} samples · {typeLabel}
            </p>
            <p className="mt-0.5 text-xs text-slate-600">
              k = {kMin}–{kMax} · bootstrap = {ds.data.meta.bs}
            </p>
            <p className="mt-2 text-[11px] text-slate-600">{ds.description}</p>
          </motion.button>
        )
      })}
    </div>
  )
}
