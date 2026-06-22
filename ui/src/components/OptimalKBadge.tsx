import { motion } from 'framer-motion'
import { Badge } from '@/components/ui/badge'

interface OptimalKBadgeProps {
  optimalK: number
  currentK: number
}

export function OptimalKBadge({ optimalK, currentK }: OptimalKBadgeProps) {
  const isOptimal = currentK === optimalK

  return (
    <motion.div
      className={[
        'flex flex-wrap items-start gap-2 gap-x-3 rounded-lg border px-4 py-3 transition-colors duration-300',
        isOptimal
          ? 'border-emerald-500/40 bg-emerald-500/10'
          : 'border-amber-500/20 bg-amber-500/5',
      ].join(' ')}
      initial={{ opacity: 0, scale: 0.95 }}
      animate={{ opacity: 1, scale: 1 }}
      transition={{ duration: 0.4 }}
    >
      <Badge variant={isOptimal ? 'success' : 'outline'} className={isOptimal ? '' : 'text-amber-400 border-amber-500/40'}>
        Recommended k = {optimalK}
      </Badge>
      <p className={['text-sm transition-colors', isOptimal ? 'text-emerald-300' : 'text-slate-400'].join(' ')}>
        {isOptimal
          ? 'You are viewing the optimal k — stability and quality both peak here.'
          : `Current k = ${currentK} deviates from the evaluomeR recommendation.`}
      </p>
    </motion.div>
  )
}
