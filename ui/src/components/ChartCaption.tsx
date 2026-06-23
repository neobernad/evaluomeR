import type { LucideIcon } from 'lucide-react'
import type { ReactNode } from 'react'
import { Badge } from '@/components/ui/badge'
import { cn } from '@/lib/utils'

export interface Highlight {
  label: string
  color: 'blue' | 'emerald' | 'amber' | 'red' | 'cyan' | 'slate'
}

interface ChartCaptionProps {
  icon: LucideIcon
  text: ReactNode
  highlights?: Highlight[]
}

const HIGHLIGHT_STYLES: Record<Exclude<Highlight['color'], 'slate'>, string> = {
  blue: 'bg-blue-500/15 text-blue-300 ring-blue-500/30',
  emerald: 'bg-emerald-500/15 text-emerald-300 ring-emerald-500/30',
  amber: 'bg-amber-500/15 text-amber-300 ring-amber-500/30',
  red: 'bg-red-500/15 text-red-300 ring-red-500/30',
  cyan: 'bg-cyan-500/15 text-cyan-300 ring-cyan-500/30',
}

export function ChartCaption({ icon: Icon, text, highlights }: ChartCaptionProps) {
  return (
    <div className="mt-3 border-t border-slate-800 pt-3">
      <div className="flex flex-col gap-2">
        <div className="flex gap-2.5">
          <Icon className="mt-0.5 h-4 w-4 shrink-0 text-slate-500" aria-hidden />
          <p className="text-xs leading-relaxed text-slate-400">{text}</p>
        </div>
        {highlights && highlights.length > 0 && (
          <div className="flex flex-wrap gap-1.5 pl-6.5">
            {highlights.map((h) => (
              <Badge
                key={h.label}
                variant={h.color === 'slate' ? 'outline' : 'default'}
                className={cn(
                  'px-2 py-0.5 text-[10px] font-medium',
                  h.color !== 'slate' && HIGHLIGHT_STYLES[h.color],
                )}
              >
                {h.label}
              </Badge>
            ))}
          </div>
        )}
      </div>
    </div>
  )
}
