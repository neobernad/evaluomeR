import { cn } from '@/lib/utils'

export function Badge({
  className,
  variant = 'default',
  style,
  ...props
}: React.HTMLAttributes<HTMLSpanElement> & { variant?: 'default' | 'success' | 'outline' }) {
  return (
    <span
      className={cn(
        'inline-flex items-center rounded-full px-3 py-1 text-xs font-semibold',
        variant === 'success'
          ? 'bg-emerald-500/15 text-emerald-300 ring-1 ring-emerald-500/30'
          : variant === 'outline'
          ? 'ring-1 ring-slate-600 text-slate-400'
          : 'bg-blue-500/15 text-blue-300 ring-1 ring-blue-500/30',
        className,
      )}
      style={style}
      {...props}
    />
  )
}
