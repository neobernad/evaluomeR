interface SectionDividerProps {
  label: string
}

export function SectionDivider({ label }: SectionDividerProps) {
  return (
    <div
      className="relative py-8"
      role="separator"
      aria-label={label}
    >
      <div
        aria-hidden
        className="pointer-events-none absolute inset-x-0 top-1/2 h-px -translate-y-1/2 bg-gradient-to-r from-transparent via-slate-600/80 to-transparent"
      />
      <div
        aria-hidden
        className="pointer-events-none absolute inset-x-12 top-1/2 h-px -translate-y-1/2 bg-gradient-to-r from-transparent via-blue-500/25 to-transparent blur-[0.5px]"
      />
      <div className="relative flex items-center justify-center gap-3">
        <span
          aria-hidden
          className="h-px w-8 bg-gradient-to-r from-transparent to-slate-600/60 sm:w-12"
        />
        <span
          aria-hidden
          className="h-1.5 w-1.5 rotate-45 rounded-[1px] bg-gradient-to-br from-blue-400/80 to-teal-400/60 shadow-[0_0_8px_rgba(59,130,246,0.35)]"
        />
        <span className="bg-slate-900/90 px-3 text-[10px] font-semibold uppercase tracking-[0.22em] text-slate-500">
          {label}
        </span>
        <span
          aria-hidden
          className="h-1.5 w-1.5 rotate-45 rounded-[1px] bg-gradient-to-br from-blue-400/80 to-teal-400/60 shadow-[0_0_8px_rgba(59,130,246,0.35)]"
        />
        <span
          aria-hidden
          className="h-px w-8 bg-gradient-to-l from-transparent to-slate-600/60 sm:w-12"
        />
      </div>
    </div>
  )
}
