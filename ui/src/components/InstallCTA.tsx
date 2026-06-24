import { useCallback, useEffect, useRef, useState } from 'react'
import { Check, Copy } from 'lucide-react'
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from '@/components/ui/card'

const BIOC_CODE = `if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("evaluomeR")`

const GITHUB_CODE = `devtools::install_github("neobernad/evaluomeR")`

type CopyKey = 'bioc' | 'github'

function CodeBlock({ label, code, copyKey, copied, onCopy }: {
  label: string
  code: string
  copyKey: CopyKey
  copied: CopyKey | null
  onCopy: (key: CopyKey, text: string) => void
}) {
  const isCopied = copied === copyKey

  return (
    <div>
      <p className="mb-2 text-xs font-medium uppercase tracking-wide text-slate-500">{label}</p>
      <div className="relative">
        <button
          type="button"
          onClick={() => onCopy(copyKey, code)}
          className="absolute right-2 top-2 rounded-md border border-slate-700 bg-slate-900/80 p-1.5 text-slate-400 transition hover:border-slate-600 hover:text-white"
          aria-label={isCopied ? 'Copied' : `Copy ${label} install command`}
        >
          {isCopied ? <Check className="h-3.5 w-3.5 text-emerald-400" /> : <Copy className="h-3.5 w-3.5" />}
        </button>
        <pre className="overflow-x-auto rounded-lg border border-slate-800 bg-slate-950 p-4 pr-12 font-mono text-sm text-slate-200">
          {code}
        </pre>
      </div>
    </div>
  )
}

export function InstallCTA() {
  const [copied, setCopied] = useState<CopyKey | null>(null)
  const timerRef = useRef<ReturnType<typeof setTimeout> | null>(null)

  useEffect(() => {
    return () => {
      if (timerRef.current) clearTimeout(timerRef.current)
    }
  }, [])

  const handleCopy = useCallback((key: CopyKey, text: string) => {
    void navigator.clipboard.writeText(text).then(() => {
      setCopied(key)
      if (timerRef.current) clearTimeout(timerRef.current)
      timerRef.current = setTimeout(() => setCopied(null), 2000)
    })
  }, [])

  return (
    <section className="mx-auto max-w-5xl px-6 pb-24">
      <Card className="border-blue-500/20 bg-gradient-to-br from-slate-900 to-blue-950/40">
        <CardHeader>
          <CardTitle>Install evaluomeR</CardTitle>
          <CardDescription>
            Run stability, quality, and optimal-k analyses on your own metric datasets —
            including cancer cell lines, ontology metrics, and more.
          </CardDescription>
        </CardHeader>
        <CardContent className="space-y-4">
          <CodeBlock
            label="Bioconductor"
            code={BIOC_CODE}
            copyKey="bioc"
            copied={copied}
            onCopy={handleCopy}
          />
          <CodeBlock
            label="GitHub (development)"
            code={GITHUB_CODE}
            copyKey="github"
            copied={copied}
            onCopy={handleCopy}
          />
          <p className="text-sm text-slate-400">
            <a
              href="https://neobernad.github.io/evaluomeR/"
              className="text-blue-400 underline-offset-4 hover:underline"
            >
              View full documentation →
            </a>
          </p>
        </CardContent>
      </Card>
    </section>
  )
}
