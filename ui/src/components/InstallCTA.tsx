import { Card, CardContent, CardDescription, CardHeader, CardTitle } from '@/components/ui/card'

export function InstallCTA() {
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
          <div>
            <p className="mb-2 text-xs font-medium uppercase tracking-wide text-slate-500">
              Bioconductor
            </p>
            <pre className="overflow-x-auto rounded-lg border border-slate-800 bg-slate-950 p-4 font-mono text-sm text-slate-200">
{`if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("evaluomeR")`}
            </pre>
          </div>
          <div>
            <p className="mb-2 text-xs font-medium uppercase tracking-wide text-slate-500">
              GitHub (development)
            </p>
            <pre className="overflow-x-auto rounded-lg border border-slate-800 bg-slate-950 p-4 font-mono text-sm text-slate-200">
{`devtools::install_github("neobernad/evaluomeR")`}
            </pre>
          </div>
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
