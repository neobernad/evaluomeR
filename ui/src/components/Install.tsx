import { useState } from 'react'

const SNIPPETS = [
  'BiocManager::install("evaluomeR")',
  'stab <- stabilityRange(df, k.range = c(2, 6), bs = 100)\nqual <- qualityRange(df, k.range = c(2, 6))\ngetOptimalKValue(stab, qual)',
]

function Code({ text }: { text: string }) {
  const [copied, setCopied] = useState(false)
  const copy = () => navigator.clipboard.writeText(text).then(() => { setCopied(true); setTimeout(() => setCopied(false), 1500) })
  return (
    <div className="code">
      <pre>{text}</pre>
      <button className="copy" onClick={copy}>{copied ? 'Copied' : 'Copy'}</button>
    </div>
  )
}

export function Install() {
  return (
    <>
      <section className="run" id="run">
        <div className="wrap">
          <h2>Run it on your own metrics</h2>
          <div className="cols">
            <div>
              <p>evaluomeR takes a table of samples by metrics, resamples it, and scores every k for stability and separation.</p>
              <div className="links">
                <a href="https://neobernad.github.io/evaluomeR/articles/manual.html">Read the manual</a>
                <a href="https://github.com/neobernad/evaluomeR">GitHub</a>
              </div>
            </div>
            <div>{SNIPPETS.map((t) => <Code key={t} text={t} />)}</div>
          </div>
        </div>
      </section>
      <footer><div className="wrap">evaluomeR is released under GPL-3. Results were precomputed in R with 100 bootstraps per k.</div></footer>
    </>
  )
}
