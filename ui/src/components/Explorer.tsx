import { useEffect, useRef, useState, type KeyboardEvent } from 'react'
import { clusterGroups, color, kValues, oval, project, shape, type Dataset } from '@/lib/demo'
import { Logo } from '@/components/Logo'

const readHash = () => new URLSearchParams(location.hash.slice(1))
const useMedia = (q: string) => {
  const [on, setOn] = useState(() => matchMedia(q).matches)
  useEffect(() => {
    const m = matchMedia(q), f = () => setOn(m.matches)
    m.addEventListener('change', f)
    return () => m.removeEventListener('change', f)
  }, [q])
  return on
}

/** Pick a k, see the verdict, the clusters on PC1/PC2, and the scores for every k. */
export function Explorer({ datasets }: { datasets: Dataset[] }) {
  const [key, setKey] = useState(() => (datasets.some((d) => d.key === readHash().get('dataset')) ? readHash().get('dataset')! : datasets[0].key))
  const ds = datasets.find((d) => d.key === key)!
  const D = ds.data
  const ks = kValues(D)
  const o = D.optimalK
  const [k, setK] = useState(() => { const h = Number(readHash().get('k')); return D.kSummary[h] ? h : o })
  const [byTruth, setByTruth] = useState(false)
  const stripRef = useRef<HTMLDivElement>(null)
  const barRef = useRef<HTMLDivElement>(null)
  const small = useMedia('(max-width: 639px)')

  useEffect(() => { history.replaceState(null, '', `#dataset=${key}&k=${k}`) }, [key, k])

  // Shadow under the k strip only while it is stuck (phone/tablet).
  useEffect(() => {
    const io = new IntersectionObserver(([e]) => barRef.current?.classList.toggle('stuck', e.intersectionRatio < 1), { threshold: 1, rootMargin: '-1px 0px 0px 0px' })
    io.observe(barRef.current!)
    return () => io.disconnect()
  }, [])

  const pickDataset = (d: Dataset) => { setKey(d.key); setK(d.data.optimalK) }

  // Roving-tabindex radiogroup: arrows move and select, Home/End jump.
  const onKey = (e: KeyboardEvent) => {
    const i = ks.indexOf(k), m = ({ ArrowRight: 1, ArrowDown: 1, ArrowLeft: -1, ArrowUp: -1 } as Record<string, number>)[e.key]
    let n: number | null = null
    if (m) n = ks[(i + m + ks.length) % ks.length]
    if (e.key === 'Home') n = ks[0]
    if (e.key === 'End') n = ks[ks.length - 1]
    if (n === null) return
    e.preventDefault()
    setK(n)
    stripRef.current?.querySelector<HTMLElement>(`[data-k="${n}"]`)?.focus()
  }

  const S = D.kSummary, s = S[k]
  const best = Math.max(...ks.map((v) => S[v].composite))
  const ari = (v: number) => (S[v].ari ?? NaN).toFixed(2)
  const verdict = k === o
    ? <>{k} groups. This is the answer <span className="pkg">evaluomeR</span> gives.</>
    : `${k} groups works, but less well than ${o}.`
  const why = k === o
    ? `Most stable under resampling and best separated of k ${ks[0]}–${ks[ks.length - 1]}${s.ari != null ? `, and it matches the known labels with an adjusted Rand index of ${ari(k)}` : ''}.`
    : `Its combined score is ${((s.composite / S[o].composite) * 100).toFixed(0)}% of k = ${o}. Switch to known labels to see which groups get split or merged.`

  // Scatter
  const W = 640, H = small ? 640 : 400, rad = small ? 10 : 7, pad = small ? 18 : 22
  const at = project(D, [pad, W - pad, pad, H - pad])
  const types = D.meta.cancerTypes
  const cmap = clusterGroups(D, k) // cluster c is drawn like the label it mostly holds
  const group = (i: number) => (byTruth ? types.indexOf(D.samples[i].cancerType) : cmap[D.clusters[k][i] - 1])
  const legend = byTruth ? types.map((t, i) => [t, i] as const) : cmap.map((g, c) => [`Cluster ${c + 1}`, g] as const)
  const cg = (i: number) => cmap[D.clusters[k][i] - 1]
  const ovals = cmap.map((g) => ({ g, ...oval(at.filter((_, i) => cg(i) === g), rad * 1.8) }))
  // Scale and centre the scene so the cluster ovals fill the plot (same in both views, so points don't jump).
  const box = ovals.map(({ cx, cy, rx, ry, deg }) => {
    const c = Math.cos((deg * Math.PI) / 180), sn = Math.sin((deg * Math.PI) / 180)
    const hx = Math.hypot(rx * c, ry * sn), hy = Math.hypot(rx * sn, ry * c)
    return [cx - hx, cx + hx, cy - hy, cy + hy]
  })
  const [bx0, bx1, by0, by1] = [Math.min(...box.map((v) => v[0])), Math.max(...box.map((v) => v[1])), Math.min(...box.map((v) => v[2])), Math.max(...box.map((v) => v[3]))]
  const fit = Math.min((W - 4) / (bx1 - bx0), (H - 4) / (by1 - by0))
  const scene = `translate(${W / 2} ${H / 2}) scale(${fit}) translate(${-(bx0 + bx1) / 2} ${-(by0 + by1) / 2})`
  const pct = (x: number) => `${Math.max(0, Math.min(1, x)) * 100}%`
  const cell = (x: number | undefined) => (
    <td><span className="num">{(x ?? NaN).toFixed(2)}</span><span className="meter" aria-hidden="true"><i style={{ width: pct(x ?? 0) }} /></span></td>
  )
  const [v1, v2] = D.pca.varExplained

  return (
      <div className="wrap" id="playground">
        <header className="top">
          <span className="brand"><Logo />evaluomeR</span>
          <div className="pills" role="group" aria-label="Dataset">
            {datasets.map((d) => (
              <button key={d.key} aria-pressed={d.key === key} onClick={() => pickDataset(d)}>
                <span className="short">{d.short}</span><span className="long">{d.label}</span>
              </button>
            ))}
          </div>
          <nav aria-label="Project">
            <a className="hide-sm" href="https://neobernad.github.io/evaluomeR/">Documentation</a>
            <a className="hide-sm" href="https://github.com/neobernad/evaluomeR">GitHub</a>
            <a href="#run">Install</a>
          </nav>
        </header>

        <section className="head">
          <div>
            <h2 className="ask">How many groups are in these {D.meta.nSamples} {ds.noun}?</h2>
            <p className="sub">Pick a k. Each was scored on {D.meta.nMetrics} metrics with {D.meta.bs} bootstraps.</p>
          </div>
          <div className="stripbar" ref={barRef}>
            <div className="strip" ref={stripRef} role="radiogroup" aria-label="Number of clusters k" onKeyDown={onKey}>
              {ks.map((v) => (
                <button key={v} role="radio" data-k={v} aria-checked={v === k} tabIndex={v === k ? 0 : -1}
                  aria-label={`k = ${v}${v === o ? ', recommended' : ''}`} onClick={() => setK(v)}>
                  <span className="nwrap"><span className="n" aria-hidden="true">{v}</span>{v === o && <span className="best" aria-hidden="true" />}</span>
                  <span className="bar" aria-hidden="true"><i style={{ width: `${(S[v].composite / best) * 100}%` }} /></span>
                </button>
              ))}
            </div>
            <p className="stripcap"><span>Bar: combined score</span><span className="dot">Recommended</span></p>
          </div>
        </section>

        <main className="body">
          <div className="v" aria-live="polite" aria-atomic="true">
            <p className="verdict">{verdict}</p>
            <p className="why">{why}</p>
          </div>
          <div className="plotbox">
            <div className="plothead">
              <h3>{byTruth ? 'Coloured by known label' : `The ${k} clusters found`}</h3>
              <div className="seg" role="group" aria-label="Colour points by">
                <button aria-pressed={!byTruth} onClick={() => setByTruth(false)}>Clusters</button>
                <button aria-pressed={byTruth} onClick={() => setByTruth(true)}>Known labels</button>
              </div>
            </div>
            <svg className="sc" viewBox={`0 0 ${W} ${H}`} role="img" aria-describedby="scdesc">
              <g transform={scene}>
              {!byTruth && ovals.map(({ g, cx, cy, rx, ry, deg }) => (
                <ellipse key={g} className="oval" cx={cx} cy={cy} rx={rx} ry={ry} transform={`rotate(${deg} ${cx} ${cy})`} stroke={color(g)} fill={color(g)} />
              ))}
              {at.map((p, i) => (
                <path key={i} d={shape(group(i), rad / fit)} fill={color(group(i))} stroke="rgba(11,22,64,.65)" strokeWidth={1.2 / fit} transform={`translate(${p.x},${p.y})`}>
                  <title>{`${D.samples[i].id}: ${D.samples[i].cancerType}, cluster ${D.clusters[k][i]}`}</title>
                </path>
              ))}
              </g>
            </svg>
            <p className="sr" id="scdesc">Scatter of samples on the first two principal components. The table lists the same scores for every k.</p>
            <div className="legend">
              {legend.map(([t, i]) => (
                <span key={t}><svg viewBox="-8 -8 16 16" aria-hidden="true"><path d={shape(i, 6)} fill={color(i)} /></svg>{t}</span>
              ))}
            </div>
            <p className="axes">PC1 explains {(v1 * 100).toFixed(0)}% of variance, PC2 {(v2 * 100).toFixed(0)}%.</p>
          </div>
          <div className="t">
            <table>
              <caption>Scores for every k</caption>
              <thead><tr>
                <th scope="col">k</th><th scope="col">Stability</th><th scope="col">Silhouette</th>
                <th scope="col"><span className="lbl-long">Matches labels</span><span className="lbl-short">Labels</span></th>
              </tr></thead>
              <tbody>
                {ks.map((v) => (
                  <tr key={v} className={v === k ? 'sel' : undefined}>
                    <th scope="row">{v}{v === o && <><span className="best" aria-hidden="true" /><span className="sr"> (recommended)</span></>}</th>
                    {cell(S[v].avgStability)}{cell(S[v].avgSilhouette)}{cell(S[v].ari)}
                  </tr>
                ))}
              </tbody>
            </table>
          </div>
        </main>
      </div>
  )
}
