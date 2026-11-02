import { useEffect, useLayoutEffect, useRef, useState, type KeyboardEvent, type RefObject } from 'react'
import { clusterGroups, color, kValues, oval, project, shape, type Dataset } from '@/lib/demo'
import { Logo } from '@/components/Logo'
import gsap from 'gsap'

const still = () => matchMedia('(prefers-reduced-motion: reduce)').matches
// What fades out and back in when the dataset changes.
const SWAP = '.head .ask, .head .sub, .strip, .stripcap, .body > *'

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

/**
 * Keeps a `.thumb` inside `ref` over the element matching `sel` (via --x/--y/--w/--h), so CSS can
 * slide it there. Re-measures on resize; the `ready` class (added after first paint) turns on the slide.
 */
function useThumb(ref: RefObject<HTMLElement | null>, sel: string, deps: unknown[]) {
  useLayoutEffect(() => {
    const el = ref.current!
    const place = () => {
      const b = el.querySelector(sel)
      if (!b) return
      const r = b.getBoundingClientRect(), c = el.getBoundingClientRect()
      Object.entries({ x: r.left - c.left, y: r.top - c.top, w: r.width, h: r.height }).forEach(([k, v]) => el.style.setProperty(`--${k}`, `${v}px`))
    }
    place()
    const raf = requestAnimationFrame(() => el.classList.add('ready'))
    const ro = new ResizeObserver(place)
    ro.observe(el)
    return () => { cancelAnimationFrame(raf); ro.disconnect() }
  }, deps)
}

/** Pick a k, see the verdict, the clusters on PC1/PC2, and the scores for every k. */
export function Explorer({ datasets }: { datasets: Dataset[] }) {
  const [key, setKey] = useState(() => (datasets.some((d) => d.key === readHash().get('dataset')) ? readHash().get('dataset')! : datasets[0].key))
  const [sel, setSel] = useState(key) // picker selection; moves at once, while `key` changes after the fade-out
  const ds = datasets.find((d) => d.key === key)!
  const D = ds.data
  const ks = kValues(D)
  const o = D.optimalK
  const [k, setK] = useState(() => { const h = Number(readHash().get('k')); return D.kSummary[h] ? h : o })
  const [byTruth, setByTruth] = useState(false)
  const stripRef = useRef<HTMLDivElement>(null)
  const barRef = useRef<HTMLDivElement>(null)
  const small = useMedia('(max-width: 639px)')
  const pillsRef = useRef<HTMLDivElement>(null)
  const rootRef = useRef<HTMLDivElement>(null)
  const svgRef = useRef<SVGSVGElement>(null)
  const tableRef = useRef<HTMLDivElement>(null)

  // Sliding highlights: dataset picker, k strip, selected table row.
  useThumb(pillsRef, '[aria-pressed=true]', [sel])
  useThumb(stripRef, '[aria-checked=true]', [key, k])
  useThumb(tableRef, 'tr.sel', [key, k])

  useEffect(() => { history.replaceState(null, '', `#dataset=${key}&k=${k}`) }, [key, k])

  // Shadow under the k strip only while it is stuck (phone/tablet).
  useEffect(() => {
    const io = new IntersectionObserver(([e]) => barRef.current?.classList.toggle('stuck', e.intersectionRatio < 1), { threshold: 1, rootMargin: '-1px 0px 0px 0px' })
    io.observe(barRef.current!)
    return () => io.disconnect()
  }, [])

  // Dataset switch: fade the old results out, swap, fade the new ones in (layout effect below).
  const pickDataset = (d: Dataset) => {
    if (d.key === sel) return
    setSel(d.key)
    const swap = () => { setKey(d.key); setK(d.data.optimalK) }
    if (still()) return swap()
    const els = rootRef.current!.querySelectorAll(SWAP)
    gsap.killTweensOf(els)
    gsap.to(els, { autoAlpha: 0, y: -8, duration: 0.18, ease: 'power1.in', stagger: 0.02, onComplete: swap })
  }
  const firstKey = useRef(true)
  useLayoutEffect(() => {
    if (firstKey.current) { firstKey.current = false; return }
    const els = rootRef.current!.querySelectorAll(SWAP)
    gsap.fromTo(els, { autoAlpha: 0, y: 12 }, { autoAlpha: 1, y: 0, duration: 0.5, ease: 'power2.out', stagger: 0.06, clearProps: 'transform,opacity,visibility' })
  }, [key])

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
  const view = { fit, mx: (bx0 + bx1) / 2, my: (by0 + by1) / 2 }
  const sceneOf = (v: typeof view) => `translate(${W / 2} ${H / 2}) scale(${v.fit}) translate(${-v.mx} ${-v.my})`
  const scene = sceneOf(view)

  // Points whose colour/shape changes (view toggle or new k) pop in. On a new k each cluster oval grows
  // out of the old oval most of its samples came from, and the zoom tweens instead of jumping.
  // All tweens drive plain numbers and write attributes (no GSAP transforms), so nothing is left over when they end.
  const shown = at.map((_, i) => group(i))
  const prev = useRef<{ key: string; k: number; shown: number[]; ovals: typeof ovals; view: typeof view } | null>(null)
  const tweens = useRef<gsap.core.Tween[]>([])
  const pops = useRef<gsap.core.Tween[]>([])
  useLayoutEffect(() => {
    const p = prev.current, svg = svgRef.current!
    prev.current = { key, k, shown, ovals, view }
    if (!p || p.key !== key || still()) return
    // Pop by redrawing the shape at a growing size, not by a GSAP transform: GSAP's SVG transform origin
    // put points ~7px off their centre mid-pop, so they jumped when it ended.
    // A click mid-pop: stop the old pops and redraw every point at full size first.
    const r = rad / fit
    pops.current.forEach((t) => t.kill())
    pops.current = []
    svg.querySelectorAll<SVGPathElement>('.pt path').forEach((e, i) => {
      const o = { s: p.shown[i] === shown[i] ? 1 : 0.2 }, draw = () => e.setAttribute('d', shape(shown[i], r * o.s))
      draw()
      if (o.s < 1) pops.current.push(gsap.to(o, { s: 1, duration: 0.55, ease: 'back.out(2.2)', delay: Math.random() * 0.15, onUpdate: draw }))
    })
    if (p.k === k) return

    tweens.current.forEach((t) => t.kill())
    const ease = 'power3.inOut', duration = 0.7
    const g = svg.firstElementChild as SVGGElement
    const v = { ...p.view }
    g.setAttribute('transform', sceneOf(v)) // start from the old zoom on this frame
    tweens.current = [gsap.to(v, { ...view, duration, ease, onUpdate: () => g.setAttribute('transform', sceneOf(v)) })]
    const prevCl = D.clusters[p.k], cl = D.clusters[k]
    svg.querySelectorAll<SVGEllipseElement>('.oval').forEach((e, c) => {
      // Parent: the old cluster holding most of this cluster's samples.
      const n: Record<number, number> = {}
      cl.forEach((x, i) => { if (x === c + 1) n[prevCl[i]] = (n[prevCl[i]] ?? 0) + 1 })
      const parent = +Object.entries(n).sort((a, b) => b[1] - a[1])[0][0]
      const to = ovals[c], from = { ...p.ovals[parent - 1] }
      from.deg += Math.round((to.deg - from.deg) / 180) * 180 // an ellipse turned 180° is the same: take the short way
      const o = { cx: from.cx, cy: from.cy, rx: from.rx, ry: from.ry, deg: from.deg }
      const set = () => {
        Object.entries({ cx: o.cx, cy: o.cy, rx: o.rx, ry: o.ry }).forEach(([a, x]) => e.setAttribute(a, String(x)))
        e.setAttribute('transform', `rotate(${o.deg} ${o.cx} ${o.cy})`)
      }
      set()
      tweens.current.push(gsap.to(o, { cx: to.cx, cy: to.cy, rx: to.rx, ry: to.ry, deg: to.deg, duration, ease, onUpdate: set }))
    })
    gsap.fromTo(rootRef.current!.querySelector('.v'), { opacity: 0.25 }, { opacity: 1, duration: 0.4, ease: 'sine.out', clearProps: 'opacity' })
  }, [key, k, byTruth])

  const pct = (x: number) => `${Math.max(0, Math.min(1, x)) * 100}%`
  const cell = (x: number | undefined) => (
    <td><span className="num">{(x ?? NaN).toFixed(2)}</span><span className="meter" aria-hidden="true"><i style={{ width: pct(x ?? 0) }} /></span></td>
  )
  const [v1, v2] = D.pca.varExplained

  return (
      <div className="wrap" id="playground" ref={rootRef}>
        <header className="top">
          <span className="brand"><Logo />evaluomeR</span>
          <div className="pick">
            <span className="lbl" id="ds-lbl"><span className="long">Example data</span><span className="short">Data</span></span>
            <div className="pills" ref={pillsRef} role="group" aria-labelledby="ds-lbl">
              <span className="thumb" aria-hidden="true" />
              {datasets.map((d) => (
                <button key={d.key} aria-pressed={d.key === sel} onClick={() => pickDataset(d)}>
                  <b><span className="short">{d.short}</span><span className="long">{d.label}</span></b>
                  <small>{d.data.meta.nSamples} samples, {d.data.meta.cancerTypes.length} types</small>
                </button>
              ))}
            </div>
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
              <span className="thumb" aria-hidden="true" />
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
              <h3>{byTruth ? `Known labels, ovals show the ${k} clusters found` : `The ${k} clusters found`}</h3>
              <div className="seg" role="group" aria-label="Colour points by">
                <button aria-pressed={!byTruth} onClick={() => setByTruth(false)}>Clusters</button>
                <button aria-pressed={byTruth} onClick={() => setByTruth(true)}>Known labels</button>
              </div>
            </div>
            <svg className="sc" ref={svgRef} viewBox={`0 0 ${W} ${H}`} role="img" aria-describedby="scdesc">
              <g transform={scene}>
              {ovals.map(({ g, cx, cy, rx, ry, deg }) => (
                <ellipse key={g} className="oval" cx={cx} cy={cy} rx={rx} ry={ry} transform={`rotate(${deg} ${cx} ${cy})`} stroke={color(g)} fill={color(g)} />
              ))}
              {at.map((p, i) => (
                <g key={i} className="pt" transform={`translate(${p.x},${p.y})`}><path d={shape(group(i), rad / fit)} fill={color(group(i))} stroke="rgba(11,22,64,.65)" strokeWidth={1.2 / fit}>
                  <title>{`${D.samples[i].id}: ${D.samples[i].cancerType}, cluster ${D.clusters[k][i]}`}</title>
                </path></g>
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
          <div className="t" ref={tableRef}>
            <span className="thumb" aria-hidden="true" />
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
