import { useEffect, useRef } from 'react'
import gsap from 'gsap'
import type { DemoData } from '@/types/demo'
import { HAZE, color, project, shape } from '@/lib/demo'

const NS = 'http://www.w3.org/2000/svg'

/**
 * Full-screen intro. In the background the real Golub samples drift in from across the screen,
 * settle into the clusters evaluomeR finds, and an oval fades in around each cluster.
 */
export function Intro({ data }: { data: DemoData }) {
  const svgRef = useRef<SVGSVGElement>(null)

  useEffect(() => {
    const svg = svgRef.current!
    const still = matchMedia('(prefers-reduced-motion: reduce)').matches
    const k = data.optimalK
    const groups = data.clusters[k].map((c) => c - 1)
    let ctx: gsap.Context | undefined

    const run = () => {
      ctx?.revert()
      svg.replaceChildren()
      const W = svg.clientWidth, H = svg.clientHeight
      svg.setAttribute('viewBox', `0 0 ${W} ${H}`)
      // Desktop: right side, clear of the text. Smaller screens: whole section, faded behind the text.
      const wide = W >= 1024
      const R = wide ? Math.max(7, W / 170) : Math.max(6, W / 70)
      const at = project(data, wide ? [W * 0.56, W * 0.93, H * 0.2, H * 0.84] : [W * 0.06, W * 0.94, H * 0.1, H * 0.9])

      const ovals = Array.from({ length: k }, (_, g) => {
        const m = at.filter((_, i) => groups[i] === g)
        const n = m.length
        const cx = m.reduce((s, p) => s + p.x, 0) / n, cy = m.reduce((s, p) => s + p.y, 0) / n
        const a = m.reduce((s, p) => s + (p.x - cx) ** 2, 0) / n
        const b = m.reduce((s, p) => s + (p.x - cx) * (p.y - cy), 0) / n
        const c = m.reduce((s, p) => s + (p.y - cy) ** 2, 0) / n
        // Tilt and proportions from the cluster covariance, scaled so every sample sits inside.
        const t = (a + c) / 2, d = Math.sqrt(((a - c) / 2) ** 2 + b * b)
        const th = Math.atan2(t + d - a, b || 1e-9), cs = Math.cos(th), sn = Math.sin(th)
        const s1 = Math.sqrt(t + d) || 1, s2 = Math.sqrt(Math.max(t - d, 1e-6))
        const sc = Math.max(...m.map((p) => Math.hypot(((p.x - cx) * cs + (p.y - cy) * sn) / s1, (-(p.x - cx) * sn + (p.y - cy) * cs) / s2)))
        const e = document.createElementNS(NS, 'ellipse')
        Object.entries({ class: 'oval', cx, cy, rx: s1 * sc + R * 1.8, ry: s2 * sc + R * 1.8, transform: `rotate(${(th * 180) / Math.PI} ${cx} ${cy})`, stroke: color(g), fill: color(g) })
          .forEach(([k2, v]) => e.setAttribute(k2, String(v)))
        const w = document.createElementNS(NS, 'g')
        w.appendChild(e)
        svg.appendChild(w)
        return w
      })
      const dots = at.map((p, i) => {
        const e = document.createElementNS(NS, 'path')
        e.setAttribute('class', 'pt')
        e.setAttribute('d', shape(groups[i], R))
        e.setAttribute('fill', color(groups[i]))
        svg.appendChild(e)
        gsap.set(e, { x: p.x, y: p.y, transformOrigin: '50% 50%' })
        return e
      })
      if (still) return

      ctx = gsap.context(() => {
        gsap.set(ovals, { opacity: 0, scale: 0.94, transformOrigin: '50% 50%' })
        dots.forEach((e) => {
          e.setAttribute('fill', HAZE)
          gsap.set(e, { x: gsap.utils.random(0, W), y: gsap.utils.random(0, H), scale: 0.6 })
        })
        gsap.timeline()
          .to(dots, { x: (i) => at[i].x, y: (i) => at[i].y, scale: 1, duration: 2.2, ease: 'expo.inOut', stagger: { each: 0.012, from: 'random' } })
          .to(dots, { attr: { fill: (i: number) => color(groups[i]) }, duration: 0.6, stagger: { each: 0.008, from: 'random' } }, '-=1')
          .to(ovals, { opacity: 1, scale: 1, duration: 1.2, ease: 'sine.inOut', stagger: 0.15 }, 1.5)
          .add(() => dots.forEach((e, i) => gsap.to(e, {
            x: at[i].x + gsap.utils.random(-5, 5), y: at[i].y + gsap.utils.random(-5, 5),
            duration: gsap.utils.random(2, 4), ease: 'sine.inOut', repeat: -1, yoyo: true,
          })))
      }, svg)
    }

    run()
    let lastW = innerWidth
    const onResize = () => { if (Math.abs(innerWidth - lastW) > 80) { lastW = innerWidth; run() } }
    addEventListener('resize', onResize)
    return () => { removeEventListener('resize', onResize); ctx?.revert() }
  }, [data])

  return (
    <section className="intro" aria-labelledby="intro-h">
      <svg className="field" ref={svgRef} aria-hidden="true" />
      <div className="wrap inner">
        <p className="kind">R package</p>
        <h1 id="intro-h">evaluomeR</h1>
        <p className="lead">
          Discover the optimal number of clusters for your bioinformatics metrics, using bootstrap stability and
          silhouette quality on precomputed datasets.
        </p>
        <div className="ctas">
          <a className="btn primary" href="#playground">
            Open playground
            <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2.2" strokeLinecap="round" strokeLinejoin="round" aria-hidden="true"><path d="M12 5v14M5 12l7 7 7-7" /></svg>
          </a>
          <a className="btn ghost" href="https://neobernad.github.io/evaluomeR/">
            <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2" strokeLinecap="round" strokeLinejoin="round" aria-hidden="true"><path d="M2 4h7a3 3 0 0 1 3 3v13a2 2 0 0 0-2-2H2zM22 4h-7a3 3 0 0 0-3 3v13a2 2 0 0 1 2-2h8z" /></svg>
            Read docs
          </a>
        </div>
      </div>
      <p className="wrap note">Precomputed demo on the Golub and NCI-60 datasets.</p>
      <p className="cap" aria-hidden="true">{data.meta.nSamples} leukaemia samples, grouped into k = {data.optimalK}</p>
    </section>
  )
}
