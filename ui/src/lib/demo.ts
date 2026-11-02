import type { DemoData } from '@/types/demo'

// Cluster colours, readable on cobalt. Each group also gets a shape so colour is never the only cue.
export const COLORS = ['#FFD84D', '#FF8FA3', '#7CE0C3', '#C9B6FF', '#FFB067', '#8FD3FF', '#F2F2F2', '#B5E36B', '#FF6B6B', '#9AA6FF', '#E0A3FF', '#66D9E8']
export const HAZE = '#C3CEF8'

const SHAPES = [
  (r: number) => `M${-r},0a${r},${r} 0 1,0 ${2 * r},0a${r},${r} 0 1,0 ${-2 * r},0`,
  (r: number) => `M${-r * 0.9},${-r * 0.9}h${r * 1.8}v${r * 1.8}h${-r * 1.8}z`,
  (r: number) => `M0,${-r * 1.1}L${r * 1.05},${r * 0.8}H${-r * 1.05}z`,
  (r: number) => `M0,${-r * 1.2}L${r * 1.1},0L0,${r * 1.2}L${-r * 1.1},0z`,
  (r: number) => `M0,${r * 1.1}L${r * 1.05},${-r * 0.8}H${-r * 1.05}z`,
  (r: number) => `M${-r * 0.4},${-r}h${r * 0.8}v${r * 0.6}h${r * 0.6}v${r * 0.8}h${-r * 0.6}v${r * 0.6}h${-r * 0.8}v${-r * 0.6}h${-r * 0.6}v${-r * 0.8}h${r * 0.6}z`,
]
export const shape = (group: number, r: number) => SHAPES[group % SHAPES.length](r)
export const color = (group: number) => COLORS[group % COLORS.length]

export interface Dataset {
  key: string
  label: string
  short: string
  noun: string
  data: DemoData
}

/**
 * Colour/shape index for each sample's cluster at this k. Each cluster takes the index of the
 * known label it holds most of (largest overlaps first, one cluster per label), so a cluster
 * looks the same in both views. Leftover clusters take the lowest unused index.
 */
export function clusterGroups(d: DemoData, k: number) {
  const types = d.meta.cancerTypes, cl = d.clusters[k]
  const count = Array.from({ length: k }, () => types.map(() => 0))
  cl.forEach((c, i) => count[c - 1][types.indexOf(d.samples[i].cancerType)]++)
  const pairs = count.flatMap((row, c) => row.map((v, t) => [v, c, t])).sort((a, b) => b[0] - a[0])
  const map: number[] = Array(k).fill(-1), used = new Set<number>()
  for (const [v, c, t] of pairs) if (v > 0 && map[c] < 0 && !used.has(t)) { map[c] = t; used.add(t) }
  for (let c = 0, t = 0; c < k; c++) if (map[c] < 0) { while (used.has(t)) t++; map[c] = t; used.add(t) }
  return map
}

/** Ellipse around points: tilt and proportions from their covariance, grown so every point sits inside, plus `pad`. */
export function oval(m: { x: number; y: number }[], pad: number) {
  const n = m.length
  const cx = m.reduce((s, p) => s + p.x, 0) / n, cy = m.reduce((s, p) => s + p.y, 0) / n
  const a = m.reduce((s, p) => s + (p.x - cx) ** 2, 0) / n
  const b = m.reduce((s, p) => s + (p.x - cx) * (p.y - cy), 0) / n
  const c = m.reduce((s, p) => s + (p.y - cy) ** 2, 0) / n
  const t = (a + c) / 2, d = Math.sqrt(((a - c) / 2) ** 2 + b * b)
  const th = Math.atan2(t + d - a, b || 1e-9), cs = Math.cos(th), sn = Math.sin(th)
  const s1 = Math.sqrt(t + d) || 1, s2 = Math.sqrt(Math.max(t - d, 1e-6))
  const sc = Math.max(...m.map((p) => Math.hypot(((p.x - cx) * cs + (p.y - cy) * sn) / s1, (-(p.x - cx) * sn + (p.y - cy) * cs) / s2)))
  return { cx, cy, rx: s1 * sc + pad, ry: s2 * sc + pad, deg: (th * 180) / Math.PI }
}

export const kValues =(d: DemoData) => Object.keys(d.kSummary).map(Number).sort((a, b) => a - b)

/** Maps PC1/PC2 into a [x0,x1]×[y0,y1] box (y up). */
export function project(d: DemoData, box: [number, number, number, number]) {
  const xs = d.pca.coords.map((c) => c.pc1)
  const ys = d.pca.coords.map((c) => c.pc2)
  const [mx, Mx, my, My] = [Math.min(...xs), Math.max(...xs), Math.min(...ys), Math.max(...ys)]
  const [x0, x1, y0, y1] = box
  return d.pca.coords.map((c) => ({
    x: x0 + ((c.pc1 - mx) / (Mx - mx)) * (x1 - x0),
    y: y1 - ((c.pc2 - my) / (My - my)) * (y1 - y0),
  }))
}
